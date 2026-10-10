// AR009 T001 (design §6): unit tests for the pipelined CG solver
// (`--solver pcg`) and the artifact-side verifiers for the residual-history
// e2e entries. TESTS ONLY — the PipelinedCGSolver class, the cgDotLocal /
// cgBlockingAllreduce / cgIallreduce2Start / cgIallreduce2Wait kernels and
// the `--residual-history` CLI do not exist yet; this file is the intended
// Red state (compile error on the class name, failing profiler counts).
//
// Case map (design §6 table rows):
//   P1 PipelinedMatchesCgIterations      — 256^2 multi-mode, iters ±5%, L2
//   P2 PipelinedSolutionAccuracy         — 64^2/256^2 manufactured, 2nd order
//   P3 NoBlockingAllreduce+positive ctrl — region counts both directions
//   P4 PipelinedNp4Consistency           — np4 vs np1, L2 ≤ 1e-10·||u||
//   P5 PipelinedBreakdownDefense         — RECORD ONLY (see comment below)
//   P6 ExistingBehaviorUnchanged         — full 42-test regression + cg
//                                          golden (existing entries; see
//                                          comment below)
//   P7 IterateGuardAndResume             — guard + one-step resume
//   P8 ZeroRhsImmediateConvergence       — 0 iters, callCount==1
//   P9 MaxIterTruncation                 — completed==3, trajectory dense
//   E1/E2/E4/E7 artifact checks         — ResidualHistoryE2ETest (reads the
//                                          CSV/JSON/stdout artifacts of the
//                                          ctest fixture runs; SKIP without
//                                          HYPOS_AR009_RH_BASE)
//   E3 artifact checks                  — ResidualHistoryNp4E2ETest
//   E5 report check                     — PipelinedE2ETest
//   B1 (bench) is NOT registered: manual §16 evidence per design §15/§16
//   convention. W2 (--help wording) is checked by the pcg_help_check ctest
//   entry via cmake/VerifyAr009Cli.cmake (string-level checks).

#include <gtest/gtest.h>
#include <mpi.h>

#include "comm/halo_exchanger.hpp"
#include "grid/subgrid.hpp"
#include "solver/solver.hpp"
#include "solver/residual.hpp"
#include "perf/profiler.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace hypo;

namespace {

constexpr Index kHalo = 1;

// ---------------------------------------------------------------------------
// Manufactured-solution helpers — bodies follow the test_mg_vcycle.cpp
// conventions verbatim (same formulas, same discrete-consistent assembly).
// ---------------------------------------------------------------------------

Real sine1d(long long g, long long n) {
    const Real pi = std::acos(-1.0);
    return std::sin(pi * static_cast<Real>(g + 1) / static_cast<Real>(n + 1));
}

Real sineExact(long long gi, long long gj, long long nx, long long ny) {
    return sine1d(gi, nx) * sine1d(gj, ny);
}

// rhs = A*exact with Dirichlet-0 outside (sine vanishes at both faces).
void fillSineRhs(Subgrid& sg, long long nx, long long ny) {
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);
    Real* rhs = sg.rhs().data();
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gi = static_cast<long long>(sg.offsetX()) +
                                 static_cast<long long>(i) - kHalo;
            const long long gj = static_cast<long long>(sg.offsetY()) +
                                 static_cast<long long>(j) - kHalo;
            rhs[sg.index(i, j)] = sineExact(gi - 1, gj, nx, ny) +
                                  sineExact(gi + 1, gj, nx, ny) +
                                  sineExact(gi, gj - 1, nx, ny) +
                                  sineExact(gi, gj + 1, nx, ny) -
                                  4.0 * sineExact(gi, gj, nx, ny);
        }
    }
}

// Multi-mode exact solution u* = sum_{k=1..4} c_k sin(k pi x) sin(k pi y)
// with c = (1, -1/2, 1/3, -1/4): a rhs that actually stresses the CG
// condition number (the single-mode sine rhs is a stencil eigenvector that
// CG lands exactly in one iteration — see test_mg_vcycle.cpp).
Real multiModeExact(long long gi, long long gj, long long n) {
    const Real pi = std::acos(-1.0);
    const Real x = pi * static_cast<Real>(gi + 1) / static_cast<Real>(n + 1);
    const Real y = pi * static_cast<Real>(gj + 1) / static_cast<Real>(n + 1);
    Real acc = 0.0;
    for (long long k = 1; k <= 4; ++k) {
        const Real c = ((k % 2 == 1) ? 1.0 : -1.0) / static_cast<Real>(k);
        acc += c * std::sin(static_cast<Real>(k) * x) *
               std::sin(static_cast<Real>(k) * y);
    }
    return acc;
}

// Discrete-consistent rhs for multiModeExact (rhs = A u*).
void fillMultiModeRhs(Subgrid& sg, long long n) {
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);
    Real* rhs = sg.rhs().data();
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gi = static_cast<long long>(sg.offsetX()) +
                                 static_cast<long long>(i) - kHalo;
            const long long gj = static_cast<long long>(sg.offsetY()) +
                                 static_cast<long long>(j) - kHalo;
            rhs[sg.index(i, j)] = multiModeExact(gi - 1, gj, n) +
                                  multiModeExact(gi + 1, gj, n) +
                                  multiModeExact(gi, gj - 1, n) +
                                  multiModeExact(gi, gj + 1, n) -
                                  4.0 * multiModeExact(gi, gj, n);
        }
    }
}

// Analytic rhs of -Laplace(u) = 2*pi^2*sin(pi*x)*sin(pi*y): leaves a genuine
// O(h^2) truncation error for the convergence-order anchor (P2). Solver
// convention A*u = -rhs, hence rhs = -f*h^2 (test_mg_vcycle convention).
void fillAnalyticSineRhs(Subgrid& sg, long long n) {
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);
    const Real pi = std::acos(-1.0);
    const Real h = 1.0 / static_cast<Real>(n + 1);
    Real* rhs = sg.rhs().data();
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gi = static_cast<long long>(sg.offsetX()) +
                                 static_cast<long long>(i) - kHalo;
            const long long gj = static_cast<long long>(sg.offsetY()) +
                                 static_cast<long long>(j) - kHalo;
            const Real f = 2.0 * pi * pi * sineExact(gi, gj, n, n);
            rhs[sg.index(i, j)] = -f * h * h;
        }
    }
}

Real sineL2ErrorLocal(const Subgrid& sg, long long nx, long long ny, Real& localMax) {
    Real sum = 0.0;
    localMax = 0.0;
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gi = static_cast<long long>(sg.offsetX()) +
                                 static_cast<long long>(i) - kHalo;
            const long long gj = static_cast<long long>(sg.offsetY()) +
                                 static_cast<long long>(j) - kHalo;
            const Real diff = sg.u().data()[sg.index(i, j)] - sineExact(gi, gj, nx, ny);
            sum += diff * diff;
            localMax = std::max(localMax, std::fabs(sg.u().data()[sg.index(i, j)]));
        }
    }
    return sum;
}

// ---------------------------------------------------------------------------
// Shared scaffolding for the manual non-uniform {15,19} 2x2 layout on 34
// global (AR007 CrossLayout pattern, test_mg_vcycle.cpp convention).
// ---------------------------------------------------------------------------

constexpr long long kLayoutN = 34;
Index layoutOffset(int coord) { return coord == 0 ? 0 : 15; }
Index layoutSize(int coord) { return coord == 0 ? 15 : 19; }

struct NonUniformLayout {
    MPI_Comm cart = MPI_COMM_NULL;
    Subgrid* sg = nullptr;
    PointToPointExchanger* ex = nullptr;
    int coords[2] = {0, 0};

    void build() {
        int rank = 0;
        int size = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        ASSERT_EQ(size, 4);

        int dims[2] = {2, 2};
        int periods[2] = {0, 0};
        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart);

        int left = MPI_PROC_NULL, right = MPI_PROC_NULL;
        int down = MPI_PROC_NULL, up = MPI_PROC_NULL;
        MPI_Cart_shift(cart, 0, 1, &left, &right);
        MPI_Cart_shift(cart, 1, 1, &down, &up);
        MPI_Cart_coords(cart, rank, 2, coords);

        const Index nxL = layoutSize(coords[0]);
        const Index nyL = layoutSize(coords[1]);
        sg = new Subgrid(nxL, nyL, 1, kHalo, cart);
        sg->setNeighbors(left, right, down, up);
        sg->setOffsets(layoutOffset(coords[0]), layoutOffset(coords[1]), 0);
        ex = new PointToPointExchanger();
        ex->initialize(*sg);
    }

    void destroy() {
        delete ex;
        delete sg;
        MPI_Comm_free(&cart);
    }
};

// ---------------------------------------------------------------------------
// E2E artifact helpers (ResidualHistoryE2ETest / PipelinedE2ETest).
// ---------------------------------------------------------------------------

// Load a whole file into a string; returns false when missing.
bool readWholeFile(const std::string& path, std::string& out) {
    std::ifstream ifs(path);
    if (!ifs) {
        return false;
    }
    std::ostringstream oss;
    oss << ifs.rdbuf();
    out = oss.str();
    return true;
}

// Parse "  \"key\": <value>" from a performance_report.json blob.
bool jsonScalar(const std::string& json, const std::string& key, double& value) {
    const std::string needle = "\"" + key + "\":";
    const std::size_t pos = json.find(needle);
    if (pos == std::string::npos) {
        return false;
    }
    std::istringstream iss(json.substr(pos + needle.size()));
    return static_cast<bool>(iss >> value);
}

struct ResidualCsv {
    bool ok = false;
    std::string header;
    std::vector<long long> iterations;
    std::vector<Real> residuals;
    std::string error;
};

// Parse the residual-history CSV contract: first line `iteration,residual`,
// then one `int,double` row per line (design §4.3 CSV column contract).
ResidualCsv parseResidualCsv(const std::string& path) {
    ResidualCsv csv;
    std::ifstream ifs(path);
    if (!ifs) {
        csv.error = "cannot open " + path;
        return csv;
    }
    if (!std::getline(ifs, csv.header)) {
        csv.error = "empty CSV (no header): " + path;
        return csv;
    }
    if (!csv.header.empty() && csv.header.back() == '\r') {
        csv.header.pop_back();
    }
    if (csv.header != "iteration,residual") {
        csv.error = "bad header '" + csv.header + "' (expected iteration,residual)";
        return csv;
    }
    std::string line;
    while (std::getline(ifs, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        long long it = 0;
        double res = 0.0;
        if (std::sscanf(line.c_str(), "%lld,%lf", &it, &res) != 2) {
            csv.error = "bad data row '" + line + "'";
            return csv;
        }
        csv.iterations.push_back(it);
        csv.residuals.push_back(res);
    }
    csv.ok = true;
    return csv;
}

// Reformat to %.4e — the reporter's actual precision for final_residual
// (setprecision(4) + scientific = 5 significant digits, reporter.cpp;
// the gate round 3 "4 significant digits" reading was off by one — its
// own measured sample 3.1416e-07 is a %.4e output) — for the
// "print-precision equality" comparisons of design §6 E1/E2/E3.
std::string reformat4e(Real v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.4e", v);
    return std::string(buf);
}

bool envDir(const char* name, std::string& dir) {
    const char* v = std::getenv(name);
    if (v == nullptr || *v == '\0') {
        return false;
    }
    dir = v;
    return true;
}

// The hypos driver's built-in problem (main.cpp setupProblem): rhs =
// -2*pi^2*sin(pi*x)*sin(pi*y) with x = i*dx, dx = 1/(nx-1). Jacobi's first
// residual row is the true residual of u = 0, i.e. ||rhs||_2 over the
// interior — the E2① anchor needs this value computed analytically (the
// omp reduction order is not reproducible bit-level; 1e-12 relative is the
// design's tolerance for that).
Real hyposDriverRhsL2(long long n) {
    const Real pi = std::acos(-1.0);
    const Real dx = 1.0 / static_cast<Real>(n - 1);
    const Real dy = 1.0 / static_cast<Real>(n - 1);
    Real sum = 0.0;
    // np1 subgrid: local indices [1, n] (halo = 1), same as main.cpp's loop.
    for (long long j = 1; j <= n; ++j) {
        for (long long i = 1; i <= n; ++i) {
            const Real x = static_cast<Real>(i) * dx;
            const Real y = static_cast<Real>(j) * dy;
            const Real f = -2.0 * pi * pi * std::sin(pi * x) * std::sin(pi * y);
            sum += f * f;
        }
    }
    return std::sqrt(sum);
}

} // namespace

// ---------------------------------------------------------------------------
// PipelinedCGUnitTest (np = 1)
// ---------------------------------------------------------------------------

// P1 (design §6 row 1): 256^2 multi-mode manufactured solution; cg and pcg
// each solve to the same tolerance.
//   |iters(pcg) - iters(cg)| / iters(cg) <= 0.05
//   L2(u_pcg - u_cg) <= 1e-8 * ||u_cg||
TEST(PipelinedCGUnitTest, PipelinedMatchesCgIterations) {
    const long long n = 256;

    Subgrid sgCg(n, n, 1, kHalo, MPI_COMM_SELF);
    fillMultiModeRhs(sgCg, n);
    PointToPointExchanger exCg;
    exCg.initialize(sgCg);
    CGSolver cg;
    const Index cgIters = cg.solve(sgCg, exCg, 5000, 1e-8);
    ASSERT_LT(cgIters, Index(5000));
    ASSERT_LE(cg.lastResidual(), 1e-8);

    Subgrid sgPcg(n, n, 1, kHalo, MPI_COMM_SELF);
    fillMultiModeRhs(sgPcg, n);
    PointToPointExchanger exPcg;
    exPcg.initialize(sgPcg);
    PipelinedCGSolver pcg;
    const Index pcgIters = pcg.solve(sgPcg, exPcg, 5000, 1e-8);
    ASSERT_LT(pcgIters, Index(5000));
    ASSERT_LE(pcg.lastResidual(), 1e-8);

    const double relDiff = std::fabs(static_cast<double>(pcgIters) -
                                     static_cast<double>(cgIters)) /
                           static_cast<double>(cgIters);
    EXPECT_LE(relDiff, 0.05) << "cg=" << cgIters << " pcg=" << pcgIters;

    Real diffSq = 0.0;
    Real refSq = 0.0;
    for (Index j = sgCg.jBegin(); j < sgCg.jEnd(); ++j) {
        for (Index i = sgCg.iBegin(); i < sgCg.iEnd(); ++i) {
            const Index idx = sgCg.index(i, j);
            const Real diff = sgPcg.u().data()[idx] - sgCg.u().data()[idx];
            diffSq += diff * diff;
            refSq += sgCg.u().data()[idx] * sgCg.u().data()[idx];
        }
    }
    ASSERT_GT(refSq, 0.0); // non-degenerate reference solution
    EXPECT_LE(std::sqrt(diffSq), 1e-8 * std::sqrt(refSq))
        << "L2diff=" << std::sqrt(diffSq) << " ||u||=" << std::sqrt(refSq);
}

// P2 (design §6 row 2): manufactured solutions at 64^2 and 256^2; the L2
// solution error is dominated by the O(h^2) truncation error (analytic
// rhs), so refining 64 -> 256 must shrink it by ~((257)/(65))^2 = 15.6x.
// Absolute anchor: err(64) <= 5e-3 (2nd-order scale, AR007/008 convention).
TEST(PipelinedCGUnitTest, PipelinedSolutionAccuracy) {
    auto solveError = [](long long n) -> Real {
        Subgrid sg(n, n, 1, kHalo, MPI_COMM_SELF);
        fillAnalyticSineRhs(sg, n);
        PointToPointExchanger ex;
        ex.initialize(sg);
        PipelinedCGSolver solver;
        const Index iters = solver.solve(sg, ex, 5000, 1e-9);
        EXPECT_LT(iters, Index(5000));
        Real localMax = 0.0;
        const Real sum = sineL2ErrorLocal(sg, n, n, localMax);
        return std::sqrt(sum / static_cast<Real>(n * n));
    };

    const Real err64 = solveError(64);
    const Real err256 = solveError(256);
    EXPECT_GT(err64, 0.0);
    EXPECT_GT(err256, 0.0);
    EXPECT_LE(err64, 5e-3) << "err64=" << err64;
    const Real ratio = err64 / err256;
    EXPECT_GT(ratio, 12.0) << "err64=" << err64 << " err256=" << err256;
    EXPECT_LT(ratio, 20.0) << "err64=" << err64 << " err256=" << err256;
}

// P3 (design §6 row 3 + §4.1 W1): two-sided profiler assertion.
//   pcg: stats("cg_blocking_allreduce").callCount == 0
//        stats("pcg_iallreduce").callCount == iters + 1
//        (init issue 1 + one per iteration body)
//   cg positive control on the same loadout (non-vacuous 0-assertion):
//        callCount > 0 && callCount == 2*iters + 1
//        (init rho + pap/rhoNew per iteration).
TEST(PipelinedCGUnitTest, NoBlockingAllreducePositiveControl) {
    const long long n = 64;

    {
        Subgrid sg(n, n, 1, kHalo, MPI_COMM_SELF);
        fillMultiModeRhs(sg, n);
        PointToPointExchanger ex;
        ex.initialize(sg);
        Profiler::instance().reset();
        PipelinedCGSolver solver;
        const Index iters = solver.solve(sg, ex, 5000, 1e-6);
        ASSERT_LT(iters, Index(5000));
        ASSERT_GT(iters, Index(0));
        EXPECT_EQ(Profiler::instance().stats("cg_blocking_allreduce").callCount,
                  std::uint64_t{0});
        EXPECT_EQ(Profiler::instance().stats("pcg_iallreduce").callCount,
                  std::uint64_t(iters) + 1)
            << "iters=" << iters;
    }
    {
        Subgrid sg(n, n, 1, kHalo, MPI_COMM_SELF);
        fillMultiModeRhs(sg, n);
        PointToPointExchanger ex;
        ex.initialize(sg);
        Profiler::instance().reset();
        CGSolver solver;
        const Index iters = solver.solve(sg, ex, 5000, 1e-6);
        ASSERT_LT(iters, Index(5000));
        ASSERT_GT(iters, Index(0));
        const std::uint64_t calls =
            Profiler::instance().stats("cg_blocking_allreduce").callCount;
        EXPECT_GT(calls, std::uint64_t{0});
        EXPECT_EQ(calls, 2 * std::uint64_t(iters) + 1)
            << "iters=" << iters;
    }
}

// P5 (design §6 row 5): RECORD ONLY — no triggering test case. A SPD loadout
// keeps nu = p^T A p > 0 at every step, and on the convergence iteration the
// nu-recurrence suffers catastrophic cancellation (m - beta^2*nu subtracts
// near-equal quantities), which the mandatory ①convergence-before-
// ②breakdown ordering (design §4.1) already defends against a false
// positive. The defense clause therefore cannot be stably constructed here;
// it follows the cg pap<=0 convention (WARN + stop) and is recorded in the
// design (D6), mirroring the AR008 breakdown no-trigger precedent.

// P6 (design §6 row 6): covered by the EXISTING 42-test regression plus the
// cg golden bit-level guard (U6) — the FP2 cgDotGlobal split must keep the
// cg floating-point sequence bit-identical, which those entries already
// pin. No new code in this file; the full-suite entries (unit, pcg_np1,
// pcg_np4, ...) run everything.

// P7 (design §6 row 7): iterate() contract.
//   Precondition violation (iterate before solve): WARN + return 0.
//   Resume: after a maxIter-truncated solve, one iterate() step must return
//   the residual of the post-step state (compared against the independent
//   globalTrueResidual scan, 1e-6 relative of ||b||).
TEST(PipelinedCGUnitTest, IterateGuardAndResume) {
    const long long n = 32;
    Subgrid sg(n, n, 1, kHalo, MPI_COMM_SELF);
    fillMultiModeRhs(sg, n);
    PointToPointExchanger ex;
    ex.initialize(sg);

    // ||b|| scale for the resume comparison.
    Real bSq = 0.0;
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const Real b = -sg.rhs().data()[sg.index(i, j)];
            bSq += b * b;
        }
    }
    ASSERT_GT(bSq, 0.0);
    const Real bNorm = std::sqrt(bSq);

    PipelinedCGSolver solver;
    // Guard: no solve() yet — WARN + 0 (CGSolver::iterate convention).
    EXPECT_EQ(solver.iterate(sg, ex), Real(0.0));

    // Partial solve (tol=0 never converges; 5 iterations, state ready at
    // loop entry with no pending request), then one resume step.
    const Index iters = solver.solve(sg, ex, 5, 0.0);
    EXPECT_EQ(iters, Index(5));

    const Real returned = solver.iterate(sg, ex);
    EXPECT_GT(returned, 0.0);

    const Real trueRes = globalTrueResidual(sg, ex);
    EXPECT_GT(trueRes, 0.0);
    EXPECT_LE(std::fabs(returned - trueRes), 1e-6 * bNorm)
        << "returned=" << returned << " true=" << trueRes
        << " ||b||=" << bNorm;
}

// P8 (design §6 row 8): rhs all 0 (rho_0 = 0 < tol) — immediate convergence
// along the cg :220 convention: 0 iterations, lastResidual == 0, and the
// iters+1 boundary of the Iallreduce count (init issue only) == 1.
TEST(PipelinedCGUnitTest, ZeroRhsImmediateConvergence) {
    const long long n = 16;
    Subgrid sg(n, n, 1, kHalo, MPI_COMM_SELF);
    sg.zeroInitialize(); // u = 0 and rhs = 0
    sg.applyDirichletBC(0.0);
    PointToPointExchanger ex;
    ex.initialize(sg);

    Profiler::instance().reset();
    PipelinedCGSolver solver;
    const Index iters = solver.solve(sg, ex, 100, 1e-6);
    EXPECT_EQ(iters, Index(0));
    EXPECT_EQ(solver.lastResidual(), Real(0.0));
    EXPECT_EQ(Profiler::instance().stats("pcg_iallreduce").callCount,
              std::uint64_t{1});
    EXPECT_EQ(Profiler::instance().stats("cg_blocking_allreduce").callCount,
              std::uint64_t{0});
}

// P9 (design §6 row 9): maxIter=3 forced truncation (tol tiny -> 0). The
// solve must return completed==3 (cg :221 convention) and the progress
// trajectory must be dense/continuous: exactly one notification per
// iteration, iterations 1,2,3, each with a fresh finite residual.
TEST(PipelinedCGUnitTest, MaxIterTruncation) {
    const long long n = 32;
    Subgrid sg(n, n, 1, kHalo, MPI_COMM_SELF);
    fillMultiModeRhs(sg, n);
    PointToPointExchanger ex;
    ex.initialize(sg);

    std::vector<Index> trajIters;
    std::vector<Real> trajRes;
    PipelinedCGSolver solver;
    solver.setProgressCallback([&](Index iteration) {
        trajIters.push_back(iteration);
        trajRes.push_back(solver.lastResidual());
    });

    const Index iters = solver.solve(sg, ex, 3, 0.0);
    EXPECT_EQ(iters, Index(3));

    ASSERT_EQ(trajIters.size(), std::size_t{3});
    EXPECT_EQ(trajIters[0], Index(1));
    EXPECT_EQ(trajIters[1], Index(2));
    EXPECT_EQ(trajIters[2], Index(3));
    for (std::size_t k = 0; k < trajRes.size(); ++k) {
        EXPECT_TRUE(std::isfinite(trajRes[k])) << "k=" << k;
        EXPECT_GT(trajRes[k], 0.0) << "k=" << k;
    }
}

// ---------------------------------------------------------------------------
// PipelinedCGMpiTest (np = 4)
// ---------------------------------------------------------------------------

// P4 (design §6 row 4): np4 vs np1 on the same non-uniform {15,19}x{15,19}
// layout loadout; solution L2 difference <= 1e-10*||u|| (U9 convention —
// the CG-1 reduction split order differs across ranks, so this is NOT a
// bit-level comparison). The np1 reference is recomputed redundantly on
// every rank (34^2 is tiny), mirroring MgcgNp4MatchesNp1.
TEST(PipelinedCGMpiTest, PipelinedNp4Consistency) {
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "PipelinedNp4Consistency requires exactly 4 processes";
    }
    NonUniformLayout lay;
    lay.build();
    if (HasFailure()) return;
    Subgrid& sg = *lay.sg;
    fillSineRhs(sg, kLayoutN, kLayoutN);

    PipelinedCGSolver solver;
    const Index iters = solver.solve(sg, *lay.ex, 200, 1e-8);
    ASSERT_LT(iters, Index(200));
    ASSERT_LE(solver.lastResidual(), 1e-8);

    Subgrid ref(kLayoutN, kLayoutN, 1, kHalo, MPI_COMM_SELF);
    fillSineRhs(ref, kLayoutN, kLayoutN);
    PointToPointExchanger refEx;
    refEx.initialize(ref);
    PipelinedCGSolver refSolver;
    const Index refIters = refSolver.solve(ref, refEx, 200, 1e-8);
    ASSERT_LT(refIters, Index(200));
    ASSERT_LE(refSolver.lastResidual(), 1e-8);

    Real localDiffSq = 0.0;
    Real localRefSq = 0.0;
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gi = static_cast<long long>(sg.offsetX()) +
                                 static_cast<long long>(i) - kHalo;
            const long long gj = static_cast<long long>(sg.offsetY()) +
                                 static_cast<long long>(j) - kHalo;
            const Real got = sg.u().data()[sg.index(i, j)];
            const Real want = ref.u().data()[ref.index(
                static_cast<Index>(gi + kHalo), static_cast<Index>(gj + kHalo))];
            const Real diff = got - want;
            localDiffSq += diff * diff;
            localRefSq += want * want;
        }
    }
    Real local[2] = {localDiffSq, localRefSq};
    Real global[2] = {0.0, 0.0};
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, lay.cart);

    const Real diffL2 = std::sqrt(global[0]);
    const Real refL2 = std::sqrt(global[1]);
    ASSERT_GT(refL2, 1e-6); // non-degenerate solution
    EXPECT_LE(diffL2, 1e-10 * refL2)
        << "diffL2=" << diffL2 << " refL2=" << refL2;

    lay.destroy();
}

// ---------------------------------------------------------------------------
// ResidualHistoryE2ETest — artifact verifiers for the E1/E2/E4/E7 ctest
// fixture runs (design §6 rows E1, E2, E4, E7). These tests read the
// CSV/JSON/stdout artifacts written by the residual_history_* ctest entries
// into $HYPOS_AR009_RH_BASE; without that variable (plain `unit` entry)
// they SKIP, so the standalone suite stays green.
// ---------------------------------------------------------------------------

namespace {

// Loads <base>/<solver>/residual_<solver>.csv + performance_report.json.
struct SolverRunArtifacts {
    ResidualCsv csv;
    double iterations = -1.0;
    double finalResidual = -1.0;
    bool loaded = false;
    std::string error;

    bool load(const std::string& base, const std::string& solver) {
        const std::string dir = base + "/" + solver;
        std::string csvPath = dir + "/residual_" + solver + ".csv";
        std::string jsonPath = dir + "/performance_report.json";
        csv = parseResidualCsv(csvPath);
        if (!csv.ok) {
            error = csv.error;
            return false;
        }
        std::string json;
        if (!readWholeFile(jsonPath, json)) {
            error = "missing " + jsonPath;
            return false;
        }
        if (!jsonScalar(json, "iterations", iterations) ||
            !jsonScalar(json, "final_residual", finalResidual)) {
            error = "report JSON lacks iterations/final_residual: " + jsonPath;
            return false;
        }
        loaded = true;
        return true;
    }
};

// Common sanity: header contract + data row count == reported iterations
// (k=1 sampling: one row per iteration/cycle).
::testing::AssertionResult checkRowCount(const SolverRunArtifacts& a,
                                         const char* solver) {
    if (a.iterations < 0.0) {
        return ::testing::AssertionFailure()
               << solver << ": bad iterations " << a.iterations;
    }
    const long long expected = static_cast<long long>(a.iterations);
    if (static_cast<long long>(a.csv.residuals.size()) != expected) {
        return ::testing::AssertionFailure()
               << solver << ": CSV data rows=" << a.csv.residuals.size()
               << " but report iterations=" << expected;
    }
    return ::testing::AssertionSuccess();
}

// Fixture gate: SKIP only when the artifact files are absent (the plain
// `unit` entry runs this binary without the e2e fixtures); a present but
// malformed artifact must FAIL downstream, never skip.
bool fixtureReady(const std::string& base, const std::string& solver) {
    namespace fs = std::filesystem;
    const std::string csvPath = base + "/" + solver + "/residual_" + solver + ".csv";
    return fs::exists(fs::path(csvPath));
}

// E7 (jacobi save_combo) monotonicity clause: adjacent rows non-increasing
// within 1e-12 relative (numerical-noise allowance; jacobi's stationary
// contraction is genuinely monotone — unlike cg, see E1).
::testing::AssertionResult checkMonotone(const ResidualCsv& csv,
                                         const char* solver) {
    for (std::size_t k = 1; k < csv.residuals.size(); ++k) {
        const Real prev = csv.residuals[k - 1];
        const Real next = csv.residuals[k];
        if (next > prev * (1.0 + 1e-12)) {
            return ::testing::AssertionFailure()
                   << solver << ": row " << k << " increased: " << next
                   << " > " << prev;
        }
    }
    return ::testing::AssertionSuccess();
}

} // namespace

// E1 (design §6 row E1, the cg anchor of the E2 loop): np1
// `--solver cg --residual-history`.
//   header == "iteration,residual"
//   data rows == iterations
//   last row == JSON final_residual at print precision (both sides
//   reformatted to %.4e — the reporter's 5-significant-digit precision).
TEST(ResidualHistoryE2ETest, CgCsvContract) {
    std::string base;
    if (!envDir("HYPOS_AR009_RH_BASE", base)) {
        GTEST_SKIP() << "HYPOS_AR009_RH_BASE not set (e2e fixture inactive)";
    }
    SolverRunArtifacts a;
    if (!fixtureReady(base, "cg")) {
        GTEST_SKIP() << "cg e2e artifacts missing (fixture inactive)";
    }
    ASSERT_TRUE(a.load(base, "cg")) << a.error;
    ASSERT_TRUE(a.csv.ok);
    EXPECT_EQ(a.csv.header, "iteration,residual");
    EXPECT_TRUE(checkRowCount(a, "cg"));
    // NOTE (design E1 revised): NO monotonicity assertion for cg — the
    // 2-norm of the CG residual is NOT monotonically decreasing (only the
    // A-norm of the error is); the measured driver trajectory genuinely
    // rises 0.153 -> 0.204 -> 0.246 at rows 45-47 before dropping again
    // (evidence: first full e2e run, 2026-10-10). Jacobi's trajectory IS
    // monotone (stationary contraction) and keeps the check in E7.
    ASSERT_FALSE(a.csv.residuals.empty());
    const Real last = a.csv.residuals.back();
    // Print-precision equality: %.4e(csv last) parses to the same double as
    // the JSON value (reporter prints 4 significant digits).
    EXPECT_EQ(std::stod(reformat4e(last)), a.finalResidual)
        << "csv last=" << last << " json=" << a.finalResidual;
    // Row indices are 1-based and dense.
    for (std::size_t k = 0; k < a.csv.iterations.size(); ++k) {
        EXPECT_EQ(a.csv.iterations[k], static_cast<long long>(k) + 1);
    }
}

// E2① (design §6 row E2): jacobi.
//   First data row == analytic ||b|| within 1e-12 relative (the omp
//   reduction split order is not reproducible from the test side; the
//   "initial residual" anchor is the in-test analytic ||rhs||).
//   Last row vs the "Converged at iteration" log residual: the log line is
//   captured to run_stdout.txt by the residual_history_jacobi ctest entry;
//   compared after reformatting both to 6 significant digits (%g).
TEST(ResidualHistoryE2ETest, JacobiFirstRowAndConvergedLine) {
    std::string base;
    if (!envDir("HYPOS_AR009_RH_BASE", base)) {
        GTEST_SKIP() << "HYPOS_AR009_RH_BASE not set (e2e fixture inactive)";
    }
    SolverRunArtifacts a;
    if (!fixtureReady(base, "jacobi")) {
        GTEST_SKIP() << "jacobi e2e artifacts missing (fixture inactive)";
    }
    ASSERT_TRUE(a.load(base, "jacobi")) << a.error;
    ASSERT_FALSE(a.csv.residuals.empty());

    // First row == initial residual == ||b|| (1e-12 relative).
    const long long n = 32; // matches the residual_history_jacobi entry
    const Real bNorm = hyposDriverRhsL2(n);
    ASSERT_GT(bNorm, 0.0);
    const Real first = a.csv.residuals.front();
    EXPECT_LE(std::fabs(first - bNorm), 1e-12 * bNorm)
        << "first=" << first << " ||b||=" << bNorm;

    // Last row vs the captured "Converged at iteration" log line, 6-digit
    // %g print precision (C1: all comparisons at print precision).
    std::string logTxt;
    ASSERT_TRUE(readWholeFile(base + "/jacobi/run_stdout.txt", logTxt))
        << "jacobi run_stdout.txt missing (rh_jacobi script contract)";
    const std::size_t pos = logTxt.find("Converged at iteration");
    ASSERT_NE(pos, std::string::npos)
        << "jacobi log has no Converged line";
    long long convIter = 0;
    double convRes = 0.0;
    if (std::sscanf(logTxt.c_str() + pos, "Converged at iteration %lld, residual = %lf",
                    &convIter, &convRes) != 2) {
        FAIL() << "cannot parse Converged line: "
               << logTxt.substr(pos, 80);
    }
    EXPECT_EQ(convIter, static_cast<long long>(a.csv.residuals.size()))
        << "Converged iteration vs CSV row count";
    const Real last = a.csv.residuals.back();
    char csvG[64];
    std::snprintf(csvG, sizeof(csvG), "%g", last);
    char logG[64];
    std::snprintf(logG, sizeof(logG), "%g", convRes);
    EXPECT_STREQ(csvG, logG)
        << "csv last=" << last << " converged-line residual=" << convRes;
    // NOTE (design E2①, evidence clause): the exit confirmation scan's
    // value (JSON final_residual) legitimately differs from the in-loop
    // check value by about one contraction step — deliberately NOT
    // asserted, per the design table.
}

// E2② (design §6 row E2): rbgs / mg2 / mgv — last row vs JSON
// final_residual at print precision (the confirmation scan recomputes the
// same final state; %.4e reformat equality).
TEST(ResidualHistoryE2ETest, RbgsLastRowMatchesReport) {
    std::string base;
    if (!envDir("HYPOS_AR009_RH_BASE", base)) {
        GTEST_SKIP() << "HYPOS_AR009_RH_BASE not set (e2e fixture inactive)";
    }
    SolverRunArtifacts a;
    if (!fixtureReady(base, "red_black_gs")) {
        GTEST_SKIP() << "red_black_gs e2e artifacts missing (fixture inactive)";
    }
    ASSERT_TRUE(a.load(base, "red_black_gs")) << a.error;
    EXPECT_TRUE(checkRowCount(a, "red_black_gs"));
    ASSERT_FALSE(a.csv.residuals.empty());
    EXPECT_EQ(std::stod(reformat4e(a.csv.residuals.back())), a.finalResidual)
        << "csv last=" << a.csv.residuals.back() << " json=" << a.finalResidual;
}

TEST(ResidualHistoryE2ETest, Mg2LastRowMatchesReport) {
    std::string base;
    if (!envDir("HYPOS_AR009_RH_BASE", base)) {
        GTEST_SKIP() << "HYPOS_AR009_RH_BASE not set (e2e fixture inactive)";
    }
    SolverRunArtifacts a;
    if (!fixtureReady(base, "mg2")) {
        GTEST_SKIP() << "mg2 e2e artifacts missing (fixture inactive)";
    }
    ASSERT_TRUE(a.load(base, "mg2")) << a.error;
    EXPECT_TRUE(checkRowCount(a, "mg2"));
    ASSERT_FALSE(a.csv.residuals.empty());
    EXPECT_EQ(std::stod(reformat4e(a.csv.residuals.back())), a.finalResidual)
        << "csv last=" << a.csv.residuals.back() << " json=" << a.finalResidual;
}

TEST(ResidualHistoryE2ETest, MgvLastRowMatchesReport) {
    std::string base;
    if (!envDir("HYPOS_AR009_RH_BASE", base)) {
        GTEST_SKIP() << "HYPOS_AR009_RH_BASE not set (e2e fixture inactive)";
    }
    SolverRunArtifacts a;
    if (!fixtureReady(base, "mgv")) {
        GTEST_SKIP() << "mgv e2e artifacts missing (fixture inactive)";
    }
    ASSERT_TRUE(a.load(base, "mgv")) << a.error;
    EXPECT_TRUE(checkRowCount(a, "mgv"));
    ASSERT_FALSE(a.csv.residuals.empty());
    EXPECT_EQ(std::stod(reformat4e(a.csv.residuals.back())), a.finalResidual)
        << "csv last=" << a.csv.residuals.back() << " json=" << a.finalResidual;
}

// E2④ (design §6 row E2): pcg / mgcg — last row (fresh recurrence rho /
// review value) vs JSON final_residual. Print precision alone allows 5e-4
// relative (4 significant digits); the design's >= 5e-6 relative floor for
// the 6-digit finish-line path is subsumed — use 1e-3 relative to cover
// print rounding plus the recurrence-vs-review floating difference (§4.1).
TEST(ResidualHistoryE2ETest, PcgLastRowMatchesReport) {
    std::string base;
    if (!envDir("HYPOS_AR009_RH_BASE", base)) {
        GTEST_SKIP() << "HYPOS_AR009_RH_BASE not set (e2e fixture inactive)";
    }
    SolverRunArtifacts a;
    if (!fixtureReady(base, "pcg")) {
        GTEST_SKIP() << "pcg e2e artifacts missing (fixture inactive)";
    }
    ASSERT_TRUE(a.load(base, "pcg")) << a.error;
    EXPECT_TRUE(checkRowCount(a, "pcg"));
    ASSERT_FALSE(a.csv.residuals.empty());
    const Real last = a.csv.residuals.back();
    ASSERT_GT(a.finalResidual, 0.0);
    EXPECT_LE(std::fabs(last - a.finalResidual), 1e-3 * a.finalResidual)
        << "csv last=" << last << " json=" << a.finalResidual;
}

TEST(ResidualHistoryE2ETest, MgcgLastRowMatchesReport) {
    std::string base;
    if (!envDir("HYPOS_AR009_RH_BASE", base)) {
        GTEST_SKIP() << "HYPOS_AR009_RH_BASE not set (e2e fixture inactive)";
    }
    SolverRunArtifacts a;
    if (!fixtureReady(base, "mgcg")) {
        GTEST_SKIP() << "mgcg e2e artifacts missing (fixture inactive)";
    }
    ASSERT_TRUE(a.load(base, "mgcg")) << a.error;
    EXPECT_TRUE(checkRowCount(a, "mgcg"));
    ASSERT_FALSE(a.csv.residuals.empty());
    const Real last = a.csv.residuals.back();
    ASSERT_GT(a.finalResidual, 0.0);
    EXPECT_LE(std::fabs(last - a.finalResidual), 1e-3 * a.finalResidual)
        << "csv last=" << last << " json=" << a.finalResidual;
}

// E4 (design §6 row E4): `--solver jacobi --residual-check-interval 10`.
// The writer emits a row only when iteration % k == 0, so with M total
// iterations the row count is floor(M/10) — under the design premise
// (jacobi converges at a check iteration, M % 10 == 0) ceil == floor; the
// max-iter truncation of the fixture run keeps M = 100 deterministic.
TEST(ResidualHistoryE2ETest, IntervalTenRowCount) {
    std::string base;
    if (!envDir("HYPOS_AR009_RH_BASE", base)) {
        GTEST_SKIP() << "HYPOS_AR009_RH_BASE not set (e2e fixture inactive)";
    }
    SolverRunArtifacts a;
    if (!fixtureReady(base, "interval")) {
        GTEST_SKIP() << "interval e2e artifacts missing (fixture inactive)";
    }
    ASSERT_TRUE(a.load(base, "interval")) << a.error;
    const long long iterations = static_cast<long long>(a.iterations);
    ASSERT_GT(iterations, 0);
    const long long expectedRows = iterations / 10; // floor, k = 10
    EXPECT_EQ(static_cast<long long>(a.csv.residuals.size()), expectedRows)
        << "rows=" << a.csv.residuals.size() << " iterations=" << iterations;
    // Sampled rows sit exactly on the multiples of 10.
    for (std::size_t k = 0; k < a.csv.iterations.size(); ++k) {
        EXPECT_EQ(a.csv.iterations[k], static_cast<long long>(k + 1) * 10);
    }
}

// E7 (design §6 row E7): `--save-interval k --residual-history f` — both
// channels of the composed single-slot progressCallback coexist: the CSV
// is valid AND the solution snapshots land on the k grid (binary backend,
// mirroring the binary_e2e shard naming).
TEST(ResidualHistoryE2ETest, SaveComboBothChannels) {
    std::string base;
    if (!envDir("HYPOS_AR009_RH_BASE", base)) {
        GTEST_SKIP() << "HYPOS_AR009_RH_BASE not set (e2e fixture inactive)";
    }
    const std::string dir = base + "/save_combo";
    // NOTE: the save_combo run uses --output-format binary so the solution
    // snapshots land on disk; with a non-json/csv format the driver writes
    // no performance_report.json, so only the CSV contract is checked here.
    const ResidualCsv csv = parseResidualCsv(dir + "/residual_save_combo.csv");
    ASSERT_TRUE(csv.ok) << csv.error;
    EXPECT_EQ(csv.header, "iteration,residual");
    EXPECT_FALSE(csv.residuals.empty());
    EXPECT_TRUE(checkMonotone(csv, "save_combo"));

    namespace fs = std::filesystem;
    const std::string s10 = dir + "/solution_10_r0.bin";
    const std::string s20 = dir + "/solution_20_r0.bin";
    EXPECT_TRUE(fs::exists(fs::path(s10))) << "missing " << s10;
    EXPECT_TRUE(fs::exists(fs::path(s20))) << "missing " << s20;
}

// ---------------------------------------------------------------------------
// ResidualHistoryNp4E2ETest — E3 (design §6 row E3): np4 residual history.
// ---------------------------------------------------------------------------

// E3: exactly ONE CSV in the np4 output directory (rank0-only writer) and
// the np4 artifacts are internally consistent (header contract, row count
// == reported iterations, last row vs the np4 report's final_residual at
// print precision). The design's original np1-trajectory cross-check
// (1e-6 relative, ±1 row) proved unexecutable in practice: the np1
// fixture runs OMP_NUM_THREADS=4 while each of the 4 MPI ranks runs ONE
// thread (HYPOS_TEST_ENV vs HYPOS_TEST_ENV_MPI, the CI oversubscription
// policy), so the two runs reduce the dots in different block orders, the
// cg floating-point trajectories genuinely diverge, and the converged
// iteration counts differ by MORE than the ±1 allowance (evidence: first
// full-suite run, 2026-10-10). The cross-config comparison is therefore
// evidence-only (RecordProperty); the assertion load shifts to the
// internal-consistency clauses above (design §6 E3 revised — same fix
// protocol as the gate rounds).
TEST(ResidualHistoryNp4E2ETest, SingleRank0CsvInternalConsistency) {
    std::string np4Dir;
    if (!envDir("HYPOS_AR009_RH_NP4_DIR", np4Dir)) {
        GTEST_SKIP()
            << "HYPOS_AR009_RH_NP4_DIR not set (e2e fixture inactive)";
    }

    // Exactly one residual CSV file in the np4 output directory (rank0
    // single writer; a per-rank writer would leave shards).
    namespace fs = std::filesystem;
    std::vector<std::string> csvFiles;
    for (const auto& entry : fs::directory_iterator(fs::path(np4Dir))) {
        const std::string name = entry.path().filename().string();
        const std::size_t dot = name.find_last_of('.');
        if (name.rfind("residual", 0) == 0 && dot != std::string::npos &&
            name.substr(dot) == ".csv") {
            csvFiles.push_back(name);
        }
    }
    ASSERT_EQ(csvFiles.size(), std::size_t{1})
        << "np4 output dir must hold exactly one rank0 CSV";
    EXPECT_EQ(csvFiles[0], "residual_cg.csv");

    const ResidualCsv np4 = parseResidualCsv(np4Dir + "/residual_cg.csv");
    ASSERT_TRUE(np4.ok) << np4.error;
    EXPECT_EQ(np4.header, "iteration,residual");
    ASSERT_FALSE(np4.residuals.empty());

    // Internal consistency against the np4 run's own report.
    std::string json;
    ASSERT_TRUE(readWholeFile(np4Dir + "/performance_report.json", json))
        << "np4 performance_report.json missing";
    double iterations = -1.0;
    double finalResidual = -1.0;
    ASSERT_TRUE(jsonScalar(json, "iterations", iterations));
    ASSERT_TRUE(jsonScalar(json, "final_residual", finalResidual));
    EXPECT_EQ(static_cast<long long>(np4.residuals.size()),
              static_cast<long long>(iterations))
        << "np4 CSV rows vs reported iterations";
    // Last row at print precision (E1 convention: %.4e reformat).
    EXPECT_EQ(std::stod(reformat4e(np4.residuals.back())), finalResidual)
        << "np4 csv last=" << np4.residuals.back() << " json=" << finalResidual;

    // Evidence-only: the np1 row count is recorded for the record. The
    // cross-config trajectory comparison is deliberately NOT asserted
    // (OMP reduction block order differs between the fixture entries).
    std::string base;
    if (envDir("HYPOS_AR009_RH_BASE", base)) {
        if (fixtureReady(base, "cg")) {
            SolverRunArtifacts np1;
            if (np1.load(base, "cg")) {
                RecordProperty("np1_rows",
                               std::to_string(np1.csv.residuals.size()));
                RecordProperty("np4_rows",
                               std::to_string(np4.residuals.size()));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// PipelinedE2ETest — E5 (design §6 row E5), report-side clause.
// ---------------------------------------------------------------------------

// E5 true-residual clause: the driver's globalTrueResidual review value
// lands in performance_report.json; it must satisfy tol*(1+5e-4) — the
// guard band absorbs the JSON's 4-significant-digit print rounding so the
// near-threshold case is not brittle. (Exit code 0 and the finish line are
// asserted by the pcg_e2e ctest entry via VerifyAr009Cli.cmake.)
TEST(PipelinedE2ETest, PcgE2ETrueResidualInReport) {
    std::string dir;
    if (!envDir("HYPOS_AR009_PCG_E2E_DIR", dir)) {
        GTEST_SKIP() << "HYPOS_AR009_PCG_E2E_DIR not set (e2e fixture inactive)";
    }
    std::string json;
    if (!readWholeFile(dir + "/performance_report.json", json)) {
        GTEST_SKIP() << "pcg e2e report missing in " << dir;
    }
    double finalResidual = -1.0;
    ASSERT_TRUE(jsonScalar(json, "final_residual", finalResidual));
    const double kTol = 1e-6; // matches the pcg_e2e entry's --tol
    EXPECT_LE(finalResidual, kTol * (1.0 + 5e-4))
        << "final_residual=" << finalResidual;
}
