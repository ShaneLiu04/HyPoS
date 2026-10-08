#include <gtest/gtest.h>
#include "solver/solver.hpp"
// AR004 T001: real-residual kernel; the header is created by the implementer.
// Until then this include fails compilation, which is the intended Red state.
#include "solver/residual.hpp"
#include "grid/subgrid.hpp"
#include "comm/halo_exchanger.hpp"
#include "perf/profiler.hpp"
#include <mpi.h>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

using namespace hypo;

TEST(SolverTest, JacobiConvergence) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    (void)rank;

    // Create a small 2D subgrid (no neighbors, so use MPI_COMM_SELF)
    Subgrid subgrid(32, 32, 1, 1, MPI_COMM_SELF);
    subgrid.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);

    // Set up a simple problem: u = 0 on boundary, lap(u) = -1 everywhere.
    // The solver convention is lap(u) = rhs, so rhs = -1 yields a positive
    // interior solution (maximum principle).
    subgrid.zeroInitialize();
    subgrid.applyDirichletBC(0.0);

    Real* rhs = subgrid.rhs().data();
    for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
        for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
            rhs[subgrid.index(i, j)] = -1.0;
        }
    }

    PointToPointExchanger exchanger;
    exchanger.initialize(subgrid);

    JacobiSolver solver;
    Index iters = solver.solve(subgrid, exchanger, 5000, 1e-5);

    EXPECT_GT(iters, 0);
    EXPECT_LT(iters, 5000); // Should converge before max iterations

    // Check that interior values are positive (since lap(u) = -1 and boundary = 0)
    const Real* u = subgrid.u().data();
    for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
        for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
            EXPECT_GT(u[subgrid.index(i, j)], 0.0);
        }
    }
}

TEST(SolverApiTest, OverlapFlagDefaultsAndGetter) {
    JacobiSolver defaultSolver;
    EXPECT_FALSE(defaultSolver.overlapEnabled());

    JacobiSolver overlapSolver(true);
    EXPECT_TRUE(overlapSolver.overlapEnabled());

    EXPECT_EQ(overlapSolver.name(), "jacobi");
}

TEST(SolverApiTest, LastResidualZeroBeforeSolve) {
    JacobiSolver solver;
    EXPECT_DOUBLE_EQ(solver.lastResidual(), 0.0);
}

TEST(SolverApiTest, ZeroMaxIterations) {
    Subgrid sg(16, 16, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);

    PointToPointExchanger ex;
    ex.initialize(sg);

    JacobiSolver solver;
    Index iters = solver.solve(sg, ex, 0, 1e-6);

    EXPECT_EQ(iters, Index(0));
    EXPECT_DOUBLE_EQ(solver.lastResidual(), 0.0);
}

TEST(SolverApiTest, FirstIterationConvergenceSemantics) {
    Subgrid sg(8, 8, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    sg.zeroInitialize(); // u = 0 and rhs = 0, so the first update is already converged

    PointToPointExchanger ex;
    ex.initialize(sg);

    JacobiSolver solver;
    // solve() returns the number of executed iterations: convergence detected
    // after the first update must report 1, not 0.
    Index iters = solver.solve(sg, ex, 100, std::numeric_limits<Real>::infinity());

    EXPECT_EQ(iters, Index(1));
    EXPECT_TRUE(std::isfinite(solver.lastResidual()));
}

TEST(SolverApiTest, LastResidualAfterConvergence) {
    Subgrid sg(32, 32, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);

    Real* rhs = sg.rhs().data();
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            rhs[sg.index(i, j)] = -1.0;
        }
    }

    PointToPointExchanger ex;
    ex.initialize(sg);

    JacobiSolver solver;
    Index iters = solver.solve(sg, ex, 20000, 1e-5);

    EXPECT_GT(iters, Index(0));
    EXPECT_LT(iters, Index(20000));
    EXPECT_GT(solver.lastResidual(), 0.0);
    EXPECT_LE(solver.lastResidual(), 1e-5);
}

namespace {

class MinimalSolver : public PoissonSolver {
public:
    Index solve(Subgrid&, HaloExchanger&, Index, Real) override { return 0; }
    Real iterate(Subgrid&, HaloExchanger&) override { return 0.0; }
    std::string name() const override { return "minimal"; }
};

} // namespace

TEST(SolverApiTest, BaseClassLastResidualDefault) {
    // MinimalSolver does not override lastResidual(); the base-class default
    // must return 0 without crashing.
    MinimalSolver s;
    EXPECT_DOUBLE_EQ(s.lastResidual(), 0.0);
}

TEST(SolverApiTest, NonConvergenceAfterMaxIterations) {
    Subgrid sg(16, 16, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);

    Real* rhs = sg.rhs().data();
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            rhs[sg.index(i, j)] = -1.0;
        }
    }

    PointToPointExchanger ex;
    ex.initialize(sg);
    JacobiSolver solver;

    // Runs out of iteration budget before reaching the tight tolerance.
    Index iters = solver.solve(sg, ex, 1, 1e-12);
    EXPECT_EQ(iters, Index(1));
    EXPECT_GT(solver.lastResidual(), 1e-12);
}

// ============================================================================
// AR004 T001: A2-Jacobi real residual + residual check interval.
// These tests reference the new interface (trueResidualSquaredLocal,
// setResidualCheckInterval, converted iterate() return, residual_confirm
// region) before it exists: compilation failure here is the legal Red state.
// ============================================================================

namespace {

// Relative-error tolerance for residual comparisons: the OpenMP reduction in
// trueResidualSquaredLocal sums in a different order than the serial
// recomputation below, so results agree to ~1 ulp but not bit-exactly.
// The AR004 design fixes 1e-12 relative; EXPECT_DOUBLE_EQ is therefore wrong.
constexpr Real kAr004RelTol = 1e-12;

void expectCloseRelative(Real actual, Real expected, const std::string& context) {
    const Real tol = kAr004RelTol * std::fabs(expected);
    EXPECT_NEAR(actual, expected, tol) << context;
}

// Exact manufactured solution on an n x n global grid: u = sin(pi x) sin(pi y)
// with x = (gI+1)/(n+1). Ghost points gI = -1 and gI = n map to x = 0 / 1,
// where the sine vanishes, so a Dirichlet-0 halo is exact for single rank.
Real sineExact2D(long long gI, long long gJ, long long n) {
    const double x = static_cast<double>(gI + 1) / static_cast<double>(n + 1);
    const double y = static_cast<double>(gJ + 1) / static_cast<double>(n + 1);
    return std::sin(M_PI * x) * std::sin(M_PI * y);
}

// Single-rank manufactured sine setup: rhs = discrete Laplacian of the exact
// solution, following the setupManufacturedSine pattern in test_alt_solvers.
void setupManufacturedSine(Subgrid& sg, Index globalN) {
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);
    Real* rhs = sg.rhs().data();
    const long long hw = static_cast<long long>(sg.haloWidth());
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gI = static_cast<long long>(i) - hw;
            const long long gJ = static_cast<long long>(j) - hw;
            rhs[sg.index(i, j)] =
                sineExact2D(gI - 1, gJ, globalN) + sineExact2D(gI + 1, gJ, globalN) +
                sineExact2D(gI, gJ - 1, globalN) + sineExact2D(gI, gJ + 1, globalN) -
                4.0 * sineExact2D(gI, gJ, globalN);
        }
    }
}

// Same uniform rhs = -1 setup used by the existing solver tests.
void setupUniformNegativeRhs(Subgrid& sg) {
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);
    Real* rhs = sg.rhs().data();
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            rhs[sg.index(i, j)] = -1.0;
        }
    }
}

// Independent serial recomputation of the local true residual sum of squares
// for the Au = -rhs convention (r = D*u - sum(neighbors) + rhs), D = 4 in 2D,
// D = 6 in 3D. Takes an explicit u buffer so callers can evaluate a snapshot
// from before an iteration. Halo values are read as-is: the caller is
// responsible for having refreshed them (exchange / applyPhysicalBoundary).
Real trueResidualSumSqSerial(const Subgrid& sg, const Real* u) {
    const Real* rhs = sg.rhs().data();
    const Index nxT = sg.nxTotal();
    const Index nyT = sg.nyTotal();
    const bool is2D = (sg.nzLocal() == 1);
    const Real denom = is2D ? 4.0 : 6.0;
    // 2D lives in the k = 0 plane of the padded array.
    const Index k0 = is2D ? 0 : sg.kBegin();
    const Index k1 = is2D ? 1 : sg.kEnd();

    Real sum = 0.0;
    for (Index k = k0; k < k1; ++k) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                const Index idx = (k * nyT + j) * nxT + i;
                const Real neighborSum = u[idx - 1] + u[idx + 1] +
                                         u[idx - nxT] + u[idx + nxT] +
                                         (is2D ? 0.0
                                               : u[idx - nxT * nyT] +
                                                 u[idx + nxT * nyT]);
                const Real r = denom * u[idx] - neighborSum + rhs[idx];
                sum += r * r;
            }
        }
    }
    return sum;
}

} // namespace

// U1 (2D): trueResidualSquaredLocal must reproduce an independent serial
// recomputation on a small hand-filled grid.
TEST(SolverAR004Test, TrueResidualSquaredLocalMatchesSerialRecomputation2D) {
    const Index kNx = 3;  // 3x3 interior on a single rank with halo width 1
    const Index kNy = 3;

    Subgrid sg(kNx, kNy, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);

    // Fill the whole padded buffer (all k planes, halo included) with
    // asymmetric values: any stencil/indexing slip changes the result, and
    // non-zero halo contents prove the kernel reads halo values instead of
    // assuming zeros. Directly written halo is a valid post-exchange state,
    // so no applyPhysicalBoundary is needed for this variant.
    Real* u = sg.u().data();
    Real* rhs = sg.rhs().data();
    for (Index k = 0; k < sg.nzTotal(); ++k) {
        for (Index j = 0; j < sg.nyTotal(); ++j) {
            for (Index i = 0; i < sg.nxTotal(); ++i) {
                const Real fi = static_cast<Real>(i);
                const Real fj = static_cast<Real>(j);
                const Real fk = static_cast<Real>(k);
                const Index idx = sg.index(i, j, k);
                u[idx] = 0.30 * fi - 0.70 * fj + 0.11 * fi * fj + 0.05 * fk;
                rhs[idx] = 0.23 * fi + 0.41 * fj - 0.17 * fi * fj + 0.09 * fk;
            }
        }
    }

    const Real expected = trueResidualSumSqSerial(sg, sg.u().data());
    SCOPED_TRACE("2D: kernel vs serial recomputation");
    ASSERT_GT(expected, 0.0);  // guard: a zero expected value would void the check
    expectCloseRelative(trueResidualSquaredLocal(sg), expected,
                        "trueResidualSquaredLocal(2D)");
}

// U1 (3D): same contract on a 3x3x3 interior, exercising the realistic caller
// path where the halo comes from applyPhysicalBoundary with a non-zero
// Dirichlet value (so halo contributions are non-trivial).
TEST(SolverAR004Test, TrueResidualSquaredLocalMatchesSerialRecomputation3D) {
    const Index kNx = 3;  // 3x3x3 interior on a single rank with halo width 1
    const Index kNy = 3;
    const Index kNz = 3;
    const Real kDirichletValue = 0.25;  // non-zero so halo terms are non-trivial

    Subgrid sg(kNx, kNy, kNz, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
                    MPI_PROC_NULL, MPI_PROC_NULL);

    Real* u = sg.u().data();
    Real* rhs = sg.rhs().data();
    for (Index k = 0; k < sg.nzTotal(); ++k) {
        for (Index j = 0; j < sg.nyTotal(); ++j) {
            for (Index i = 0; i < sg.nxTotal(); ++i) {
                const Real fi = static_cast<Real>(i);
                const Real fj = static_cast<Real>(j);
                const Real fk = static_cast<Real>(k);
                const Index idx = sg.index(i, j, k);
                u[idx] = 0.13 * fi - 0.29 * fj + 0.37 * fk + 0.07 * fi * fj * fk;
                rhs[idx] = 0.19 * fi + 0.43 * fj - 0.11 * fk + 0.05 * fi * fj;
            }
        }
    }
    // Caller contract: refresh the PROC_NULL halo faces before the call.
    sg.applyPhysicalBoundary(kDirichletValue);

    const Real expected = trueResidualSumSqSerial(sg, sg.u().data());
    SCOPED_TRACE("3D: kernel vs serial recomputation");
    ASSERT_GT(expected, 0.0);
    expectCloseRelative(trueResidualSquaredLocal(sg), expected,
                        "trueResidualSquaredLocal(3D)");
}

// U2: after convergence, lastResidual() must be the true residual of the
// final iterate (confirm-scan semantics), not a diff-based proxy.
TEST(SolverAR004Test, LastResidualAfterConvergenceIsTrueResidual) {
    const Index kGridN = 32;      // manufactured sine domain 32x32
    const Real kTol = 1e-6;       // convergence threshold on the true residual
    const Index kMaxIter = 20000; // generous budget; Jacobi on 32^2 needs ~2.5k

    Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sg, kGridN);

    PointToPointExchanger ex;
    ex.initialize(sg);

    JacobiSolver solver;
    const Index iters = solver.solve(sg, ex, kMaxIter, kTol);

    SCOPED_TRACE("converged solve");
    ASSERT_GT(iters, Index(0));
    ASSERT_LT(iters, kMaxIter);  // must converge, not run out of budget
    EXPECT_LE(solver.lastResidual(), kTol);

    // Recompute the true residual of the final u. Single rank with all
    // PROC_NULL neighbors: halo IS the physical boundary, so refreshing it
    // via applyPhysicalBoundary is the correct pre-computation step.
    sg.applyPhysicalBoundary();
    const Real expected =
        std::sqrt(trueResidualSumSqSerial(sg, sg.u().data()));
    ASSERT_GT(expected, 0.0);
    expectCloseRelative(solver.lastResidual(), expected,
                        "lastResidual vs recomputed true residual");
}

// U3: iterate() must return denom * sqrt(sum(diff^2)), which by the algebraic
// identity r_t = D * diff equals the exact true residual ||r_t|| of the state
// BEFORE the update (the snapshot).
TEST(SolverAR004Test, IterateReturnsTrueResidualOfPreviousState) {
    const Index kGridN = 8; // small single-rank grid, one direct iterate() call

    Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sg, kGridN);

    // Perturb u away from zero with asymmetric values so the diff-based
    // return value is non-trivial and any indexing slip shows up.
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const Real fi = static_cast<Real>(i);
            const Real fj = static_cast<Real>(j);
            sg.at(i, j) = 0.10 * fi - 0.05 * fj + 0.02 * fi * fj;
        }
    }
    // Refresh the PROC_NULL halo faces: iterate() consumes the halo as-is.
    sg.applyPhysicalBoundary();

    // Snapshot of u before the update; iterate() swaps buffers afterwards,
    // so the pre-iteration state must be captured now.
    const Real* uBefore = sg.u().data();
    const std::vector<Real> snapshot(uBefore, uBefore + sg.totalCells());

    PointToPointExchanger ex;
    ex.initialize(sg);

    JacobiSolver solver;
    const Real returned = solver.iterate(sg, ex);

    const Real expected =
        std::sqrt(trueResidualSumSqSerial(sg, snapshot.data()));
    SCOPED_TRACE("iterate() return vs true residual of the snapshot");
    ASSERT_GT(expected, 0.0);
    expectCloseRelative(returned, expected,
                        "iterate() return vs ||r_t|| of pre-update state");
}

// U9: setResidualCheckInterval must clamp values < 1 to 1. The clamp cannot
// be observed through the private member, so we assert its behavioral
// consequence: with 0 (and a negative int, which wraps to a huge unsigned
// Index), the residual_allreduce count stays at the every-iteration default.
TEST(SolverAR004Test, ResidualCheckIntervalClampedToMinimumOne) {
    const Index kNx = 8;
    const Index kNy = 8;
    const Index kIters = 6;   // fixed iteration count, never converges below tol 0
    const Real kNoConvergeTol = 0.0; // residual >= 0 never drops below 0

    // Each run gets a fresh subgrid/exchanger so runs stay independent.
    const auto runFixedIters = [&](JacobiSolver& solver) {
        Subgrid sg(kNx, kNy, 1, 1, MPI_COMM_SELF);
        sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
        setupUniformNegativeRhs(sg);
        PointToPointExchanger ex;
        ex.initialize(sg);
        solver.solve(sg, ex, kIters, kNoConvergeTol);
    };

    // Baseline: default interval (1) checks after every iteration.
    Profiler::instance().reset();
    {
        JacobiSolver solver;
        runFixedIters(solver);
    }
    const auto baseline =
        Profiler::instance().stats("residual_allreduce").callCount;
    EXPECT_EQ(baseline, std::uint64_t{kIters});
    EXPECT_EQ(Profiler::instance().stats("jacobi_iteration").callCount,
              std::uint64_t{kIters});

    // interval = 0 must be clamped to 1: same allreduce count as baseline.
    Profiler::instance().reset();
    {
        JacobiSolver solver;
        solver.setResidualCheckInterval(0);
        runFixedIters(solver);
    }
    SCOPED_TRACE("setResidualCheckInterval(0)");
    EXPECT_EQ(Profiler::instance().stats("residual_allreduce").callCount, baseline);
    EXPECT_EQ(Profiler::instance().stats("jacobi_iteration").callCount,
              std::uint64_t{kIters});

    // A negative int becomes a huge unsigned Index; still clamped to 1.
    Profiler::instance().reset();
    {
        JacobiSolver solver;
        solver.setResidualCheckInterval(static_cast<Index>(-5));
        runFixedIters(solver);
    }
    SCOPED_TRACE("setResidualCheckInterval(-5)");
    EXPECT_EQ(Profiler::instance().stats("residual_allreduce").callCount, baseline);
    EXPECT_EQ(Profiler::instance().stats("jacobi_iteration").callCount,
              std::uint64_t{kIters});
}

// U10: when the loop exits at maxIter, exactly one confirmation scan runs and
// lastResidual() is the true residual of the final iterate.
TEST(SolverAR004Test, MaxIterExitRunsSingleConfirmationScan) {
    const Index kGridN = 32;      // manufactured sine domain
    const Index kMaxIter = 5;     // exhaust the budget
    const Real kNoConvergeTol = 0.0;

    Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sg, kGridN);

    PointToPointExchanger ex;
    ex.initialize(sg);

    Profiler::instance().reset();
    JacobiSolver solver;
    const Index iters = solver.solve(sg, ex, kMaxIter, kNoConvergeTol);

    SCOPED_TRACE("maxIter exit");
    EXPECT_EQ(iters, kMaxIter);
    EXPECT_EQ(Profiler::instance().stats("residual_confirm").callCount,
              std::uint64_t{1});

    sg.applyPhysicalBoundary();
    const Real expected =
        std::sqrt(trueResidualSumSqSerial(sg, sg.u().data()));
    ASSERT_GT(expected, 0.0);
    expectCloseRelative(solver.lastResidual(), expected,
                        "lastResidual vs recomputed true residual at exit");
}

// I1s (solver level): with interval k, the in-loop Allreduce check fires only
// every k-th iteration: maxIter=20, interval=10 -> 2 checks (at 10 and 20).
TEST(SolverAR004Test, IntervalTenHalvesAllreduceChecks) {
    const Index kGridN = 32;   // manufactured sine domain
    const Index kMaxIter = 20;
    const Index kInterval = 10;      // floor(20/10) = 2 in-loop checks
    const Real kNoConvergeTol = 0.0;

    Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sg, kGridN);

    PointToPointExchanger ex;
    ex.initialize(sg);

    Profiler::instance().reset();
    JacobiSolver solver;
    solver.setResidualCheckInterval(kInterval);
    solver.solve(sg, ex, kMaxIter, kNoConvergeTol);

    SCOPED_TRACE("interval=10, maxIter=20");
    EXPECT_EQ(Profiler::instance().stats("residual_allreduce").callCount,
              std::uint64_t{2});
    EXPECT_EQ(Profiler::instance().stats("residual_confirm").callCount,
              std::uint64_t{1});
    EXPECT_EQ(Profiler::instance().stats("jacobi_iteration").callCount,
              std::uint64_t{kMaxIter});
}

// I5s: default interval (k=1) keeps the every-iteration check.
TEST(SolverAR004Test, DefaultIntervalChecksEveryIteration) {
    const Index kGridN = 32;
    const Index kMaxIter = 20;
    const Real kNoConvergeTol = 0.0;

    Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sg, kGridN);

    PointToPointExchanger ex;
    ex.initialize(sg);

    Profiler::instance().reset();
    JacobiSolver solver;  // no setResidualCheckInterval: default must be 1
    solver.solve(sg, ex, kMaxIter, kNoConvergeTol);

    SCOPED_TRACE("default interval, maxIter=20");
    EXPECT_EQ(Profiler::instance().stats("residual_allreduce").callCount,
              std::uint64_t{kMaxIter});
    EXPECT_EQ(Profiler::instance().stats("residual_confirm").callCount,
              std::uint64_t{1});
    EXPECT_EQ(Profiler::instance().stats("jacobi_iteration").callCount,
              std::uint64_t{kMaxIter});
}

// E4s: interval >= maxIter means no in-loop check at all; the exit still gets
// exactly one confirmation scan and a true lastResidual.
TEST(SolverAR004Test, IntervalAboveMaxIterSkipsInLoopChecks) {
    const Index kGridN = 32;
    const Index kMaxIter = 5;
    const Index kInterval = 100;    // no completed iteration is a multiple of 100
    const Real kNoConvergeTol = 0.0;

    Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sg, kGridN);

    PointToPointExchanger ex;
    ex.initialize(sg);

    Profiler::instance().reset();
    JacobiSolver solver;
    solver.setResidualCheckInterval(kInterval);
    const Index iters = solver.solve(sg, ex, kMaxIter, kNoConvergeTol);

    SCOPED_TRACE("interval=100, maxIter=5");
    EXPECT_EQ(iters, kMaxIter);
    EXPECT_EQ(Profiler::instance().stats("residual_allreduce").callCount,
              std::uint64_t{0});
    EXPECT_EQ(Profiler::instance().stats("residual_confirm").callCount,
              std::uint64_t{1});
    EXPECT_EQ(Profiler::instance().stats("jacobi_iteration").callCount,
              std::uint64_t{kMaxIter});

    sg.applyPhysicalBoundary();
    const Real expected =
        std::sqrt(trueResidualSumSqSerial(sg, sg.u().data()));
    ASSERT_GT(expected, 0.0);
    expectCloseRelative(solver.lastResidual(), expected,
                        "lastResidual vs recomputed true residual at exit");
}
