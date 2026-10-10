#include <gtest/gtest.h>
#include <mpi.h>

#include "comm/halo_exchanger.hpp"
#include "grid/subgrid.hpp"
#include "solver/solver.hpp"
#include "solver/residual.hpp"
#include "solver/mg_operators.hpp"

#include <cmath>
#include <vector>

using namespace hypo;

namespace {

constexpr Index kHalo = 1;

// AR007 design §4.2: manual non-uniform 2x2 layout {15,19} on 34 global
// (34 = 17+17 would be uniform; 15/19 exercises odd ends and all corners).
constexpr long long kLayoutN = 34;
Index layoutOffset(int coord) { return coord == 0 ? 0 : 15; }
Index layoutSize(int coord) { return coord == 0 ? 15 : 19; }

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

// Analytic right-hand side of -Laplace(u) = 2*pi^2*sin(pi*x)*sin(pi*y) with
// u = 0 on the boundary, discretized as (4u - sum(nb)) = f*h^2. Unlike
// fillSineRhs (the discrete-consistent rhs, whose exact discrete solution is
// the sine samples themselves and therefore has ZERO discretization error),
// this rhs leaves a genuine O(h^2) truncation error for the convergence-order
// test to measure. Solver convention: A*u = -rhs, hence rhs = -f*h^2.
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

} // namespace

// ---------------------------------------------------------------------------
// TwoLevelMGUnitTest (np = 1)
// ---------------------------------------------------------------------------

TEST(TwoLevelMGUnitTest, CoarseRangePureFunction) {
    // design §4.1 hand-computed table: even/odd offsets, odd nLocal, empty.
    EXPECT_EQ(coarseRangeBegin(0, 4), Index(0));
    EXPECT_EQ(coarseRangeCount(0, 4), Index(2));
    EXPECT_EQ(coarseRangeBegin(3, 3), Index(2));
    EXPECT_EQ(coarseRangeCount(3, 3), Index(1));
    EXPECT_EQ(coarseRangeBegin(1, 1), Index(1));
    EXPECT_EQ(coarseRangeCount(1, 1), Index(0)); // degenerate empty piece
    EXPECT_EQ(coarseRangeBegin(0, 15), Index(0));
    EXPECT_EQ(coarseRangeCount(0, 15), Index(8));
    EXPECT_EQ(coarseRangeBegin(15, 19), Index(8));
    EXPECT_EQ(coarseRangeCount(15, 19), Index(9));
}

TEST(TwoLevelMGUnitTest, RestrictionKnownCoefficients) {
    const Index nx = 8;
    const Index nxH = 4;
    Subgrid fine(nx, nx, 1, kHalo, MPI_COMM_SELF);
    AlignedBuffer<Real> buf(fine.totalCells());
    buf.fill(0.0);

    // constant field 1.0 in the interior; physical halo -> 0 (outside).
    for (Index j = fine.jBegin(); j < fine.jEnd(); ++j) {
        for (Index i = fine.iBegin(); i < fine.iEnd(); ++i) {
            buf[fine.index(i, j)] = 1.0;
        }
    }
    fine.applyPhysicalBoundary(buf.data(), 0.0);

    AlignedBuffer<Real> packed(nxH * nxH);
    restrictResidual(fine, buf.data(), 0, 0, nxH, nxH, packed.data());

    // vertex-coincident FW: only I=0 / J=0 coarse lines see the outside
    // (right/top fine points nx-1 stay inside the 9-point stencil).
    const Real corner = 9.0 / 16.0;
    const Real edge = 12.0 / 16.0;
    for (Index J = 0; J < nxH; ++J) {
        for (Index I = 0; I < nxH; ++I) {
            const Real expected =
                (I == 0 && J == 0) ? corner : (I == 0 || J == 0) ? edge : 1.0;
            EXPECT_DOUBLE_EQ(packed[J * nxH + I], expected)
                << "I=" << I << " J=" << J;
        }
    }

    // linear field r(gi,gj) = gi + 2*gj against an in-test reference
    // (outside = 0, weights 4/2/1 over 16).
    for (Index j = fine.jBegin(); j < fine.jEnd(); ++j) {
        for (Index i = fine.iBegin(); i < fine.iEnd(); ++i) {
            const long long gi = static_cast<long long>(i) - kHalo;
            const long long gj = static_cast<long long>(j) - kHalo;
            buf[fine.index(i, j)] = static_cast<Real>(gi + 2 * gj);
        }
    }
    fine.applyPhysicalBoundary(buf.data(), 0.0);
    restrictResidual(fine, buf.data(), 0, 0, nxH, nxH, packed.data());

    auto rAt = [&](long long gi, long long gj) -> Real {
        if (gi < 0 || gj < 0 || gi >= static_cast<long long>(nx) ||
            gj >= static_cast<long long>(nx)) {
            return 0.0;
        }
        return static_cast<Real>(gi + 2 * gj);
    };
    for (Index J = 0; J < nxH; ++J) {
        for (Index I = 0; I < nxH; ++I) {
            const long long ci = 2 * static_cast<long long>(I);
            const long long cj = 2 * static_cast<long long>(J);
            const Real expected =
                (4.0 * rAt(ci, cj) +
                 2.0 * (rAt(ci - 1, cj) + rAt(ci + 1, cj) + rAt(ci, cj - 1) + rAt(ci, cj + 1)) +
                 rAt(ci - 1, cj - 1) + rAt(ci - 1, cj + 1) + rAt(ci + 1, cj - 1) +
                 rAt(ci + 1, cj + 1)) / 16.0;
            EXPECT_DOUBLE_EQ(packed[J * nxH + I], expected)
                << "I=" << I << " J=" << J;
        }
    }
}

TEST(TwoLevelMGUnitTest, ProlongationKnownCoefficients) {
    const Index nx = 8;
    const Index nxH = 4;
    Subgrid fine(nx, nx, 1, kHalo, MPI_COMM_SELF);
    Subgrid coarse(nxH, nxH, 1, kHalo, MPI_COMM_SELF);
    fine.zeroInitialize();

    // constant coarse field: even fine points and odd fine points with both
    // parents inside keep the constant; fine nx-1 (odd, parent+1 outside)
    // gets the recorded 1/2-weight asymmetry (design D3).
    for (Index j = coarse.jBegin(); j < coarse.jEnd(); ++j) {
        for (Index i = coarse.iBegin(); i < coarse.iEnd(); ++i) {
            coarse.u().data()[coarse.index(i, j)] = 5.0;
        }
    }
    prolongateCorrection(fine, coarse);
    for (Index j = fine.jBegin(); j < fine.jEnd(); ++j) {
        for (Index i = fine.iBegin(); i < fine.iEnd(); ++i) {
            const long long gi = static_cast<long long>(i) - kHalo;
            const long long gj = static_cast<long long>(j) - kHalo;
            Real expected = 5.0;
            if (gi % 2 == 1 && gi == static_cast<long long>(nx) - 1) expected = 2.5;
            if (gj % 2 == 1 && gj == static_cast<long long>(nx) - 1) expected *= 0.5;
            EXPECT_DOUBLE_EQ(fine.u().data()[fine.index(i, j)], expected)
                << "gi=" << gi << " gj=" << gj;
        }
    }

    // unit impulse at coarse (2,1): tensor 1 / 1/2 / 1/4 weights.
    fine.zeroInitialize();
    for (Index k = 0; k < coarse.totalCells(); ++k) {
        coarse.u().data()[k] = 0.0;
    }
    coarse.u().data()[coarse.index(2 + kHalo, 1 + kHalo)] = 1.0;
    prolongateCorrection(fine, coarse);
    EXPECT_DOUBLE_EQ(fine.u().data()[fine.index(4 + kHalo, 2 + kHalo)], 1.0);
    EXPECT_DOUBLE_EQ(fine.u().data()[fine.index(3 + kHalo, 2 + kHalo)], 0.5);
    EXPECT_DOUBLE_EQ(fine.u().data()[fine.index(4 + kHalo, 1 + kHalo)], 0.5);
    EXPECT_DOUBLE_EQ(fine.u().data()[fine.index(3 + kHalo, 1 + kHalo)], 0.25);
    EXPECT_DOUBLE_EQ(fine.u().data()[fine.index(5 + kHalo, 3 + kHalo)], 0.25);
    EXPECT_DOUBLE_EQ(fine.u().data()[fine.index(0 + kHalo, 0 + kHalo)], 0.0);
}

TEST(TwoLevelMGUnitTest, ProlongRestrictProjectionIdentity) {
    const Index nx = 8;
    const Index nxH = 4;
    Subgrid fine(nx, nx, 1, kHalo, MPI_COMM_SELF);
    Subgrid coarse(nxH, nxH, 1, kHalo, MPI_COMM_SELF);
    AlignedBuffer<Real> buf(fine.totalCells());
    buf.fill(0.0);

    for (Index j = fine.jBegin(); j < fine.jEnd(); ++j) {
        for (Index i = fine.iBegin(); i < fine.iEnd(); ++i) {
            const long long gi = static_cast<long long>(i) - kHalo;
            const long long gj = static_cast<long long>(j) - kHalo;
            buf[fine.index(i, j)] =
                static_cast<Real>(gi * gi + 3 * gj + 1) / 7.0;
        }
    }
    fine.applyPhysicalBoundary(buf.data(), 0.0);

    AlignedBuffer<Real> packed(nxH * nxH);
    restrictResidual(fine, buf.data(), 0, 0, nxH, nxH, packed.data());

    coarse.zeroInitialize();
    for (Index j = coarse.jBegin(); j < coarse.jEnd(); ++j) {
        for (Index i = coarse.iBegin(); i < coarse.iEnd(); ++i) {
            coarse.u().data()[coarse.index(i, j)] =
                packed[(j - kHalo) * nxH + (i - kHalo)];
        }
    }
    fine.zeroInitialize();
    prolongateCorrection(fine, coarse);

    // strict identity at the coarse-coincident (even, even) fine points.
    for (Index J = 0; J < nxH; ++J) {
        for (Index I = 0; I < nxH; ++I) {
            EXPECT_DOUBLE_EQ(fine.u().data()[fine.index(2 * I + kHalo, 2 * J + kHalo)],
                             packed[J * nxH + I])
                << "I=" << I << " J=" << J;
        }
    }
}

TEST(TwoLevelMGUnitTest, ResidualFieldMatchesTrueResidualScan) {
    Subgrid sg(8, 8, 1, kHalo, MPI_COMM_SELF);
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gi = static_cast<long long>(i) - kHalo;
            const long long gj = static_cast<long long>(j) - kHalo;
            sg.u().data()[sg.index(i, j)] = std::cos(0.3 * gi + 0.7 * gj);
            sg.rhs().data()[sg.index(i, j)] = 0.05 * (gi - gj);
        }
    }
    sg.applyPhysicalBoundary();

    AlignedBuffer<Real> field(sg.totalCells());
    field.fill(0.0);
    residualFieldLocal(sg, field.data());

    Real fieldSq = 0.0;
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            fieldSq += field[sg.index(i, j)] * field[sg.index(i, j)];
        }
    }
    const Real scanSq = trueResidualSquaredLocal(sg);
    EXPECT_GT(scanSq, 0.0);
    EXPECT_NEAR(std::sqrt(fieldSq), std::sqrt(scanSq),
                1e-15 * std::sqrt(scanSq) + 1e-300);
}

TEST(TwoLevelMGUnitTest, SmoothMatchesIterateBitExact) {
    auto init = [](Subgrid& sg) {
        sg.zeroInitialize();
        sg.applyDirichletBC(0.0);
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                const long long gi = static_cast<long long>(i) - kHalo;
                const long long gj = static_cast<long long>(j) - kHalo;
                sg.u().data()[sg.index(i, j)] = 0.01 * (gi * gj % 13) - 0.06;
                sg.rhs().data()[sg.index(i, j)] = -0.02 * ((gi + gj) % 7);
            }
        }
        sg.applyPhysicalBoundary();
    };

    Subgrid sgA(16, 16, 1, kHalo, MPI_COMM_SELF);
    Subgrid sgB(16, 16, 1, kHalo, MPI_COMM_SELF);
    init(sgA);
    init(sgB);

    PointToPointExchanger exA;
    exA.initialize(sgA);
    PointToPointExchanger exB;
    exB.initialize(sgB);

    RedBlackGSSolver smoother;
    smoother.smooth(sgA, exA, 3);

    RedBlackGSSolver viaIterate;
    for (Index k = 0; k < 3; ++k) {
        viaIterate.iterate(sgB, exB);
    }

    for (Index idx = 0; idx < sgA.totalCells(); ++idx) {
        ASSERT_EQ(sgA.u().data()[idx], sgB.u().data()[idx]) << "idx=" << idx;
    }
}

namespace {

// design D4 defense contract: HYPOS_ERROR log + 0 iterations + u untouched.
void expectRejected(Subgrid& sg, HaloExchanger& ex) {
    for (Index idx = 0; idx < sg.totalCells(); ++idx) {
        sg.u().data()[idx] = 7.5;
    }
    TwoLevelMGSolver solver;
    const Index iters = solver.solve(sg, ex, 10, 1e-6);
    EXPECT_EQ(iters, Index(0));
    EXPECT_EQ(solver.lastResidual(), 0.0);
    for (Index idx = 0; idx < sg.totalCells(); ++idx) {
        ASSERT_EQ(sg.u().data()[idx], 7.5) << "idx=" << idx;
    }
}

} // namespace

TEST(TwoLevelMGUnitTest, DefenseRejectsUnsupportedConfigs) {
    {
        Subgrid sg(8, 8, 2, kHalo, MPI_COMM_SELF); // global 3D
        PointToPointExchanger ex;
        ex.initialize(sg);
        expectRejected(sg, ex);
    }
    {
        Subgrid sg(5, 5, 1, kHalo, MPI_COMM_SELF); // odd global nx
        PointToPointExchanger ex;
        ex.initialize(sg);
        expectRejected(sg, ex);
    }
    {
        Subgrid sg(8, 8, 1, kHalo, MPI_COMM_SELF);
        sg.setBoundaryCondition(BoundaryCondition::Neumann);
        PointToPointExchanger ex;
        ex.initialize(sg);
        expectRejected(sg, ex);
    }
    {
        Subgrid sg(6, 6, 1, kHalo, MPI_COMM_SELF); // coarse 3x3 < 4x4
        PointToPointExchanger ex;
        ex.initialize(sg);
        expectRejected(sg, ex);
    }
    {
        // A single COMM_SELF subgrid that does not tile its own coarse grid:
        // offset (2, 0) makes the inferred global grid 10 x 8 (coarse 5 x 4),
        // but this subgrid's coarse counts sum to 4 x 4 = 16 != 20 cells.
        Subgrid sg(8, 8, 1, kHalo, MPI_COMM_SELF);
        sg.setOffsets(2, 0, 0);
        PointToPointExchanger ex;
        ex.initialize(sg);
        expectRejected(sg, ex);
    }
}

TEST(TwoLevelMGUnitTest, ManufacturedSineConvergenceNp1) {
    const long long n = 64;
    Subgrid sg(n, n, 1, kHalo, MPI_COMM_SELF);
    fillSineRhs(sg, n, n);

    PointToPointExchanger ex;
    ex.initialize(sg);
    TwoLevelMGSolver solver;
    const Index iters = solver.solve(sg, ex, 500, 1e-7);
    EXPECT_LT(iters, Index(500));
    EXPECT_GT(solver.lastResidual(), 0.0);
    EXPECT_LE(solver.lastResidual(), 1e-7);

    Real localMax = 0.0;
    const Real sum = sineL2ErrorLocal(sg, n, n, localMax);
    const Real l2 = std::sqrt(sum / static_cast<Real>(n * n));
    EXPECT_LT(l2, 1e-3);
}

TEST(TwoLevelMGUnitTest, SecondOrderConvergence) {
    // Analytic rhs: the L2 solution error is dominated by the O(h^2)
    // truncation error, so refining 32 -> 64 must shrink it by ~4x.
    auto solveError = [](long long n) -> Real {
        Subgrid sg(n, n, 1, kHalo, MPI_COMM_SELF);
        fillAnalyticSineRhs(sg, n);
        PointToPointExchanger ex;
        ex.initialize(sg);
        TwoLevelMGSolver solver;
        const Index iters = solver.solve(sg, ex, 1000, 1e-8);
        EXPECT_LT(iters, Index(1000));
        Real localMax = 0.0;
        const Real sum = sineL2ErrorLocal(sg, n, n, localMax);
        return std::sqrt(sum / static_cast<Real>(n * n));
    };

    const Real err32 = solveError(32);
    const Real err64 = solveError(64);
    EXPECT_GT(err32, 0.0);
    EXPECT_GT(err64, 0.0);
    const Real ratio = err32 / err64;
    EXPECT_GT(ratio, 3.5);
    EXPECT_LT(ratio, 4.5) << "err32=" << err32 << " err64=" << err64;
}

// ---------------------------------------------------------------------------
// TwoLevelMGMpiTest (np = 4)
// ---------------------------------------------------------------------------

TEST(TwoLevelMGMpiTest, ManufacturedSineConvergenceNp4) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "ManufacturedSineConvergenceNp4 requires exactly 4 processes";
    }

    const long long n = 64;
    int dims[2] = {2, 2};
    int periods[2] = {0, 0};
    MPI_Comm cart = MPI_COMM_NULL;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart);

    int left = MPI_PROC_NULL, right = MPI_PROC_NULL;
    int down = MPI_PROC_NULL, up = MPI_PROC_NULL;
    MPI_Cart_shift(cart, 0, 1, &left, &right);
    MPI_Cart_shift(cart, 1, 1, &down, &up);
    int coords[2] = {0, 0};
    MPI_Cart_coords(cart, rank, 2, coords);

    Subgrid sg(32, 32, 1, kHalo, cart);
    sg.setNeighbors(left, right, down, up);
    sg.setOffsets(static_cast<Index>(coords[0]) * 32,
                  static_cast<Index>(coords[1]) * 32, 0);
    fillSineRhs(sg, n, n);

    PointToPointExchanger ex;
    ex.initialize(sg);
    TwoLevelMGSolver solver;
    const Index iters = solver.solve(sg, ex, 200, 1e-7);
    EXPECT_LT(iters, Index(200));
    EXPECT_LE(solver.lastResidual(), 1e-7);

    Real localMax = 0.0;
    Real local = sineL2ErrorLocal(sg, n, n, localMax);
    Real global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, cart);
    const Real l2 = std::sqrt(global / static_cast<Real>(n * n));
    EXPECT_LT(l2, 1e-3);

    MPI_Comm_free(&cart);
}

namespace {

// Shared scaffolding for the manual non-uniform {15,19} 2x2 layout.
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

} // namespace

TEST(TwoLevelMGMpiTest, CoarseOwnershipPartitionsGlobalRange) {
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "CoarseOwnershipPartitionsGlobalRange requires exactly 4 processes";
    }
    NonUniformLayout lay;
    lay.build();
    if (HasFailure()) return;

    // hand-computed anchors for {15,19}: [0,8) count 8 and [8,17) count 9.
    const Index bx = coarseRangeBegin(lay.sg->offsetX(), lay.sg->nxLocal());
    const Index cx = coarseRangeCount(lay.sg->offsetX(), lay.sg->nxLocal());
    const Index by = coarseRangeBegin(lay.sg->offsetY(), lay.sg->nyLocal());
    const Index cy = coarseRangeCount(lay.sg->offsetY(), lay.sg->nyLocal());
    EXPECT_EQ(bx, lay.coords[0] == 0 ? Index(0) : Index(8));
    EXPECT_EQ(cx, lay.coords[0] == 0 ? Index(8) : Index(9));
    EXPECT_EQ(by, lay.coords[1] == 0 ? Index(0) : Index(8));
    EXPECT_EQ(cy, lay.coords[1] == 0 ? Index(8) : Index(9));

    // allgather (begin, count, column) and prove the DISTINCT column ranges
    // tile [0, 17) seam-free (ranks sharing a column must agree).
    int r12[12] = {};
    const int mine3[3] = {static_cast<int>(bx), static_cast<int>(cx),
                          lay.coords[0]};
    MPI_Allgather(mine3, 3, MPI_INT, r12, 3, MPI_INT, lay.cart);
    int distinct[2][2] = {{-1, -1}, {-1, -1}};
    for (int r = 0; r < 4; ++r) {
        const int col = r12[3 * r + 2];
        ASSERT_TRUE(col == 0 || col == 1);
        if (distinct[col][0] < 0) {
            distinct[col][0] = r12[3 * r + 0];
            distinct[col][1] = r12[3 * r + 1];
        }
        EXPECT_EQ(distinct[col][0], r12[3 * r + 0]);
        EXPECT_EQ(distinct[col][1], r12[3 * r + 1]);
    }
    bool covered[17] = {};
    Index total = 0;
    for (int col = 0; col < 2; ++col) {
        for (int v = distinct[col][0]; v < distinct[col][0] + distinct[col][1]; ++v) {
            ASSERT_FALSE(covered[v]) << "overlap at " << v;
            covered[v] = true;
        }
        total += static_cast<Index>(distinct[col][1]);
    }
    for (int v = 0; v < 17; ++v) {
        EXPECT_TRUE(covered[v]) << "hole at " << v;
    }
    EXPECT_EQ(total, Index(17));

    // degenerate empty piece stays a legal Allgatherv count.
    EXPECT_EQ(coarseRangeCount(1, 1), Index(0));

    lay.destroy();
}

TEST(TwoLevelMGMpiTest, DiagonalCornerExchangeFilled) {
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "DiagonalCornerExchangeFilled requires exactly 4 processes";
    }
    NonUniformLayout lay;
    lay.build();
    if (HasFailure()) return;
    Subgrid& sg = *lay.sg;

    // whole-field buffer with globally unique values; physical halo -> 0.
    AlignedBuffer<Real> buf(sg.totalCells());
    buf.fill(0.0);
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gi = static_cast<long long>(sg.offsetX()) +
                                 static_cast<long long>(i) - kHalo;
            const long long gj = static_cast<long long>(sg.offsetY()) +
                                 static_cast<long long>(j) - kHalo;
            buf[sg.index(i, j)] = static_cast<Real>((gi + 1) * 1000 + (gj + 1));
        }
    }
    sg.applyPhysicalBoundary(buf.data(), 0.0);

    exchangeCorners(sg, buf.data());

    const long long nxL = static_cast<long long>(sg.nxLocal());
    const long long nyL = static_cast<long long>(sg.nyLocal());
    const long long ox = static_cast<long long>(sg.offsetX());
    const long long oy = static_cast<long long>(sg.offsetY());
    const bool hasXm = sg.neighborLeft() != MPI_PROC_NULL;
    const bool hasXp = sg.neighborRight() != MPI_PROC_NULL;
    const bool hasYm = sg.neighborDown() != MPI_PROC_NULL;
    const bool hasYp = sg.neighborUp() != MPI_PROC_NULL;

    struct Corner { bool has; long long gi; long long gj; Index li; Index lj; };
    const Corner corners[4] = {
        {hasXm && hasYm, ox - 1, oy - 1, sg.iBegin() - 1, sg.jBegin() - 1},
        {hasXp && hasYm, ox + nxL, oy - 1, sg.iEnd(), sg.jBegin() - 1},
        {hasXm && hasYp, ox - 1, oy + nyL, sg.iBegin() - 1, sg.jEnd()},
        {hasXp && hasYp, ox + nxL, oy + nyL, sg.iEnd(), sg.jEnd()},
    };
    for (const Corner& c : corners) {
        const Real expected =
            c.has ? static_cast<Real>((c.gi + 1) * 1000 + (c.gj + 1)) : 0.0;
        EXPECT_DOUBLE_EQ(buf[sg.index(c.li, c.lj)], expected)
            << "corner at global (" << c.gi << "," << c.gj << ")";
    }

    lay.destroy();
}

TEST(TwoLevelMGMpiTest, CrossLayoutCorrectionConsistency) {
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "CrossLayoutCorrectionConsistency requires exactly 4 processes";
    }
    NonUniformLayout lay;
    lay.build();
    if (HasFailure()) return;
    Subgrid& sg = *lay.sg;
    fillSineRhs(sg, kLayoutN, kLayoutN);

    TwoLevelMGSolver solver;
    const Index iters = solver.solve(sg, *lay.ex, 1, 0.0);
    EXPECT_EQ(iters, Index(1));

    // Serial reference recomputed redundantly on every rank (34^2 is tiny).
    Subgrid ref(kLayoutN, kLayoutN, 1, kHalo, MPI_COMM_SELF);
    fillSineRhs(ref, kLayoutN, kLayoutN);
    PointToPointExchanger refEx;
    refEx.initialize(ref);
    TwoLevelMGSolver refSolver;
    refSolver.solve(ref, refEx, 1, 0.0);

    Real localMax = 0.0;
    Real maxDiff = 0.0;
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gi = static_cast<long long>(sg.offsetX()) +
                                 static_cast<long long>(i) - kHalo;
            const long long gj = static_cast<long long>(sg.offsetY()) +
                                 static_cast<long long>(j) - kHalo;
            const Real got = sg.u().data()[sg.index(i, j)];
            const Real want = ref.u().data()[ref.index(
                static_cast<Index>(gi + kHalo), static_cast<Index>(gj + kHalo))];
            maxDiff = std::max(maxDiff, std::fabs(got - want));
            localMax = std::max(localMax, std::fabs(want));
        }
    }
    Real globalMax[2] = {maxDiff, localMax};
    Real allMax[2] = {0.0, 0.0};
    MPI_Allreduce(globalMax, allMax, 2, MPI_DOUBLE, MPI_MAX, lay.cart);

    ASSERT_GT(allMax[1], 1e-6); // non-degenerate correction field
    EXPECT_LE(allMax[0], 1e-12 * allMax[1])
        << "maxDiff=" << allMax[0] << " scale=" << allMax[1];

    lay.destroy();
}
