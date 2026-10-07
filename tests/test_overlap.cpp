#include <gtest/gtest.h>
#include <mpi.h>

#include "comm/halo_exchanger.hpp"
#include "grid/subgrid.hpp"
#include "solver/solver.hpp"

using namespace hypo;

namespace {

void setupPoisson(Subgrid& sg) {
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);

    Real* rhs = sg.rhs().data();
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            rhs[sg.index(i, j)] = -1.0;
        }
    }
}

void expectSameSolution(const Subgrid& a, const Subgrid& b, Real tolerance) {
    const Index n = a.totalCells();
    ASSERT_EQ(b.totalCells(), n);
    const Real* ua = a.u().data();
    const Real* ub = b.u().data();
    for (Index idx = 0; idx < n; ++idx) {
        EXPECT_NEAR(ua[idx], ub[idx], tolerance) << "cell index " << idx;
    }
}

} // namespace

TEST(OverlapTest, SelfLoopOnOffConsistency) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Subgrid sgOff(16, 16, 1, 1, MPI_COMM_WORLD);
    Subgrid sgOn(16, 16, 1, 1, MPI_COMM_WORLD);
    sgOff.setNeighbors(rank, rank, rank, rank);
    sgOn.setNeighbors(rank, rank, rank, rank);

    setupPoisson(sgOff);
    setupPoisson(sgOn);

    PointToPointExchanger exOff;
    PointToPointExchanger exOn;
    exOff.initialize(sgOff);
    exOn.initialize(sgOn);

    JacobiSolver sOff(false);
    JacobiSolver sOn(true);

    Index itOff = sOff.solve(sgOff, exOff, 5000, 1e-6);
    Index itOn = sOn.solve(sgOn, exOn, 5000, 1e-6);

    EXPECT_EQ(itOff, itOn);
    expectSameSolution(sgOff, sgOn, 1e-12);
    EXPECT_NEAR(sOff.lastResidual(), sOn.lastResidual(), 1e-12);
}

TEST(OverlapTest, MpiFourRanksOnOffConsistency) {
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "MpiFourRanksOnOffConsistency requires exactly 4 processes";
    }

    int dims[2] = {2, 2};
    int periods[2] = {0, 0};
    MPI_Comm cart = MPI_COMM_NULL;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart);

    int left = MPI_PROC_NULL;
    int right = MPI_PROC_NULL;
    int down = MPI_PROC_NULL;
    int up = MPI_PROC_NULL;
    MPI_Cart_shift(cart, 0, 1, &left, &right);
    MPI_Cart_shift(cart, 1, 1, &down, &up);

    // Global 32x32 over a 2x2 rank grid: each rank owns a 16x16 block.
    Subgrid sgOff(16, 16, 1, 1, cart);
    Subgrid sgOn(16, 16, 1, 1, cart);
    sgOff.setNeighbors(left, right, down, up);
    sgOn.setNeighbors(left, right, down, up);

    setupPoisson(sgOff);
    setupPoisson(sgOn);

    PointToPointExchanger exOff;
    PointToPointExchanger exOn;
    exOff.initialize(sgOff);
    exOn.initialize(sgOn);

    JacobiSolver sOff(false);
    JacobiSolver sOn(true);

    Index itOff = sOff.solve(sgOff, exOff, 5000, 1e-6);
    Index itOn = sOn.solve(sgOn, exOn, 5000, 1e-6);

    EXPECT_EQ(itOff, itOn);
    expectSameSolution(sgOff, sgOn, 1e-12);
    EXPECT_NEAR(sOff.lastResidual(), sOn.lastResidual(), 1e-12);

    MPI_Comm_free(&cart);
}
