#include <gtest/gtest.h>
#include <mpi.h>

#include "comm/halo_exchanger.hpp"
#include "grid/subgrid.hpp"
#include "solver/solver.hpp"

#include <cmath>

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

void setupPoisson3D(Subgrid& sg) {
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);

    Real* rhs = sg.rhs().data();
    for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                rhs[sg.index(i, j, k)] = -1.0;
            }
        }
    }
}

} // namespace

TEST(SolverMpiTest, ConvergedSolutionMatchesManufacturedSine) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "ConvergedSolutionMatchesManufacturedSine requires exactly 4 processes";
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

    int coords[2] = {0, 0};
    MPI_Cart_coords(cart, rank, 2, coords);
    const long long offsetX = static_cast<long long>(coords[0]) * 32;
    const long long offsetY = static_cast<long long>(coords[1]) * 32;

    Subgrid sg(32, 32, 1, 1, cart);
    sg.setNeighbors(left, right, down, up);
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);

    const Real pi = std::acos(-1.0);
    auto exact = [pi](long long gI, long long gJ) {
        return std::sin(pi * static_cast<Real>(gI + 1) / 65.0) *
               std::sin(pi * static_cast<Real>(gJ + 1) / 65.0);
    };

    Real* rhs = sg.rhs().data();
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gI = offsetX + static_cast<long long>(i) - 1;
            const long long gJ = offsetY + static_cast<long long>(j) - 1;
            rhs[sg.index(i, j)] = exact(gI - 1, gJ) + exact(gI + 1, gJ) +
                                  exact(gI, gJ - 1) + exact(gI, gJ + 1) -
                                  4.0 * exact(gI, gJ);
        }
    }

    PointToPointExchanger ex;
    ex.initialize(sg);
    JacobiSolver solver;

    const Index iters = solver.solve(sg, ex, 30000, 1e-7);
    EXPECT_LT(iters, Index(30000));
    EXPECT_GT(solver.lastResidual(), 0.0);
    EXPECT_LE(solver.lastResidual(), 1e-7);

    Real localSum = 0.0;
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gI = offsetX + static_cast<long long>(i) - 1;
            const long long gJ = offsetY + static_cast<long long>(j) - 1;
            const Real diff = sg.u().data()[sg.index(i, j)] - exact(gI, gJ);
            localSum += diff * diff;
        }
    }

    Real globalSum = 0.0;
    MPI_Allreduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, cart);
    const Real l2Error = std::sqrt(globalSum / (64.0 * 64.0));
    EXPECT_LT(l2Error, 1e-3);

    MPI_Comm_free(&cart);
}

TEST(SolverMpiTest, FixedIterationsMatchesSerialReference) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "FixedIterationsMatchesSerialReference requires exactly 4 processes";
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

    int coords[2] = {0, 0};
    MPI_Cart_coords(cart, rank, 2, coords);
    const long long offsetX = static_cast<long long>(coords[0]) * 64;
    const long long offsetY = static_cast<long long>(coords[1]) * 64;

    Subgrid sg(64, 64, 1, 1, cart);
    sg.setNeighbors(left, right, down, up);
    setupPoisson(sg);

    PointToPointExchanger ex;
    ex.initialize(sg);
    JacobiSolver solver;

    const Index iters = solver.solve(sg, ex, 200, 0.0);
    EXPECT_EQ(iters, Index(200));

    Subgrid ref(128, 128, 1, 1, MPI_COMM_SELF);
    ref.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupPoisson(ref);

    PointToPointExchanger refEx;
    refEx.initialize(ref);
    JacobiSolver refSolver;
    refSolver.solve(ref, refEx, 200, 0.0);

    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gI = offsetX + static_cast<long long>(i) - 1;
            const long long gJ = offsetY + static_cast<long long>(j) - 1;
            EXPECT_NEAR(sg.u().data()[sg.index(i, j)],
                        ref.u().data()[ref.index(static_cast<Index>(gI + 1),
                                                 static_cast<Index>(gJ + 1))],
                        1e-12);
        }
    }

    MPI_Comm_free(&cart);
}

TEST(SolverMpiTest, ThreeDNoZDecompositionProcNullAndRun) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "ThreeDNoZDecompositionProcNullAndRun requires exactly 4 processes";
    }

    int dims[3] = {2, 2, 1};
    int periods[3] = {0, 0, 0};
    MPI_Comm cart = MPI_COMM_NULL;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart);

    int left = MPI_PROC_NULL;
    int right = MPI_PROC_NULL;
    int down = MPI_PROC_NULL;
    int up = MPI_PROC_NULL;
    int back = MPI_PROC_NULL;
    int front = MPI_PROC_NULL;
    MPI_Cart_shift(cart, 0, 1, &left, &right);
    MPI_Cart_shift(cart, 1, 1, &down, &up);
    MPI_Cart_shift(cart, 2, 1, &back, &front);

    EXPECT_EQ(back, MPI_PROC_NULL);
    EXPECT_EQ(front, MPI_PROC_NULL);

    Subgrid sg(32, 32, 8, 1, cart);
    sg.setNeighbors(left, right, down, up, back, front);
    setupPoisson3D(sg);

    PointToPointExchanger ex;
    ex.initialize(sg);
    JacobiSolver solver;

    const Index iters = solver.solve(sg, ex, 100, 0.0);
    EXPECT_EQ(iters, Index(100));

    for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                EXPECT_TRUE(std::isfinite(sg.u().data()[sg.index(i, j, k)]));
            }
        }
    }

    MPI_Comm_free(&cart);
}

TEST(SolverMpiTest, TinyDecompositionSmoke) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "TinyDecompositionSmoke requires exactly 4 processes";
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

    // Global 4x4 with 2x2 ranks: 2x2 per rank (E5 scenario).
    {
        Subgrid sg(2, 2, 1, 1, cart);
        sg.setNeighbors(left, right, down, up);
        setupPoisson(sg);

        PointToPointExchanger ex;
        ex.initialize(sg);
        JacobiSolver solver;

        const Index iters = solver.solve(sg, ex, 50, 0.0);
        EXPECT_EQ(iters, Index(50));
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                EXPECT_TRUE(std::isfinite(sg.u().data()[sg.index(i, j)]));
            }
        }
    }

    // Global 2x2 with 2x2 ranks: 1x1 per rank (degenerate inner box).
    {
        Subgrid sg(1, 1, 1, 1, cart);
        sg.setNeighbors(left, right, down, up);
        setupPoisson(sg);

        PointToPointExchanger ex;
        ex.initialize(sg);
        JacobiSolver solver;

        const Index iters = solver.solve(sg, ex, 10, 0.0);
        EXPECT_EQ(iters, Index(10));
        EXPECT_TRUE(std::isfinite(sg.u().data()[sg.index(sg.iBegin(), sg.jBegin())]));
    }

    MPI_Comm_free(&cart);
}
