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

    int coords[2] = {0, 0};
    MPI_Cart_coords(cart, rank, 2, coords);

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

        // Reference: redundant full-domain serial solve, compared per cell.
        Subgrid ref(4, 4, 1, 1, MPI_COMM_SELF);
        ref.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
        setupPoisson(ref);

        PointToPointExchanger refEx;
        refEx.initialize(ref);
        JacobiSolver refSolver;
        refSolver.solve(ref, refEx, 50, 0.0);

        const Index offsetX = static_cast<Index>(coords[0]) * 2;
        const Index offsetY = static_cast<Index>(coords[1]) * 2;
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                const Index gI = offsetX + i - 1;
                const Index gJ = offsetY + j - 1;
                EXPECT_NEAR(sg.u().data()[sg.index(i, j)],
                            ref.u().data()[ref.index(gI + 1, gJ + 1)], 1e-12);
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

        Subgrid ref(2, 2, 1, 1, MPI_COMM_SELF);
        ref.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
        setupPoisson(ref);

        PointToPointExchanger refEx;
        refEx.initialize(ref);
        JacobiSolver refSolver;
        refSolver.solve(ref, refEx, 10, 0.0);

        const Index offsetX = static_cast<Index>(coords[0]);
        const Index offsetY = static_cast<Index>(coords[1]);
        EXPECT_NEAR(sg.u().data()[sg.index(sg.iBegin(), sg.jBegin())],
                    ref.u().data()[ref.index(offsetX + 1, offsetY + 1)], 1e-12);
    }

    MPI_Comm_free(&cart);
}

namespace {

Real sineExactMpi(long long gI, long long gJ) {
    const double x = static_cast<double>(gI + 1) / 65.0;
    const double y = static_cast<double>(gJ + 1) / 65.0;
    return std::sin(M_PI * x) * std::sin(M_PI * y);
}

void setupManufacturedSineMpi(Subgrid& sg, Index offsetX, Index offsetY) {
    sg.setOffsets(offsetX, offsetY, 0);
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);
    Real* rhs = sg.rhs().data();
    const long long hw = static_cast<long long>(sg.haloWidth());
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gI = static_cast<long long>(offsetX) +
                                 static_cast<long long>(i) - hw;
            const long long gJ = static_cast<long long>(offsetY) +
                                 static_cast<long long>(j) - hw;
            rhs[sg.index(i, j)] = sineExactMpi(gI - 1, gJ) + sineExactMpi(gI + 1, gJ) +
                                  sineExactMpi(gI, gJ - 1) + sineExactMpi(gI, gJ + 1) -
                                  4.0 * sineExactMpi(gI, gJ);
        }
    }
}

Real globalL2VsSineMpi(const Subgrid& sg, MPI_Comm comm) {
    Real local = 0.0;
    const long long hw = static_cast<long long>(sg.haloWidth());
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gI = static_cast<long long>(sg.offsetX()) +
                                 static_cast<long long>(i) - hw;
            const long long gJ = static_cast<long long>(sg.offsetY()) +
                                 static_cast<long long>(j) - hw;
            const Real diff = sg.u().data()[sg.index(i, j)] - sineExactMpi(gI, gJ);
            local += diff * diff;
        }
    }
    Real global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, comm);
    return std::sqrt(global / (64.0 * 64.0));
}

} // namespace

TEST(SolverMpiTest, AltSolversConvergeFasterThanJacobi2D) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "AltSolversConvergeFasterThanJacobi2D requires exactly 4 processes";
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
    const Index offsetX = static_cast<Index>(coords[0]) * 32;
    const Index offsetY = static_cast<Index>(coords[1]) * 32;

    const int maxIter = 30000;
    const Real tol = 1e-8;

    Subgrid sgJ(32, 32, 1, 1, cart);
    sgJ.setNeighbors(left, right, down, up);
    setupManufacturedSineMpi(sgJ, offsetX, offsetY);
    PointToPointExchanger exJ;
    exJ.initialize(sgJ);
    JacobiSolver jacobi;
    const Index itersJ = jacobi.solve(sgJ, exJ, maxIter, tol);
    ASSERT_LT(itersJ, Index(maxIter));

    Subgrid sgR(32, 32, 1, 1, cart);
    sgR.setNeighbors(left, right, down, up);
    setupManufacturedSineMpi(sgR, offsetX, offsetY);
    PointToPointExchanger exR;
    exR.initialize(sgR);
    RedBlackGSSolver rbgs;
    const Index itersR = rbgs.solve(sgR, exR, maxIter, tol);
    EXPECT_LT(itersR, Index(maxIter));
    EXPECT_LE(itersR, itersJ * 60 / 100);
    EXPECT_LT(globalL2VsSineMpi(sgR, cart), 1e-3);

    Subgrid sgC(32, 32, 1, 1, cart);
    sgC.setNeighbors(left, right, down, up);
    setupManufacturedSineMpi(sgC, offsetX, offsetY);
    PointToPointExchanger exC;
    exC.initialize(sgC);
    CGSolver cg;
    const Index itersC = cg.solve(sgC, exC, maxIter, tol);
    EXPECT_LT(itersC, Index(maxIter));
    EXPECT_LE(itersC, itersJ / 10);
    EXPECT_LT(globalL2VsSineMpi(sgC, cart), 1e-3);

    Subgrid sgT(32, 32, 1, 1, cart);
    sgT.setNeighbors(left, right, down, up);
    setupManufacturedSineMpi(sgT, offsetX, offsetY);
    PointToPointExchanger exT;
    exT.initialize(sgT);
    JacobiSolver tight;
    tight.solve(sgT, exT, 100000, 1e-12);

    Subgrid sgR2(32, 32, 1, 1, cart);
    sgR2.setNeighbors(left, right, down, up);
    setupManufacturedSineMpi(sgR2, offsetX, offsetY);
    PointToPointExchanger exR2;
    exR2.initialize(sgR2);
    RedBlackGSSolver rbgsTight;
    rbgsTight.solve(sgR2, exR2, 30000, 1e-10);

    for (Index j = sgR.jBegin(); j < sgR.jEnd(); ++j) {
        for (Index i = sgR.iBegin(); i < sgR.iEnd(); ++i) {
            EXPECT_NEAR(sgR2.u().data()[sgR2.index(i, j)],
                        sgT.u().data()[sgT.index(i, j)], 1e-8);
            EXPECT_NEAR(sgC.u().data()[sgC.index(i, j)],
                        sgT.u().data()[sgT.index(i, j)], 1e-8);
        }
    }

    MPI_Comm_free(&cart);
}

TEST(SolverMpiTest, AltSolvers3DSmoke) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "AltSolvers3DSmoke requires exactly 4 processes";
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

    {
        Subgrid sg(32, 32, 8, 1, cart);
        sg.setNeighbors(left, right, down, up, back, front);
        setupPoisson3D(sg);
        PointToPointExchanger ex;
        ex.initialize(sg);
        RedBlackGSSolver rbgs;
        const Index iters = rbgs.solve(sg, ex, 100, 0.0);
        EXPECT_EQ(iters, Index(100));
        for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
            for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
                for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                    EXPECT_TRUE(std::isfinite(sg.u().data()[sg.index(i, j, k)]));
                }
            }
        }
    }

    {
        Subgrid sg(32, 32, 8, 1, cart);
        sg.setNeighbors(left, right, down, up, back, front);
        setupPoisson3D(sg);
        PointToPointExchanger ex;
        ex.initialize(sg);
        CGSolver cg;
        const Index iters = cg.solve(sg, ex, 100, 0.0);
        EXPECT_EQ(iters, Index(100));
        for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
            for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
                for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                    EXPECT_TRUE(std::isfinite(sg.u().data()[sg.index(i, j, k)]));
                }
            }
        }
    }

    MPI_Comm_free(&cart);
}
