#include <gtest/gtest.h>
#include "solver/solver.hpp"
#include "grid/subgrid.hpp"
#include "comm/halo_exchanger.hpp"
#include <mpi.h>
#include <cmath>
#include <limits>
#include <string>

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
