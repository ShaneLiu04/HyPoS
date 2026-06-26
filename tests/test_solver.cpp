#include <gtest/gtest.h>
#include "solver/solver.hpp"
#include "grid/subgrid.hpp"
#include "comm/halo_exchanger.hpp"
#include <mpi.h>
#include <cmath>

using namespace hypo;

TEST(SolverTest, JacobiConvergence) {
    int argc = 0;
    char** argv = nullptr;
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create a small 2D subgrid (no neighbors, so use MPI_COMM_SELF)
    Subgrid subgrid(32, 32, 1, 1, MPI_COMM_SELF);
    subgrid.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);

    // Set up a simple problem: u = 0 on boundary, f = 1 everywhere
    subgrid.zeroInitialize();
    subgrid.applyDirichletBC(0.0);

    Real* rhs = subgrid.rhs().data();
    for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
        for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
            rhs[subgrid.index(i, j)] = 1.0;
        }
    }

    PointToPointExchanger exchanger;
    exchanger.initialize(subgrid);

    JacobiSolver solver;
    Index iters = solver.solve(subgrid, exchanger, 5000, 1e-5);

    EXPECT_GT(iters, 0);
    EXPECT_LT(iters, 5000); // Should converge before max iterations

    // Check that interior values are positive (since f=1 and boundary=0)
    const Real* u = subgrid.u().data();
    for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
        for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
            EXPECT_GT(u[subgrid.index(i, j)], 0.0);
        }
    }

    MPI_Finalize();
}
