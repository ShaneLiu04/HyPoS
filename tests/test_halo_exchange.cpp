#include <gtest/gtest.h>
#include "comm/halo_exchanger.hpp"
#include "grid/subgrid.hpp"
#include <mpi.h>
#include <vector>

using namespace hypo;

TEST(CommTest, HaloExchangeCorrectness) {
    int argc = 0;
    char** argv = nullptr;
    MPI_Init(&argc, &argv);

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 4) {
        GTEST_SKIP() << "This test requires exactly 4 MPI processes";
    }

    // Create 2x2 Cartesian grid
    int dims[2] = {2, 2};
    int periods[2] = {0, 0};
    MPI_Comm cartComm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cartComm);

    int left, right, down, up;
    MPI_Cart_shift(cartComm, 0, 1, &left, &right);
    MPI_Cart_shift(cartComm, 1, 1, &down, &up);

    // Each process gets 4x4 local grid
    Subgrid subgrid(4, 4, 1, 1, cartComm);
    subgrid.setNeighbors(left, right, down, up);

    // Initialize interior with rank number
    Real* u = subgrid.u().data();
    for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
        for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
            u[subgrid.index(i, j)] = static_cast<Real>(rank);
        }
    }

    // Exchange halos
    PointToPointExchanger exchanger;
    exchanger.initialize(subgrid);
    exchanger.exchange(subgrid);

    // Verify halo values match neighbor ranks
    // For a 2x2 grid:
    // rank 0: neighbors right=1, up=2
    // rank 1: neighbors left=0, up=3
    // rank 2: neighbors right=3, down=0
    // rank 3: neighbors left=2, down=1

    if (right != MPI_PROC_NULL) {
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            EXPECT_DOUBLE_EQ(u[subgrid.index(subgrid.iEnd(), j)], static_cast<Real>(right));
        }
    }
    if (left != MPI_PROC_NULL) {
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            EXPECT_DOUBLE_EQ(u[subgrid.index(0, j)], static_cast<Real>(left));
        }
    }
    if (up != MPI_PROC_NULL) {
        for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
            EXPECT_DOUBLE_EQ(u[subgrid.index(i, subgrid.jEnd())], static_cast<Real>(up));
        }
    }
    if (down != MPI_PROC_NULL) {
        for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
            EXPECT_DOUBLE_EQ(u[subgrid.index(i, 0)], static_cast<Real>(down));
        }
    }

    MPI_Comm_free(&cartComm);
    MPI_Finalize();
}
