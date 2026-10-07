#include <gtest/gtest.h>
#include "grid/grid.hpp"
#include "grid/subgrid.hpp"
#include "grid/partition.hpp"
#include "utils/mpi_env.hpp"
#include <mpi.h>

using namespace hypo;

TEST(GridTest, Dimensions) {
    Grid g;
    g.nx = 100;
    g.ny = 200;
    g.nz = 1;
    EXPECT_EQ(g.totalCells(), 20000);
    EXPECT_TRUE(g.is2D());
    EXPECT_FALSE(g.is3D());
}

TEST(GridTest, PartitionUniform) {
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Grid grid;
    grid.nx = 1024;
    grid.ny = 1024;
    grid.nz = 1;

    UniformPartition partition;
    SubgridInfo info = partition.partition(grid, MPI_COMM_WORLD, rank);

    // Verify that local sizes sum to global
    Index totalNx = 0, totalNy = 0;
    MPI_Reduce(&info.nxLocal, &totalNx, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&info.nyLocal, &totalNy, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        EXPECT_EQ(totalNx, grid.nx);
        EXPECT_EQ(totalNy, grid.ny);
    }
}

TEST(SubgridTest, MemoryAllocation) {
    Subgrid sg(100, 100, 1, 1, MPI_COMM_SELF);
    EXPECT_EQ(sg.nxLocal(), 100);
    EXPECT_EQ(sg.nyLocal(), 100);
    EXPECT_EQ(sg.nzLocal(), 1);
    EXPECT_EQ(sg.nxTotal(), 102);
    EXPECT_EQ(sg.nyTotal(), 102);
    EXPECT_EQ(sg.nzTotal(), 3); // 1 + 2*halo
    EXPECT_EQ(sg.totalCells(), 102 * 102 * 3);
}
