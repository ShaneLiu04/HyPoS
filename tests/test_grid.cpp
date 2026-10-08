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

// U7 (AR004 T005, design §6): zeroInitialize must zero the ENTIRE buffer of
// all three fields — including halo cells and, for 2D subgrids, the z-halo
// planes (a 2D subgrid still allocates nzTotal = 3 planes with halo 1).
// Guards the OpenMP parallelization: a wrong loop range (e.g. covering only
// the k=0 plane) would leave parts of the buffer dirty. Regression-capture
// semantics: passes on the baseline fill()-based code and pins the contract
// through the parallel rewrite.
TEST(SubgridTest, ZeroInitializeCoversFullBuffer2DAnd3D) {
    const auto poisonCheck = [](Subgrid& sg, const char* label) {
        const Index total = sg.totalCells();
        ASSERT_GT(total, Index(0));
        Real* u = sg.u().data();
        Real* un = sg.uNext().data();
        Real* r = sg.rhs().data();
        for (Index idx = 0; idx < total; ++idx) {
            u[idx] = 1.5 + static_cast<Real>(idx % 7);
            un[idx] = -2.25;
            r[idx] = 7.75;
        }
        sg.zeroInitialize();
        for (Index idx = 0; idx < total; ++idx) {
            EXPECT_EQ(u[idx], 0.0) << label << " u idx=" << idx;
            EXPECT_EQ(un[idx], 0.0) << label << " uNext idx=" << idx;
            EXPECT_EQ(r[idx], 0.0) << label << " rhs idx=" << idx;
        }
    };

    // 2D with halo 2: nzTotal = 5 planes, halo width exercises range edges.
    Subgrid sg2(16, 16, 1, 2, MPI_COMM_SELF);
    poisonCheck(sg2, "2D halo=2");

    // 3D with halo 1.
    Subgrid sg3(8, 8, 8, 1, MPI_COMM_SELF);
    poisonCheck(sg3, "3D halo=1");
}
