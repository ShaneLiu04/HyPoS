#include <gtest/gtest.h>
#include "grid/grid.hpp"
#include "grid/subgrid.hpp"
#include "grid/partition.hpp"
#include "utils/mpi_env.hpp"
#include <mpi.h>
#include <algorithm>
#include <utility>
#include <vector>

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

    // AR004 T006: partition() now transfers cartComm ownership to the caller
    if (info.cartComm != MPI_COMM_NULL) {
        MPI_Comm_free(&info.cartComm);
    }
}

// B3 (AR004 T006, design §4.2.2-6, D5/D6): np=4 topology consistency —
// UniformPartition::partition() must expose the Cartesian communicator it
// created (SubgridInfo::cartComm, ownership transferred to the caller), and
// the topology data derived from it must be mutually consistent:
//   * per-dimension sum of local sizes equals the global size;
//   * offsetX/offsetY form a contiguous chain within every grid row/column
//     (ranks grouped by the cartComm coords of the orthogonal dimension);
//   * Cart_shift neighbors reciprocate (my right neighbor's left is me,
//     checked symmetrically in both topology dimensions);
//   * coords-offset correspondence is implied by the chain check, which
//     orders ranks by the coords obtained from info.cartComm.
TEST(PartitionTopologyTest, Np4TopologyIsConsistent) {
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "Np4TopologyIsConsistent requires exactly 4 processes";
    }

    Grid grid;
    grid.nx = 1024;
    grid.ny = 1024;
    grid.nz = 1;

    UniformPartition partition;
    SubgridInfo info = partition.partition(grid, MPI_COMM_WORLD, rank);

    // D5: partition() must hand the Cartesian communicator to the caller.
    ASSERT_NE(info.cartComm, MPI_COMM_NULL);

    // --- Cart rank/coords must be queried on cartComm (D6), never on the
    // raw communicator: with reorder=1 the two rank spaces may differ.
    int cartRank = -1;
    MPI_Comm_rank(info.cartComm, &cartRank);
    int coords[3] = {0, 0, 0};
    MPI_Cart_coords(info.cartComm, cartRank, 3, coords);

    // --- Gather (coords, offsets, local sizes) from all ranks.
    std::vector<int> allC0(size), allC1(size);
    std::vector<Index> allOffX(size), allOffY(size), allLnX(size), allLnY(size);
    MPI_Allgather(&coords[0], 1, MPI_INT, allC0.data(), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&coords[1], 1, MPI_INT, allC1.data(), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&info.offsetX, 1, MPI_UNSIGNED_LONG, allOffX.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&info.offsetY, 1, MPI_UNSIGNED_LONG, allOffY.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&info.nxLocal, 1, MPI_UNSIGNED_LONG, allLnX.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&info.nyLocal, 1, MPI_UNSIGNED_LONG, allLnY.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    // --- Grouped size sums. A plain global sum is WRONG for 2D tiling:
    // every rank contributes its full x-extent, so Σ nxLocal over all ranks
    // is grid.nx * dims[1]. The correct invariant is per grid row (fixed
    // coords[1]): Σ nxLocal == grid.nx; per column (fixed coords[0]):
    // Σ nyLocal == grid.ny; per (coords[0], coords[1]) column pair:
    // Σ nzLocal == grid.nz.
    {
        std::vector<int> uniqueC1;
        for (int c : allC1) {
            if (std::find(uniqueC1.begin(), uniqueC1.end(), c) == uniqueC1.end()) {
                uniqueC1.push_back(c);
            }
        }
        for (int c1 : uniqueC1) {
            Index rowSum = 0;
            for (int r = 0; r < size; ++r) {
                if (allC1[r] == c1) rowSum += allLnX[r];
            }
            EXPECT_EQ(rowSum, grid.nx)
                << "row coords[1]=" << c1 << " local x-extents must sum to grid.nx";
        }
        std::vector<int> uniqueC0;
        for (int c : allC0) {
            if (std::find(uniqueC0.begin(), uniqueC0.end(), c) == uniqueC0.end()) {
                uniqueC0.push_back(c);
            }
        }
        for (int c0 : uniqueC0) {
            Index colSumY = 0;
            for (int r = 0; r < size; ++r) {
                if (allC0[r] == c0) {
                    colSumY += allLnY[r];
                }
            }
            EXPECT_EQ(colSumY, grid.ny)
                << "column coords[0]=" << c0 << " local y-extents must sum to grid.ny";
        }
        // 2D grid: the z dimension is not decomposed (dims[2] == 1), so
        // every rank owns the full z extent.
        EXPECT_EQ(info.nzLocal, grid.nz);
    }

    std::vector<int> order(size);
    for (int i = 0; i < size; ++i) order[i] = i;

    // --- Offset chain along X within every grid row (fixed coords[1]):
    // ranks of a row sorted by coords[0] must tile [0, grid.nx) exactly.
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return std::make_pair(allC1[a], allC0[a]) < std::make_pair(allC1[b], allC0[b]);
    });
    for (std::size_t i = 0; i < order.size();) {
        std::size_t j = i;
        while (j < order.size() && allC1[order[j]] == allC1[order[i]]) ++j;
        EXPECT_EQ(allOffX[order[i]], Index(0))
            << "row coords[1]=" << allC1[order[i]] << " must start at offsetX 0";
        for (std::size_t k = i + 1; k < j; ++k) {
            EXPECT_EQ(allOffX[order[k]], allOffX[order[k - 1]] + allLnX[order[k - 1]])
                << "x-offset chain broken between coords[0]=" << allC0[order[k - 1]]
                << " and coords[0]=" << allC0[order[k]] << " in row coords[1]=" << allC1[order[i]];
        }
        EXPECT_EQ(allOffX[order[j - 1]] + allLnX[order[j - 1]], grid.nx)
            << "row coords[1]=" << allC1[order[i]] << " must end exactly at grid.nx";
        i = j;
    }

    // --- Analogous offset chain along Y within every grid column
    // (fixed coords[0]), sorted by coords[1], tiling [0, grid.ny).
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return std::make_pair(allC0[a], allC1[a]) < std::make_pair(allC0[b], allC1[b]);
    });
    for (std::size_t i = 0; i < order.size();) {
        std::size_t j = i;
        while (j < order.size() && allC0[order[j]] == allC0[order[i]]) ++j;
        EXPECT_EQ(allOffY[order[i]], Index(0))
            << "column coords[0]=" << allC0[order[i]] << " must start at offsetY 0";
        for (std::size_t k = i + 1; k < j; ++k) {
            EXPECT_EQ(allOffY[order[k]], allOffY[order[k - 1]] + allLnY[order[k - 1]])
                << "y-offset chain broken between coords[1]=" << allC1[order[k - 1]]
                << " and coords[1]=" << allC1[order[k]] << " in column coords[0]=" << allC0[order[i]];
        }
        EXPECT_EQ(allOffY[order[j - 1]] + allLnY[order[j - 1]], grid.ny)
            << "column coords[0]=" << allC0[order[i]] << " must end exactly at grid.ny";
        i = j;
    }

    // --- Cart_shift neighbor reciprocity on info.cartComm: each rank sends
    // its own cartComm rank to one neighbor and receives the rank that
    // neighbor reports for itself; the received value must equal the
    // Cart_shift source rank (i.e. my right neighbor's left is me, and
    // symmetrically my left neighbor's right is me). Edge-safe: skip the
    // assertion when the shift yields MPI_PROC_NULL.
    int left = MPI_PROC_NULL, right = MPI_PROC_NULL;
    int down = MPI_PROC_NULL, up = MPI_PROC_NULL;
    MPI_Cart_shift(info.cartComm, 0, 1, &left, &right);
    MPI_Cart_shift(info.cartComm, 1, 1, &down, &up);

    int token = cartRank;
    int fromLeft = -1, fromRight = -1, fromDown = -1, fromUp = -1;
    MPI_Sendrecv(&token, 1, MPI_INT, right, 0, &fromLeft, 1, MPI_INT, left, 0,
                 info.cartComm, MPI_STATUS_IGNORE);
    if (left != MPI_PROC_NULL) {
        EXPECT_EQ(fromLeft, left) << "rank heard from left is not my Cart_shift left neighbor";
    }
    MPI_Sendrecv(&token, 1, MPI_INT, left, 1, &fromRight, 1, MPI_INT, right, 1,
                 info.cartComm, MPI_STATUS_IGNORE);
    if (right != MPI_PROC_NULL) {
        EXPECT_EQ(fromRight, right) << "rank heard from right is not my Cart_shift right neighbor";
    }
    MPI_Sendrecv(&token, 1, MPI_INT, up, 2, &fromDown, 1, MPI_INT, down, 2,
                 info.cartComm, MPI_STATUS_IGNORE);
    if (down != MPI_PROC_NULL) {
        EXPECT_EQ(fromDown, down) << "rank heard from down is not my Cart_shift down neighbor";
    }
    MPI_Sendrecv(&token, 1, MPI_INT, down, 3, &fromUp, 1, MPI_INT, up, 3,
                 info.cartComm, MPI_STATUS_IGNORE);
    if (up != MPI_PROC_NULL) {
        EXPECT_EQ(fromUp, up) << "rank heard from up is not my Cart_shift up neighbor";
    }

    // D5 ownership contract: the caller frees the Cartesian communicator
    // that partition() transferred to it.
    MPI_Comm_free(&info.cartComm);
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
