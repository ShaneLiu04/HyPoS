#include <gtest/gtest.h>
#include <mpi.h>

#include "comm/halo_exchanger.hpp"
#include "grid/subgrid.hpp"

#include <vector>

using namespace hypo;

namespace {

Real expectedValue(long long gI, long long gJ) {
    return 1000.0 * static_cast<Real>(gI) + static_cast<Real>(gJ);
}

void fillPadded2D(Subgrid& sg) {
    for (Index j = 0; j < sg.nyTotal(); ++j) {
        for (Index i = 0; i < sg.nxTotal(); ++i) {
            sg.at(i, j) = expectedValue(static_cast<long long>(i), static_cast<long long>(j));
        }
    }
}

} // namespace

TEST(HaloExchangeTest, DirectionSemanticsSelfLoop2D) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Subgrid sg(4, 4, 1, 1, MPI_COMM_WORLD);
    sg.setNeighbors(rank, rank, rank, rank);
    fillPadded2D(sg);

    PointToPointExchanger ex;
    ex.initialize(sg);
    ex.exchange(sg);

    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        EXPECT_DOUBLE_EQ(sg.at(sg.iEnd(), j),
                         expectedValue(static_cast<long long>(sg.iBegin()), static_cast<long long>(j)))
            << "right halo at (iEnd, " << j << ")";
    }
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        EXPECT_DOUBLE_EQ(sg.at(0, j),
                         expectedValue(static_cast<long long>(sg.iEnd() - 1), static_cast<long long>(j)))
            << "left halo at (0, " << j << ")";
    }
    for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
        EXPECT_DOUBLE_EQ(sg.at(i, sg.jEnd()),
                         expectedValue(static_cast<long long>(i), static_cast<long long>(sg.jBegin())))
            << "up halo at (" << i << ", jEnd)";
    }
    for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
        EXPECT_DOUBLE_EQ(sg.at(i, 0),
                         expectedValue(static_cast<long long>(i), static_cast<long long>(sg.jEnd() - 1)))
            << "down halo at (" << i << ", 0)";
    }
}

TEST(HaloExchangeTest, FaceBuffersSafeSelfLoop2DMultiHalo) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Subgrid sg(6, 6, 1, 2, MPI_COMM_WORLD);
    sg.setNeighbors(rank, rank, rank, rank);
    fillPadded2D(sg);

    PointToPointExchanger ex;
    ex.initialize(sg);
    ex.exchange(sg);

    const Index hw = static_cast<Index>(sg.haloWidth());

    for (Index h = 0; h < hw; ++h) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            EXPECT_DOUBLE_EQ(sg.at(sg.iEnd() + h, j),
                             expectedValue(static_cast<long long>(sg.iBegin() + h), static_cast<long long>(j)))
                << "right halo h=" << h << " j=" << j;
        }
    }
    for (Index h = 0; h < hw; ++h) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            EXPECT_DOUBLE_EQ(sg.at(h, j),
                             expectedValue(static_cast<long long>(sg.iEnd() - hw + h), static_cast<long long>(j)))
                << "left halo h=" << h << " j=" << j;
        }
    }
    for (Index h = 0; h < hw; ++h) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            EXPECT_DOUBLE_EQ(sg.at(i, sg.jEnd() + h),
                             expectedValue(static_cast<long long>(i), static_cast<long long>(sg.jBegin() + h)))
                << "up halo i=" << i << " h=" << h;
        }
    }
    for (Index h = 0; h < hw; ++h) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            EXPECT_DOUBLE_EQ(sg.at(i, h),
                             expectedValue(static_cast<long long>(i), static_cast<long long>(sg.jEnd() - hw + h)))
                << "down halo i=" << i << " h=" << h;
        }
    }
}

TEST(HaloExchangeTest, ProcNullKeepsHaloUnchanged) {
    Subgrid sg(4, 4, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    fillPadded2D(sg);

    std::vector<Real> snapshot(sg.nxTotal() * sg.nyTotal());
    for (Index j = 0; j < sg.nyTotal(); ++j) {
        for (Index i = 0; i < sg.nxTotal(); ++i) {
            snapshot[sg.index(i, j)] = sg.at(i, j);
        }
    }

    PointToPointExchanger ex;
    ex.initialize(sg);
    ex.exchange(sg);

    for (Index j = 0; j < sg.nyTotal(); ++j) {
        for (Index i = 0; i < sg.nxTotal(); ++i) {
            const Index idx = sg.index(i, j);
            EXPECT_DOUBLE_EQ(sg.u().data()[idx], snapshot[idx]) << "cell changed at i=" << i << " j=" << j;
        }
    }
}

TEST(HaloExchangeTest, CollectiveExchangerDelegatesToP2P) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Subgrid sg(4, 4, 1, 1, MPI_COMM_WORLD);
    sg.setNeighbors(rank, rank, rank, rank);
    fillPadded2D(sg);

    CollectiveExchanger ex;
    ex.initialize(sg);
    ex.exchange(sg);

    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        EXPECT_DOUBLE_EQ(sg.at(sg.iEnd(), j),
                         expectedValue(static_cast<long long>(sg.iBegin()), static_cast<long long>(j)))
            << "right halo at (iEnd, " << j << ")";
        EXPECT_DOUBLE_EQ(sg.at(0, j),
                         expectedValue(static_cast<long long>(sg.iEnd() - 1), static_cast<long long>(j)))
            << "left halo at (0, " << j << ")";
    }
}

TEST(HaloExchangeTest, MpiNonUniform4RanksHalos) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "MpiNonUniform4RanksHalos requires exactly 4 processes";
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
    const int offsetX = coords[0] * 4;
    const int offsetY = coords[1] * 4;

    Subgrid sg(4, 4, 1, 1, cart);
    sg.setNeighbors(left, right, down, up);

    const Real sentinel = -777777.0;
    for (Index j = 0; j < sg.nyTotal(); ++j) {
        for (Index i = 0; i < sg.nxTotal(); ++i) {
            const bool interior = i >= sg.iBegin() && i < sg.iEnd() &&
                                  j >= sg.jBegin() && j < sg.jEnd();
            if (interior) {
                const long long gI = static_cast<long long>(offsetX) + static_cast<long long>(i) - 1;
                const long long gJ = static_cast<long long>(offsetY) + static_cast<long long>(j) - 1;
                sg.at(i, j) = expectedValue(gI, gJ);
            } else {
                sg.at(i, j) = sentinel;
            }
        }
    }

    PointToPointExchanger ex;
    ex.initialize(sg);
    ex.exchange(sg);

    // Exchanged halo cells must hold the neighbor's interior values; halo
    // cells owned by a physical boundary (PROC_NULL) must keep the sentinel.
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        const long long gJ = static_cast<long long>(offsetY) + static_cast<long long>(j) - 1;
        if (left != MPI_PROC_NULL) {
            EXPECT_DOUBLE_EQ(sg.at(0, j), expectedValue(offsetX - 1, gJ))
                << "rank " << rank << " left halo j=" << j;
        } else {
            EXPECT_DOUBLE_EQ(sg.at(0, j), sentinel) << "rank " << rank << " left halo j=" << j;
        }
        if (right != MPI_PROC_NULL) {
            EXPECT_DOUBLE_EQ(sg.at(sg.iEnd(), j), expectedValue(offsetX + 4, gJ))
                << "rank " << rank << " right halo j=" << j;
        } else {
            EXPECT_DOUBLE_EQ(sg.at(sg.iEnd(), j), sentinel) << "rank " << rank << " right halo j=" << j;
        }
    }
    for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
        const long long gI = static_cast<long long>(offsetX) + static_cast<long long>(i) - 1;
        if (down != MPI_PROC_NULL) {
            EXPECT_DOUBLE_EQ(sg.at(i, 0), expectedValue(gI, offsetY - 1))
                << "rank " << rank << " down halo i=" << i;
        } else {
            EXPECT_DOUBLE_EQ(sg.at(i, 0), sentinel) << "rank " << rank << " down halo i=" << i;
        }
        if (up != MPI_PROC_NULL) {
            EXPECT_DOUBLE_EQ(sg.at(i, sg.jEnd()), expectedValue(gI, offsetY + 4))
                << "rank " << rank << " up halo i=" << i;
        } else {
            EXPECT_DOUBLE_EQ(sg.at(i, sg.jEnd()), sentinel) << "rank " << rank << " up halo i=" << i;
        }
    }

    MPI_Comm_free(&cart);
}

namespace {

Real expectedValue3D(long long gI, long long gJ, long long gK) {
    return 10000.0 * static_cast<Real>(gI) + 100.0 * static_cast<Real>(gJ) + static_cast<Real>(gK);
}

void fillPadded3D(Subgrid& sg) {
    for (Index k = 0; k < sg.nzTotal(); ++k) {
        for (Index j = 0; j < sg.nyTotal(); ++j) {
            for (Index i = 0; i < sg.nxTotal(); ++i) {
                sg.at(i, j, k) = expectedValue3D(static_cast<long long>(i),
                                                 static_cast<long long>(j),
                                                 static_cast<long long>(k));
            }
        }
    }
}

} // namespace

TEST(HaloExchangeTest, SelfLoop3DDirectionSemantics) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Subgrid sg(2, 2, 2, 1, MPI_COMM_WORLD);
    sg.setNeighbors(rank, rank, rank, rank, rank, rank);
    fillPadded3D(sg);

    PointToPointExchanger ex;
    ex.initialize(sg);
    ex.exchange(sg);

    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
            EXPECT_DOUBLE_EQ(sg.at(sg.iEnd(), j, k),
                             expectedValue3D(static_cast<long long>(sg.iBegin()),
                                             static_cast<long long>(j),
                                             static_cast<long long>(k)))
                << "right halo at (iEnd, " << j << ", " << k << ")";
            EXPECT_DOUBLE_EQ(sg.at(0, j, k),
                             expectedValue3D(static_cast<long long>(sg.iEnd() - 1),
                                             static_cast<long long>(j),
                                             static_cast<long long>(k)))
                << "left halo at (0, " << j << ", " << k << ")";
        }
    }
    for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
        for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
            EXPECT_DOUBLE_EQ(sg.at(i, sg.jEnd(), k),
                             expectedValue3D(static_cast<long long>(i),
                                             static_cast<long long>(sg.jBegin()),
                                             static_cast<long long>(k)))
                << "up halo at (" << i << ", jEnd, " << k << ")";
            EXPECT_DOUBLE_EQ(sg.at(i, 0, k),
                             expectedValue3D(static_cast<long long>(i),
                                             static_cast<long long>(sg.jEnd() - 1),
                                             static_cast<long long>(k)))
                << "down halo at (" << i << ", 0, " << k << ")";
        }
    }
    for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            EXPECT_DOUBLE_EQ(sg.at(i, j, sg.kEnd()),
                             expectedValue3D(static_cast<long long>(i),
                                             static_cast<long long>(j),
                                             static_cast<long long>(sg.kBegin())))
                << "front halo at (" << i << ", " << j << ", kEnd)";
            EXPECT_DOUBLE_EQ(sg.at(i, j, 0),
                             expectedValue3D(static_cast<long long>(i),
                                             static_cast<long long>(j),
                                             static_cast<long long>(sg.kEnd() - 1)))
                << "back halo at (" << i << ", " << j << ", 0)";
        }
    }
}

TEST(HaloExchangeTest, SelfLoop3DMultiHaloSmoke) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Subgrid sg(4, 4, 4, 2, MPI_COMM_WORLD);
    sg.setNeighbors(rank, rank, rank, rank, rank, rank);
    fillPadded3D(sg);

    PointToPointExchanger ex;
    ex.initialize(sg);
    ex.exchange(sg);

    const Index hw = static_cast<Index>(sg.haloWidth());

    for (Index h = 0; h < hw; ++h) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
                EXPECT_DOUBLE_EQ(sg.at(sg.iEnd() + h, j, k),
                                 expectedValue3D(static_cast<long long>(sg.iBegin() + h),
                                                 static_cast<long long>(j),
                                                 static_cast<long long>(k)))
                    << "right halo h=" << h << " j=" << j << " k=" << k;
                EXPECT_DOUBLE_EQ(sg.at(h, j, k),
                                 expectedValue3D(static_cast<long long>(sg.iEnd() - hw + h),
                                                 static_cast<long long>(j),
                                                 static_cast<long long>(k)))
                    << "left halo h=" << h << " j=" << j << " k=" << k;
            }
        }
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
                EXPECT_DOUBLE_EQ(sg.at(i, sg.jEnd() + h, k),
                                 expectedValue3D(static_cast<long long>(i),
                                                 static_cast<long long>(sg.jBegin() + h),
                                                 static_cast<long long>(k)))
                    << "up halo h=" << h << " i=" << i << " k=" << k;
                EXPECT_DOUBLE_EQ(sg.at(i, h, k),
                                 expectedValue3D(static_cast<long long>(i),
                                                 static_cast<long long>(sg.jEnd() - hw + h),
                                                 static_cast<long long>(k)))
                    << "down halo h=" << h << " i=" << i << " k=" << k;
            }
        }
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
                EXPECT_DOUBLE_EQ(sg.at(i, j, sg.kEnd() + h),
                                 expectedValue3D(static_cast<long long>(i),
                                                 static_cast<long long>(j),
                                                 static_cast<long long>(sg.kBegin() + h)))
                    << "front halo h=" << h << " i=" << i << " j=" << j;
                EXPECT_DOUBLE_EQ(sg.at(i, j, h),
                                 expectedValue3D(static_cast<long long>(i),
                                                 static_cast<long long>(j),
                                                 static_cast<long long>(sg.kEnd() - hw + h)))
                    << "back halo h=" << h << " i=" << i << " j=" << j;
            }
        }
    }
}

TEST(HaloExchangeTest, MpiNonUniform8RanksHalos3D) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 8) {
        GTEST_SKIP() << "MpiNonUniform8RanksHalos3D requires exactly 8 processes";
    }

    int dims[3] = {2, 2, 2};
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

    int coords[3] = {0, 0, 0};
    MPI_Cart_coords(cart, rank, 3, coords);
    const int offsetX = coords[0] * 2;
    const int offsetY = coords[1] * 2;
    const int offsetZ = coords[2] * 2;

    Subgrid sg(2, 2, 2, 1, cart);
    sg.setNeighbors(left, right, down, up, back, front);

    const long long llOffsetX = static_cast<long long>(offsetX);
    const long long llOffsetY = static_cast<long long>(offsetY);
    const long long llOffsetZ = static_cast<long long>(offsetZ);

    const Real sentinel = -777777.0;
    for (Index k = 0; k < sg.nzTotal(); ++k) {
        for (Index j = 0; j < sg.nyTotal(); ++j) {
            for (Index i = 0; i < sg.nxTotal(); ++i) {
                const bool interior = i >= sg.iBegin() && i < sg.iEnd() &&
                                      j >= sg.jBegin() && j < sg.jEnd() &&
                                      k >= sg.kBegin() && k < sg.kEnd();
                if (interior) {
                    const long long gI = llOffsetX + static_cast<long long>(i) - 1;
                    const long long gJ = llOffsetY + static_cast<long long>(j) - 1;
                    const long long gK = llOffsetZ + static_cast<long long>(k) - 1;
                    sg.at(i, j, k) = expectedValue3D(gI, gJ, gK);
                } else {
                    sg.at(i, j, k) = sentinel;
                }
            }
        }
    }

    PointToPointExchanger ex;
    ex.initialize(sg);
    ex.exchange(sg);

    // Exchanged halo cells must hold the neighbor's interior values; halo
    // cells owned by a physical boundary (PROC_NULL) must keep the sentinel.
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
            const long long gJ = llOffsetY + static_cast<long long>(j) - 1;
            const long long gK = llOffsetZ + static_cast<long long>(k) - 1;
            if (left != MPI_PROC_NULL) {
                EXPECT_DOUBLE_EQ(sg.at(0, j, k), expectedValue3D(llOffsetX - 1, gJ, gK))
                    << "rank " << rank << " left halo j=" << j << " k=" << k;
            } else {
                EXPECT_DOUBLE_EQ(sg.at(0, j, k), sentinel)
                    << "rank " << rank << " left halo j=" << j << " k=" << k;
            }
            if (right != MPI_PROC_NULL) {
                EXPECT_DOUBLE_EQ(sg.at(sg.iEnd(), j, k), expectedValue3D(llOffsetX + 2, gJ, gK))
                    << "rank " << rank << " right halo j=" << j << " k=" << k;
            } else {
                EXPECT_DOUBLE_EQ(sg.at(sg.iEnd(), j, k), sentinel)
                    << "rank " << rank << " right halo j=" << j << " k=" << k;
            }
        }
    }
    for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
        for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
            const long long gI = llOffsetX + static_cast<long long>(i) - 1;
            const long long gK = llOffsetZ + static_cast<long long>(k) - 1;
            if (down != MPI_PROC_NULL) {
                EXPECT_DOUBLE_EQ(sg.at(i, 0, k), expectedValue3D(gI, llOffsetY - 1, gK))
                    << "rank " << rank << " down halo i=" << i << " k=" << k;
            } else {
                EXPECT_DOUBLE_EQ(sg.at(i, 0, k), sentinel)
                    << "rank " << rank << " down halo i=" << i << " k=" << k;
            }
            if (up != MPI_PROC_NULL) {
                EXPECT_DOUBLE_EQ(sg.at(i, sg.jEnd(), k), expectedValue3D(gI, llOffsetY + 2, gK))
                    << "rank " << rank << " up halo i=" << i << " k=" << k;
            } else {
                EXPECT_DOUBLE_EQ(sg.at(i, sg.jEnd(), k), sentinel)
                    << "rank " << rank << " up halo i=" << i << " k=" << k;
            }
        }
    }
    for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            const long long gI = llOffsetX + static_cast<long long>(i) - 1;
            const long long gJ = llOffsetY + static_cast<long long>(j) - 1;
            if (back != MPI_PROC_NULL) {
                EXPECT_DOUBLE_EQ(sg.at(i, j, 0), expectedValue3D(gI, gJ, llOffsetZ - 1))
                    << "rank " << rank << " back halo i=" << i << " j=" << j;
            } else {
                EXPECT_DOUBLE_EQ(sg.at(i, j, 0), sentinel)
                    << "rank " << rank << " back halo i=" << i << " j=" << j;
            }
            if (front != MPI_PROC_NULL) {
                EXPECT_DOUBLE_EQ(sg.at(i, j, sg.kEnd()), expectedValue3D(gI, gJ, llOffsetZ + 2))
                    << "rank " << rank << " front halo i=" << i << " j=" << j;
            } else {
                EXPECT_DOUBLE_EQ(sg.at(i, j, sg.kEnd()), sentinel)
                    << "rank " << rank << " front halo i=" << i << " j=" << j;
            }
        }
    }

    MPI_Comm_free(&cart);
}
