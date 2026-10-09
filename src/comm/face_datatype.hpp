#pragma once

#include "grid/subgrid.hpp"
#include <mpi.h>

namespace hypo {

// 面方向枚举与 tag 约定（与 p2p_exchanger.cpp 内部约定一致）：
// tag = 被打包并发往该方向的面索引。
enum FaceDirection : int {
    kFaceLeft = 0,
    kFaceRight = 1,
    kFaceDown = 2,
    kFaceUp = 3,
    kFaceBack = 4,
    kFaceFront = 5
};

inline int neighborForDirection(const Subgrid& subgrid, int direction) {
    switch (direction) {
        case kFaceLeft: return subgrid.neighborLeft();
        case kFaceRight: return subgrid.neighborRight();
        case kFaceDown: return subgrid.neighborDown();
        case kFaceUp: return subgrid.neighborUp();
        case kFaceBack: return subgrid.neighborBack();
        case kFaceFront: return subgrid.neighborFront();
        default: return MPI_PROC_NULL;
    }
}

constexpr int oppositeDirection(int direction) {
    switch (direction) {
        case kFaceLeft: return kFaceRight;
        case kFaceRight: return kFaceLeft;
        case kFaceDown: return kFaceUp;
        case kFaceUp: return kFaceDown;
        case kFaceBack: return kFaceFront;
        case kFaceFront: return kFaceBack;
        default: return direction;
    }
}

// 与 p2p faceLayerRange 相同的 k 域约定：2D 场在 k=0 平面，3D 场在 [kBegin,kEnd)。
inline void faceLayerRange(const Subgrid& subgrid, Index& kFirst, Index& kLast) {
    if (subgrid.nzLocal() == 1) {
        kFirst = 0;
        kLast = 1;
    } else {
        kFirst = subgrid.kBegin();
        kLast = subgrid.kEnd();
    }
}

// 为 direction 面构造一个 subarray：send 侧描述 interior 发送带，recv 侧描述
// halo 接收带；块在 padded 缓冲内的绝对位置由型自述（接收直落，无 unpack）。
// datatype 直传与 collective alltoallw 共用；alltoallw 路径下 displs 恒 0
// （design 第 3 轮 P1：定位职责归型，双重偏移禁用）。
inline MPI_Datatype makeFaceType(const Subgrid& subgrid, int direction, bool recvSide) {
    const int nzT = static_cast<int>(subgrid.nzTotal());
    const int nyT = static_cast<int>(subgrid.nyTotal());
    const int nxT = static_cast<int>(subgrid.nxTotal());
    const int hw = static_cast<int>(subgrid.haloWidth());
    const int niL = static_cast<int>(subgrid.nxLocal());
    const int nyL = static_cast<int>(subgrid.nyLocal());
    const int iB = static_cast<int>(subgrid.iBegin());
    const int iE = static_cast<int>(subgrid.iEnd());
    const int jB = static_cast<int>(subgrid.jBegin());
    const int jE = static_cast<int>(subgrid.jEnd());
    const int kB = static_cast<int>(subgrid.kBegin());
    const int kE = static_cast<int>(subgrid.kEnd());

    Index kFirst = 0, kLast = 0;
    faceLayerRange(subgrid, kFirst, kLast);
    const int k0 = static_cast<int>(kFirst);
    const int kCount = static_cast<int>(kLast - kFirst);

    int sizes[3] = {nzT, nyT, nxT};
    int subsizes[3] = {0, 0, 0};
    int starts[3] = {0, 0, 0};

    switch (direction) {
        case kFaceLeft:
        case kFaceRight: {
            subsizes[0] = kCount;
            subsizes[1] = nyL;
            subsizes[2] = hw;
            starts[0] = k0;
            starts[1] = jB;
            if (direction == kFaceLeft) {
                starts[2] = recvSide ? 0 : iB;
            } else {
                starts[2] = recvSide ? iE : (iE - hw);
            }
            break;
        }
        case kFaceDown:
        case kFaceUp: {
            subsizes[0] = kCount;
            subsizes[1] = hw;
            subsizes[2] = niL;
            starts[0] = k0;
            if (direction == kFaceDown) {
                starts[1] = recvSide ? 0 : jB;
            } else {
                starts[1] = recvSide ? jE : (jE - hw);
            }
            starts[2] = iB;
            break;
        }
        case kFaceBack:
        case kFaceFront: {
            subsizes[0] = hw;
            subsizes[1] = nyL;
            subsizes[2] = niL;
            if (direction == kFaceBack) {
                starts[0] = recvSide ? 0 : kB;
            } else {
                starts[0] = recvSide ? kE : (kE - hw);
            }
            starts[1] = jB;
            starts[2] = iB;
            break;
        }
        default:
            break;
    }

    MPI_Datatype type = MPI_DATATYPE_NULL;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
    return type;
}

} // namespace hypo
