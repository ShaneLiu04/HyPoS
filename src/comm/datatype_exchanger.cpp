#include "comm/halo_exchanger.hpp"
#include "core/exception.hpp"
#include "utils/logger.hpp"

namespace hypo {

namespace {

// 与 p2p_exchanger.cpp 相同的方向枚举与 tag 约定。
enum FaceDirection : int {
    kFaceLeft = 0,
    kFaceRight = 1,
    kFaceDown = 2,
    kFaceUp = 3,
    kFaceBack = 4,
    kFaceFront = 5
};

int neighborForDirection(const Subgrid& subgrid, int direction) {
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
void faceLayerRange(const Subgrid& subgrid, Index& kFirst, Index& kLast) {
    if (subgrid.nzLocal() == 1) {
        kFirst = 0;
        kLast = 1;
    } else {
        kFirst = subgrid.kBegin();
        kLast = subgrid.kEnd();
    }
}

// 为 direction 面构造一个 subarray：send 侧描述 interior 发送带，recv 侧描述
// halo 接收带；块在 padded 缓冲内的绝对位置由型自述（recv 直落，无 unpack）。
MPI_Datatype makeFaceType(const Subgrid& subgrid, int direction, bool recvSide) {
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

} // namespace

// ============================================================================
// DatatypeExchanger
// ============================================================================

DatatypeExchanger::DatatypeExchanger() {
    sendTypes_.fill(MPI_DATATYPE_NULL);
    recvTypes_.fill(MPI_DATATYPE_NULL);
}

DatatypeExchanger::~DatatypeExchanger() {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized) {
        return;
    }
    releaseResources();
}

void DatatypeExchanger::releaseResources() noexcept {
    for (MPI_Request& request : pendingReqs_) {
        if (request != MPI_REQUEST_NULL) {
            MPI_Request_free(&request);
        }
    }
    pendingReqs_.clear();
    activeDirs_.clear();
    for (MPI_Datatype& type : sendTypes_) {
        if (type != MPI_DATATYPE_NULL) {
            MPI_Type_free(&type);
            type = MPI_DATATYPE_NULL;
        }
    }
    for (MPI_Datatype& type : recvTypes_) {
        if (type != MPI_DATATYPE_NULL) {
            MPI_Type_free(&type);
            type = MPI_DATATYPE_NULL;
        }
    }
}

void DatatypeExchanger::initialize(Subgrid& subgrid) {
    comm_ = subgrid.comm();

    releaseResources();

    for (int direction = 0; direction < kNumDirections; ++direction) {
        if (neighborForDirection(subgrid, direction) == MPI_PROC_NULL) {
            continue;
        }
        const std::size_t idx = static_cast<std::size_t>(direction);
        sendTypes_[idx] = makeFaceType(subgrid, direction, false);
        recvTypes_[idx] = makeFaceType(subgrid, direction, true);
        MPI_Type_commit(&sendTypes_[idx]);
        MPI_Type_commit(&recvTypes_[idx]);
        activeDirs_.push_back(direction);
    }

    initialized_ = true;
}

void DatatypeExchanger::exchange(Subgrid& subgrid, Real* data) {
    beginExchange(subgrid, data);
    endExchange(subgrid, data);
}

void DatatypeExchanger::beginExchange(Subgrid& subgrid, Real* data) {
    if (!initialized_ || comm_ == MPI_COMM_NULL || data == nullptr) {
        HYPOS_ERROR("DatatypeExchanger::beginExchange called before initialize "
                    "or with null data buffer; ignored");
        return;
    }
    if (!pendingReqs_.empty()) {
        HYPOS_ERROR("DatatypeExchanger::beginExchange called before the previous "
                    "endExchange completed; ignored");
        return;
    }

    pendingReqs_.resize(2 * activeDirs_.size());
    std::size_t slot = 0;
    for (const int direction : activeDirs_) {
        MPI_Irecv(data, 1, recvTypes_[static_cast<std::size_t>(direction)],
                  neighborForDirection(subgrid, direction), oppositeDirection(direction),
                  comm_, &pendingReqs_[slot++]);
    }
    for (const int direction : activeDirs_) {
        MPI_Isend(data, 1, sendTypes_[static_cast<std::size_t>(direction)],
                  neighborForDirection(subgrid, direction), direction,
                  comm_, &pendingReqs_[slot++]);
    }
}

void DatatypeExchanger::endExchange(Subgrid& subgrid, Real* data) {
    if (!initialized_ || comm_ == MPI_COMM_NULL || data == nullptr) {
        HYPOS_ERROR("DatatypeExchanger::endExchange called before initialize "
                    "or with null data buffer; ignored");
        return;
    }

    if (!pendingReqs_.empty()) {
        MPI_Waitall(static_cast<int>(pendingReqs_.size()), pendingReqs_.data(),
                    MPI_STATUSES_IGNORE);
        pendingReqs_.clear();
    }
}

} // namespace hypo
