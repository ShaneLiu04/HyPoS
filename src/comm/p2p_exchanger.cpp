#include "comm/halo_exchanger.hpp"
#include "core/exception.hpp"
#include "utils/logger.hpp"
#include <cstring>

namespace hypo {

namespace {

// Face directions; the MPI tag equals the face index that is packed and sent
// toward that direction.
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

// Buffer element counts; each face carries only interior ranges of the
// remaining axes so the receive side can unpack it into its halo band.
Index faceSize(const Subgrid& subgrid, int direction) {
    const Index hw = static_cast<Index>(subgrid.haloWidth());
    switch (direction) {
        case kFaceLeft:
        case kFaceRight:
            return hw * subgrid.nyLocal() * subgrid.nzLocal();
        case kFaceDown:
        case kFaceUp:
            return hw * subgrid.nxLocal() * subgrid.nzLocal();
        case kFaceBack:
        case kFaceFront:
            return hw * subgrid.nxLocal() * subgrid.nyLocal();
        default:
            return 0;
    }
}

// Plane range exchanged by the x/y faces. The 2D case (nzLocal == 1) keeps
// its field in the k=0 plane, while 3D keeps interior data in [kBegin, kEnd).
void faceLayerRange(const Subgrid& subgrid, Index& kFirst, Index& kLast) {
    if (subgrid.nzLocal() == 1) {
        kFirst = 0;
        kLast = 1;
    } else {
        kFirst = subgrid.kBegin();
        kLast = subgrid.kEnd();
    }
}

} // namespace

// ============================================================================
// PointToPointExchanger
// ============================================================================

PointToPointExchanger::~PointToPointExchanger() {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized) {
        return;
    }
    for (MPI_Request& request : activeReqs_) {
        if (request != MPI_REQUEST_NULL) {
            MPI_Request_free(&request);
        }
    }
    activeReqs_.clear();
    activeDirs_.clear();
}

void PointToPointExchanger::initialize(Subgrid& subgrid) {
    comm_ = subgrid.comm();

    activeDirs_.clear();
    activeReqs_.clear();

    for (int direction = 0; direction < kNumDirections; ++direction) {
        const Index size = faceSize(subgrid, direction);
        const std::size_t bufIdx = static_cast<std::size_t>(direction);
        sendBufs_[bufIdx].resize(size);
        recvBufs_[bufIdx].resize(size);

        const int neighbor = neighborForDirection(subgrid, direction);
        if (neighbor == MPI_PROC_NULL) {
            continue;
        }

        MPI_Request sendRequest = MPI_REQUEST_NULL;
        MPI_Request recvRequest = MPI_REQUEST_NULL;
        MPI_Send_init(sendBufs_[bufIdx].data(), static_cast<int>(size), MPI_DOUBLE,
                      neighbor, direction, comm_, &sendRequest);
        MPI_Recv_init(recvBufs_[bufIdx].data(), static_cast<int>(size), MPI_DOUBLE,
                      neighbor, oppositeDirection(direction), comm_, &recvRequest);

        activeDirs_.push_back(direction);
        activeReqs_.push_back(sendRequest);
        activeReqs_.push_back(recvRequest);
    }

    initialized_ = true;
}

void PointToPointExchanger::packFace(const Subgrid& subgrid, int direction,
                                     const Real* u, std::vector<Real>& buf) const {
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    const Index hw = static_cast<Index>(subgrid.haloWidth());
    Index idx = 0;

    switch (direction) {
        case kFaceLeft:
        case kFaceRight: {
            // Layout [k][j][h]: h is contiguous (row-wise memcpy).
            const Index srcCol = (direction == kFaceLeft) ? subgrid.iBegin()
                                                          : (subgrid.iEnd() - hw);
            Index kFirst = 0, kLast = 0;
            faceLayerRange(subgrid, kFirst, kLast);
            for (Index k = kFirst; k < kLast; ++k) {
                for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                    std::memcpy(&buf[idx], &u[(k * nyT + j) * nxT + srcCol], hw * sizeof(Real));
                    idx += hw;
                }
            }
            break;
        }
        case kFaceDown:
        case kFaceUp: {
            // Layout [h][k][i]: i is contiguous.
            const Index srcRow = (direction == kFaceDown) ? subgrid.jBegin()
                                                          : (subgrid.jEnd() - hw);
            const Index ni = subgrid.iEnd() - subgrid.iBegin();
            Index kFirst = 0, kLast = 0;
            faceLayerRange(subgrid, kFirst, kLast);
            for (Index h = 0; h < hw; ++h) {
                for (Index k = kFirst; k < kLast; ++k) {
                    std::memcpy(&buf[idx], &u[(k * nyT + (srcRow + h)) * nxT + subgrid.iBegin()],
                                ni * sizeof(Real));
                    idx += ni;
                }
            }
            break;
        }
        case kFaceBack:
        case kFaceFront: {
            // Layout [h][j][i]: i is contiguous.
            const Index srcLayer = (direction == kFaceBack) ? subgrid.kBegin()
                                                            : (subgrid.kEnd() - hw);
            const Index ni = subgrid.iEnd() - subgrid.iBegin();
            for (Index h = 0; h < hw; ++h) {
                for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                    std::memcpy(&buf[idx], &u[((srcLayer + h) * nyT + j) * nxT + subgrid.iBegin()],
                                ni * sizeof(Real));
                    idx += ni;
                }
            }
            break;
        }
        default:
            break;
    }
    HYPOS_ASSERT(idx == buf.size());
}

void PointToPointExchanger::unpackFace(Subgrid& subgrid, int direction, Real* u,
                                       const std::vector<Real>& buf) const {
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    const Index hw = static_cast<Index>(subgrid.haloWidth());
    Index idx = 0;

    switch (direction) {
        case kFaceLeft:
        case kFaceRight: {
            const Index dstCol = (direction == kFaceLeft) ? 0 : subgrid.iEnd();
            Index kFirst = 0, kLast = 0;
            faceLayerRange(subgrid, kFirst, kLast);
            for (Index k = kFirst; k < kLast; ++k) {
                for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                    std::memcpy(&u[(k * nyT + j) * nxT + dstCol], &buf[idx], hw * sizeof(Real));
                    idx += hw;
                }
            }
            break;
        }
        case kFaceDown:
        case kFaceUp: {
            const Index dstRow = (direction == kFaceDown) ? 0 : subgrid.jEnd();
            const Index ni = subgrid.iEnd() - subgrid.iBegin();
            Index kFirst = 0, kLast = 0;
            faceLayerRange(subgrid, kFirst, kLast);
            for (Index h = 0; h < hw; ++h) {
                for (Index k = kFirst; k < kLast; ++k) {
                    std::memcpy(&u[(k * nyT + (dstRow + h)) * nxT + subgrid.iBegin()],
                                &buf[idx], ni * sizeof(Real));
                    idx += ni;
                }
            }
            break;
        }
        case kFaceBack:
        case kFaceFront: {
            const Index dstLayer = (direction == kFaceBack) ? 0 : subgrid.kEnd();
            const Index ni = subgrid.iEnd() - subgrid.iBegin();
            for (Index h = 0; h < hw; ++h) {
                for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                    std::memcpy(&u[((dstLayer + h) * nyT + j) * nxT + subgrid.iBegin()],
                                &buf[idx], ni * sizeof(Real));
                    idx += ni;
                }
            }
            break;
        }
        default:
            break;
    }
    HYPOS_ASSERT(idx == buf.size());
}

void PointToPointExchanger::exchange(Subgrid& subgrid, Real* data) {
    beginExchange(subgrid, data);
    endExchange(subgrid, data);
}

void PointToPointExchanger::beginExchange(Subgrid& subgrid, Real* data) {
    if (!initialized_ || comm_ == MPI_COMM_NULL || data == nullptr) {
        HYPOS_ERROR("PointToPointExchanger::beginExchange called before initialize "
                    "or with null data buffer; ignored");
        return;
    }

    for (std::size_t i = 0; i < activeDirs_.size(); ++i) {
        const int direction = activeDirs_[i];
        packFace(subgrid, direction, data, sendBufs_[static_cast<std::size_t>(direction)]);
    }

    if (!activeReqs_.empty()) {
        MPI_Startall(static_cast<int>(activeReqs_.size()), activeReqs_.data());
    }
}

void PointToPointExchanger::endExchange(Subgrid& subgrid, Real* data) {
    if (!initialized_ || comm_ == MPI_COMM_NULL || data == nullptr) {
        HYPOS_ERROR("PointToPointExchanger::endExchange called before initialize "
                    "or with null data buffer; ignored");
        return;
    }

    if (!activeReqs_.empty()) {
        MPI_Waitall(static_cast<int>(activeReqs_.size()), activeReqs_.data(),
                    MPI_STATUSES_IGNORE);
    }

    for (int direction : activeDirs_) {
        unpackFace(subgrid, direction, data, recvBufs_[static_cast<std::size_t>(direction)]);
    }
}

} // namespace hypo
