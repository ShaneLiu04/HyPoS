#include "comm/halo_exchanger.hpp"
#include "core/exception.hpp"

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

void PointToPointExchanger::initialize(Subgrid& subgrid) {
    comm_ = subgrid.comm();
    for (int direction = 0; direction < kNumDirections; ++direction) {
        const Index size = faceSize(subgrid, direction);
        const std::size_t bufIdx = static_cast<std::size_t>(direction);
        sendBufs_[bufIdx].resize(size);
        recvBufs_[bufIdx].resize(size);
    }
    requests_.resize(2 * static_cast<std::size_t>(kNumDirections), MPI_REQUEST_NULL);
}

void PointToPointExchanger::packFace(Subgrid& subgrid, int direction, std::vector<Real>& buf) const {
    const Real* u = subgrid.u().data();
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    const Index hw = static_cast<Index>(subgrid.haloWidth());
    Index idx = 0;

    switch (direction) {
        case kFaceLeft:
        case kFaceRight: {
            const Index srcCol = (direction == kFaceLeft) ? subgrid.iBegin()
                                                          : (subgrid.iEnd() - hw);
            Index kFirst = 0, kLast = 0;
            faceLayerRange(subgrid, kFirst, kLast);
            for (Index h = 0; h < hw; ++h) {
                for (Index k = kFirst; k < kLast; ++k) {
                    for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                        buf[idx++] = u[(k * nyT + j) * nxT + (srcCol + h)];
                    }
                }
            }
            break;
        }
        case kFaceDown:
        case kFaceUp: {
            const Index srcRow = (direction == kFaceDown) ? subgrid.jBegin()
                                                          : (subgrid.jEnd() - hw);
            Index kFirst = 0, kLast = 0;
            faceLayerRange(subgrid, kFirst, kLast);
            for (Index h = 0; h < hw; ++h) {
                for (Index k = kFirst; k < kLast; ++k) {
                    for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                        buf[idx++] = u[(k * nyT + (srcRow + h)) * nxT + i];
                    }
                }
            }
            break;
        }
        case kFaceBack:
        case kFaceFront: {
            const Index srcLayer = (direction == kFaceBack) ? subgrid.kBegin()
                                                            : (subgrid.kEnd() - hw);
            for (Index h = 0; h < hw; ++h) {
                for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                    for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                        buf[idx++] = u[((srcLayer + h) * nyT + j) * nxT + i];
                    }
                }
            }
            break;
        }
        default:
            break;
    }
    HYPOS_ASSERT(idx == buf.size());
}

void PointToPointExchanger::unpackFace(Subgrid& subgrid, int direction, const std::vector<Real>& buf) const {
    Real* u = subgrid.u().data();
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
            for (Index h = 0; h < hw; ++h) {
                for (Index k = kFirst; k < kLast; ++k) {
                    for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                        u[(k * nyT + j) * nxT + (dstCol + h)] = buf[idx++];
                    }
                }
            }
            break;
        }
        case kFaceDown:
        case kFaceUp: {
            const Index dstRow = (direction == kFaceDown) ? 0 : subgrid.jEnd();
            Index kFirst = 0, kLast = 0;
            faceLayerRange(subgrid, kFirst, kLast);
            for (Index h = 0; h < hw; ++h) {
                for (Index k = kFirst; k < kLast; ++k) {
                    for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                        u[(k * nyT + (dstRow + h)) * nxT + i] = buf[idx++];
                    }
                }
            }
            break;
        }
        case kFaceBack:
        case kFaceFront: {
            const Index dstLayer = (direction == kFaceBack) ? 0 : subgrid.kEnd();
            for (Index h = 0; h < hw; ++h) {
                for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                    for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                        u[((dstLayer + h) * nyT + j) * nxT + i] = buf[idx++];
                    }
                }
            }
            break;
        }
        default:
            break;
    }
    HYPOS_ASSERT(idx == buf.size());
}

void PointToPointExchanger::exchange(Subgrid& subgrid) {
    beginExchange(subgrid);
    endExchange(subgrid);
}

void PointToPointExchanger::beginExchange(Subgrid& subgrid) {
    for (auto& request : requests_) {
        request = MPI_REQUEST_NULL;
    }

    for (int direction = 0; direction < kNumDirections; ++direction) {
        const int neighbor = neighborForDirection(subgrid, direction);
        if (neighbor == MPI_PROC_NULL) {
            continue;
        }
        const std::size_t bufIdx = static_cast<std::size_t>(direction);

        // Send this direction's own face (tag = direction).
        packFace(subgrid, direction, sendBufs_[bufIdx]);
        MPI_Isend(sendBufs_[bufIdx].data(),
                  static_cast<int>(sendBufs_[bufIdx].size()), MPI_DOUBLE,
                  neighbor, direction, comm_, &requests_[2 * bufIdx]);

        // Receive the face the neighbor sends toward this direction, which is
        // its opposite face (e.g. the left neighbor sends its right face).
        MPI_Irecv(recvBufs_[bufIdx].data(),
                  static_cast<int>(recvBufs_[bufIdx].size()), MPI_DOUBLE,
                  neighbor, oppositeDirection(direction), comm_,
                  &requests_[2 * bufIdx + 1]);
    }
}

void PointToPointExchanger::endExchange(Subgrid& subgrid) {
    MPI_Waitall(static_cast<int>(requests_.size()), requests_.data(), MPI_STATUSES_IGNORE);

    for (int direction = 0; direction < kNumDirections; ++direction) {
        if (neighborForDirection(subgrid, direction) == MPI_PROC_NULL) {
            continue;
        }
        unpackFace(subgrid, direction, recvBufs_[static_cast<std::size_t>(direction)]);
    }
}


} // namespace hypo
