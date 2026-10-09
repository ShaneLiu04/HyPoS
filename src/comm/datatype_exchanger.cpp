#include "comm/halo_exchanger.hpp"
#include "comm/face_datatype.hpp"
#include "core/exception.hpp"
#include "utils/logger.hpp"

namespace hypo {

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
