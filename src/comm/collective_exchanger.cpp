#include "comm/halo_exchanger.hpp"
#include "comm/face_datatype.hpp"
#include "core/exception.hpp"
#include "utils/logger.hpp"

namespace hypo {

// ============================================================================
// CollectiveExchanger
// ============================================================================

CollectiveExchanger::CollectiveExchanger() = default;

CollectiveExchanger::~CollectiveExchanger() {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized) {
        return;
    }
    releaseResources();
}

void CollectiveExchanger::releaseResources() noexcept {
    if (pendingReq_ != MPI_REQUEST_NULL) {
        MPI_Request_free(&pendingReq_);
        pendingReq_ = MPI_REQUEST_NULL;
    }
    for (MPI_Datatype& type : sendTypes_) {
        if (type != MPI_DATATYPE_NULL) {
            MPI_Type_free(&type);
        }
    }
    for (MPI_Datatype& type : recvTypes_) {
        if (type != MPI_DATATYPE_NULL) {
            MPI_Type_free(&type);
        }
    }
    sendTypes_.clear();
    recvTypes_.clear();
    sendCounts_.clear();
    recvCounts_.clear();
    sendDispls_.clear();
    recvDispls_.clear();
    activeDirs_.clear();
    if (graphComm_ != MPI_COMM_NULL) {
        MPI_Comm_free(&graphComm_);
        graphComm_ = MPI_COMM_NULL;
    }
}

void CollectiveExchanger::initialize(Subgrid& subgrid) {
    comm_ = subgrid.comm();

    releaseResources();

    int myRank = 0;
    MPI_Comm_rank(comm_, &myRank);

    // 邻居图（D4）：sources 对称构造（= destinations 同序同值）。入邻居多重集
    // 恒等——共享面被两端各自声明；非周期笛卡尔边界 rank 的 activeDirs 不对合
    // 自反封闭也不影响。recv 块分情形落位在下方建型处。
    std::vector<int> destinations;
    std::vector<int> sources;
    for (int direction = 0; direction < 6; ++direction) {
        const int neighbor = neighborForDirection(subgrid, direction);
        if (neighbor == MPI_PROC_NULL) {
            continue;
        }
        activeDirs_.push_back(direction);
        destinations.push_back(neighbor);
        sources.push_back(neighbor);
    }

    const int nEdges = static_cast<int>(activeDirs_.size());
    MPI_Dist_graph_create_adjacent(comm_, nEdges, sources.data(), MPI_UNWEIGHTED,
                                   nEdges, destinations.data(), MPI_UNWEIGHTED,
                                   MPI_INFO_NULL, 0 /* reorder */, &graphComm_);

    // alltoallw 参数（按 activeDirs_ 枚举序）：块位置由面型绝对定位自述，
    // displs 恒 0（design 第 3 轮 P1：双重偏移禁用）。
    const std::size_t m = activeDirs_.size();
    sendTypes_.resize(m, MPI_DATATYPE_NULL);
    recvTypes_.resize(m, MPI_DATATYPE_NULL);
    sendCounts_.assign(m, 1);
    recvCounts_.assign(m, 1);
    sendDispls_.assign(m, 0);
    recvDispls_.assign(m, 0);

    for (std::size_t i = 0; i < m; ++i) {
        const int dir = activeDirs_[i];
        sendTypes_[i] = makeFaceType(subgrid, dir, false);
        MPI_Type_commit(&sendTypes_[i]);

        // 分情形落位（D4）：唯一边（邻居 != 自身）收对端 interior(opposite(dir))
        // 落 halo(dir)；自环平行边（邻居 == 自身，非周期笛卡尔即该维 np=1）按
        // 出现序收到自己 interior(dir)，落 halo(opposite(dir))。
        const int neighbor = neighborForDirection(subgrid, dir);
        const int recvFace =
            (neighbor == myRank) ? oppositeDirection(dir) : dir;
        recvTypes_[i] = makeFaceType(subgrid, recvFace, true);
        MPI_Type_commit(&recvTypes_[i]);
    }

    initialized_ = true;
}

void CollectiveExchanger::exchange(Subgrid& subgrid, Real* data) {
    beginExchange(subgrid, data);
    endExchange(subgrid, data);
}

void CollectiveExchanger::beginExchange(Subgrid& subgrid, Real* data) {
    if (!initialized_ || graphComm_ == MPI_COMM_NULL || data == nullptr) {
        HYPOS_ERROR("CollectiveExchanger::beginExchange called before initialize "
                    "or with null data buffer; ignored");
        return;
    }
    if (pendingReq_ != MPI_REQUEST_NULL) {
        HYPOS_ERROR("CollectiveExchanger::beginExchange called before the previous "
                    "endExchange completed; ignored");
        return;
    }

    MPI_Ineighbor_alltoallw(data, sendCounts_.data(), sendDispls_.data(),
                            sendTypes_.data(), data, recvCounts_.data(),
                            recvDispls_.data(), recvTypes_.data(), graphComm_,
                            &pendingReq_);
}

void CollectiveExchanger::endExchange(Subgrid& subgrid, Real* data) {
    if (!initialized_ || graphComm_ == MPI_COMM_NULL || data == nullptr) {
        HYPOS_ERROR("CollectiveExchanger::endExchange called before initialize "
                    "or with null data buffer; ignored");
        return;
    }

    if (pendingReq_ != MPI_REQUEST_NULL) {
        MPI_Wait(&pendingReq_, MPI_STATUS_IGNORE);
        pendingReq_ = MPI_REQUEST_NULL;
    }
}

} // namespace hypo
