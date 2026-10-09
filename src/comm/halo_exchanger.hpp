#pragma once

#include "core/types.hpp"
#include "grid/subgrid.hpp"
#include <mpi.h>
#include <array>
#include <string>
#include <vector>

namespace hypo {

/**
 * @brief Abstract interface for halo (ghost cell) exchange.
 * Implementations handle boundary synchronization using different MPI patterns.
 * The primary entry points take an explicit data buffer with the subgrid
 * padded layout; the convenience overloads default to subgrid.u().
 */
class HaloExchanger {
public:
    virtual ~HaloExchanger() = default;

    /**
     * @brief Initialize the exchanger for the given subgrid.
     */
    virtual void initialize(Subgrid& subgrid) = 0;

    /**
     * @brief Perform a halo exchange on the given buffer.
     * After this call, halo cells in `data` contain neighbor interior data.
     */
    virtual void exchange(Subgrid& subgrid, Real* data) = 0;

    /**
     * @brief Non-blocking start of halo exchange (use with overlapComm).
     */
    virtual void beginExchange(Subgrid& subgrid, Real* data) = 0;

    /**
     * @brief Complete a non-blocking exchange.
     */
    virtual void endExchange(Subgrid& subgrid, Real* data) = 0;

    // Convenience overloads operating on subgrid.u().
    void exchange(Subgrid& subgrid) { exchange(subgrid, subgrid.u().data()); }
    void beginExchange(Subgrid& subgrid) { beginExchange(subgrid, subgrid.u().data()); }
    void endExchange(Subgrid& subgrid) { endExchange(subgrid, subgrid.u().data()); }

    virtual std::string name() const = 0;
};

/**
 * @brief Non-blocking point-to-point halo exchange using persistent requests.
 * 6 directions (left/right/down/up/back/front); each direction exchanges one
 * interior face layer of width haloWidth. Tags identify the face being sent
 * (0=left, 1=right, 2=down, 3=up, 4=back, 5=front).
 */
class PointToPointExchanger : public HaloExchanger {
public:
    using HaloExchanger::exchange;
    using HaloExchanger::beginExchange;
    using HaloExchanger::endExchange;

    ~PointToPointExchanger() override;

    void initialize(Subgrid& subgrid) override;
    void exchange(Subgrid& subgrid, Real* data) override;
    void beginExchange(Subgrid& subgrid, Real* data) override;
    void endExchange(Subgrid& subgrid, Real* data) override;
    std::string name() const override { return "p2p"; }

private:
    static constexpr int kNumDirections = 6;

    MPI_Comm comm_ = MPI_COMM_NULL;
    bool initialized_ = false;
    std::vector<int> activeDirs_;
    std::vector<MPI_Request> activeReqs_;
    std::array<std::vector<Real>, kNumDirections> sendBufs_;
    std::array<std::vector<Real>, kNumDirections> recvBufs_;

    void packFace(const Subgrid& subgrid, int direction, const Real* u, std::vector<Real>& buf) const;
    void unpackFace(Subgrid& subgrid, int direction, Real* u, const std::vector<Real>& buf) const;
};

/**
 * @brief Halo exchange via MPI derived datatypes (direct placement).
 * Each active face gets one MPI_Type_create_subarray describing the interior
 * send strip and the halo receive strip inside the padded buffer, so data
 * moves without pack/unpack staging buffers. Requests are non-persistent
 * because the data pointer may differ between exchanges (e.g. CG's p buffer).
 * Tags match the p2p convention: send tag = face index, receive tag =
 * opposite face index.
 */
class DatatypeExchanger : public HaloExchanger {
public:
    using HaloExchanger::exchange;
    using HaloExchanger::beginExchange;
    using HaloExchanger::endExchange;

    DatatypeExchanger();
    ~DatatypeExchanger() override;

    void initialize(Subgrid& subgrid) override;
    void exchange(Subgrid& subgrid, Real* data) override;
    void beginExchange(Subgrid& subgrid, Real* data) override;
    void endExchange(Subgrid& subgrid, Real* data) override;
    std::string name() const override { return "datatype"; }

private:
    static constexpr int kNumDirections = 6;

    MPI_Comm comm_ = MPI_COMM_NULL;
    bool initialized_ = false;
    std::vector<int> activeDirs_;
    std::vector<MPI_Request> pendingReqs_;
    std::array<MPI_Datatype, kNumDirections> sendTypes_;
    std::array<MPI_Datatype, kNumDirections> recvTypes_;

    void releaseResources() noexcept;
};

/**
 * @brief True collective halo exchange via a distributed graph communicator.
 * Builds the neighbor graph from the six face directions and moves face data
 * with MPI_(I)Neighbor_alltoallw using the same face subarray types as
 * DatatypeExchanger (block positions are self-described by the types, so all
 * displacements are zero). Receive blocks land case-wise: halo(d) for a
 * unique edge to a real neighbor, halo(opposite(d)) for a self edge (the
 * k-th occurrence pairs by order).
 */
class CollectiveExchanger : public HaloExchanger {
public:
    using HaloExchanger::exchange;
    using HaloExchanger::beginExchange;
    using HaloExchanger::endExchange;

    CollectiveExchanger();
    ~CollectiveExchanger() override;

    void initialize(Subgrid& subgrid) override;
    void exchange(Subgrid& subgrid, Real* data) override;
    void beginExchange(Subgrid& subgrid, Real* data) override;
    void endExchange(Subgrid& subgrid, Real* data) override;
    std::string name() const override { return "collective"; }

private:
    MPI_Comm comm_ = MPI_COMM_NULL;
    MPI_Comm graphComm_ = MPI_COMM_NULL;
    bool initialized_ = false;
    std::vector<int> activeDirs_;
    std::vector<MPI_Datatype> sendTypes_;
    std::vector<MPI_Datatype> recvTypes_;
    std::vector<int> sendCounts_;
    std::vector<int> recvCounts_;
    std::vector<MPI_Aint> sendDispls_;
    std::vector<MPI_Aint> recvDispls_;
    MPI_Request pendingReq_ = MPI_REQUEST_NULL;

    void releaseResources() noexcept;
};

} // namespace hypo
