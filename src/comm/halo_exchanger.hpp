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
 * @brief Collective halo exchange placeholder.
 * In this version it delegates to PointToPointExchanger; a true
 * MPI_Neighbor_allgatherv implementation remains a roadmap item.
 */
class CollectiveExchanger : public HaloExchanger {
public:
    using HaloExchanger::exchange;
    using HaloExchanger::beginExchange;
    using HaloExchanger::endExchange;

    void initialize(Subgrid& subgrid) override;
    void exchange(Subgrid& subgrid, Real* data) override;
    void beginExchange(Subgrid& subgrid, Real* data) override;
    void endExchange(Subgrid& subgrid, Real* data) override;
    std::string name() const override { return "collective"; }

private:
    PointToPointExchanger delegate_;
};

} // namespace hypo
