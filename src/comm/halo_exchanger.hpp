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
 */
class HaloExchanger {
public:
    virtual ~HaloExchanger() = default;

    /**
     * @brief Initialize the exchanger for the given subgrid.
     */
    virtual void initialize(Subgrid& subgrid) = 0;

    /**
     * @brief Perform a halo exchange.
     * After this call, halo cells in subgrid contain neighbor data.
     * @param subgrid The local subdomain (u array is read/written).
     */
    virtual void exchange(Subgrid& subgrid) = 0;

    /**
     * @brief Non-blocking start of halo exchange.
     * Use with overlapComm=true in solver.
     */
    virtual void beginExchange(Subgrid& subgrid) = 0;

    /**
     * @brief Complete a non-blocking exchange.
     */
    virtual void endExchange(Subgrid& subgrid) = 0;

    virtual std::string name() const = 0;
};

/**
 * @brief Non-blocking point-to-point halo exchange using MPI_Isend/Irecv.
 * Supports 6 directions (left/right/down/up/back/front); each direction
 * exchanges one interior face layer of width haloWidth. Tags identify the
 * face being sent (0=left, 1=right, 2=down, 3=up, 4=back, 5=front).
 */
class PointToPointExchanger : public HaloExchanger {
public:
    void initialize(Subgrid& subgrid) override;
    void exchange(Subgrid& subgrid) override;
    void beginExchange(Subgrid& subgrid) override;
    void endExchange(Subgrid& subgrid) override;
    std::string name() const override { return "p2p"; }

private:
    static constexpr int kNumDirections = 6;

    MPI_Comm comm_ = MPI_COMM_NULL;
    std::vector<MPI_Request> requests_;
    std::array<std::vector<Real>, kNumDirections> sendBufs_;
    std::array<std::vector<Real>, kNumDirections> recvBufs_;

    void packFace(Subgrid& subgrid, int direction, std::vector<Real>& buf) const;
    void unpackFace(Subgrid& subgrid, int direction, const std::vector<Real>& buf) const;
};

/**
 * @brief Collective halo exchange placeholder.
 * In this version it delegates to PointToPointExchanger; a true
 * MPI_Neighbor_allgatherv implementation remains a roadmap item.
 */
class CollectiveExchanger : public HaloExchanger {
public:
    void initialize(Subgrid& subgrid) override;
    void exchange(Subgrid& subgrid) override;
    void beginExchange(Subgrid& subgrid) override;
    void endExchange(Subgrid& subgrid) override;
    std::string name() const override { return "collective"; }

private:
    PointToPointExchanger delegate_;
};

} // namespace hypo
