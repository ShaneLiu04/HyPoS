#pragma once

#include "core/types.hpp"
#include "grid/subgrid.hpp"
#include <mpi.h>
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
 */
class PointToPointExchanger : public HaloExchanger {
public:
    void initialize(Subgrid& subgrid) override;
    void exchange(Subgrid& subgrid) override;
    void beginExchange(Subgrid& subgrid) override;
    void endExchange(Subgrid& subgrid) override;
    std::string name() const override { return "p2p"; }

private:
    MPI_Comm comm_ = MPI_COMM_NULL;
    std::vector<MPI_Request> requests_;
    std::vector<Real> sendBufLeft_;
    std::vector<Real> sendBufRight_;
    std::vector<Real> sendBufDown_;
    std::vector<Real> sendBufUp_;
    std::vector<Real> recvBufLeft_;
    std::vector<Real> recvBufRight_;
    std::vector<Real> recvBufDown_;
    std::vector<Real> recvBufUp_;

    void packSend(Subgrid& subgrid, int direction, std::vector<Real>& buf) const;
    void unpackRecv(Subgrid& subgrid, int direction, const std::vector<Real>& buf) const;
};

/**
 * @brief Collective halo exchange using MPI_Neighbor_allgather.
 */
class CollectiveExchanger : public HaloExchanger {
public:
    void initialize(Subgrid& subgrid) override;
    void exchange(Subgrid& subgrid) override;
    void beginExchange(Subgrid& subgrid) override;
    void endExchange(Subgrid& subgrid) override;
    std::string name() const override { return "collective"; }

private:
    MPI_Comm comm_ = MPI_COMM_NULL;
    MPI_Request request_ = MPI_REQUEST_NULL;
    std::vector<Real> sendBuf_;
    std::vector<Real> recvBuf_;
    std::vector<int> sendCounts_;
    std::vector<int> recvCounts_;
    std::vector<int> sendDispls_;
    std::vector<int> recvDispls_;
};

} // namespace hypo
