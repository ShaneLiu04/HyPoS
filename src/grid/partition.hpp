#pragma once

#include "core/types.hpp"
#include "grid/grid.hpp"
#include <mpi.h>
#include <string>
#include <vector>

namespace hypo {

struct SubgridInfo {
    Index nxLocal = 0;
    Index nyLocal = 0;
    Index nzLocal = 0;
    Index offsetX = 0;  ///< Global offset of this subdomain in x
    Index offsetY = 0;  ///< Global offset of this subdomain in y
    Index offsetZ = 0;  ///< Global offset of this subdomain in z
    int   rank = 0;
    /// Cartesian communicator created by the partitioner. Ownership is
    /// TRANSFERRED to the caller: the caller must MPI_Comm_free it.
    /// MPI_COMM_NULL when no topology was created.
    MPI_Comm cartComm = MPI_COMM_NULL;
};

/**
 * @brief Abstract interface for grid partitioning strategies.
 */
class GridPartition {
public:
    virtual ~GridPartition() = default;

    /**
     * @brief Given global grid and MPI communicator, compute local subdomain info.
     */
    virtual SubgridInfo partition(const Grid& grid, MPI_Comm comm, int rank) const = 0;

    /**
     * @brief Return the Cartesian topology dimensions used.
     */
    virtual void getTopologyDims(int& dimsX, int& dimsY, int& dimsZ) const = 0;

    /**
     * @brief Human-readable name of this partition strategy.
     */
    virtual std::string name() const = 0;
};

/**
 * @brief Uniform block partitioning using MPI_Cart_create.
 * Divides the grid as evenly as possible among processes.
 */
class UniformPartition : public GridPartition {
public:
    explicit UniformPartition(int dimsX = 0, int dimsY = 0, int dimsZ = 0);

    SubgridInfo partition(const Grid& grid, MPI_Comm comm, int rank) const override;
    void getTopologyDims(int& dimsX, int& dimsY, int& dimsZ) const override;
    std::string name() const override { return "uniform"; }

private:
    int dimsX_ = 0;
    int dimsY_ = 0;
    int dimsZ_ = 0;
};

} // namespace hypo
