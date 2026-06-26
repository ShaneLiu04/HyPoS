#pragma once

#include "core/types.hpp"
#include "core/aligned_buffer.hpp"
#include "grid/grid.hpp"
#include <mpi.h>

namespace hypo {

/**
 * @brief Local subdomain descriptor with halo layers.
 * Manages the local portion of the global grid, including ghost/halo cells.
 */
class Subgrid {
public:
    Subgrid() = default;

    Subgrid(Index nxLocal, Index nyLocal, Index nzLocal, Int haloWidth, MPI_Comm comm);

    // Dimensions
    Index nxLocal() const noexcept { return nxLocal_; }
    Index nyLocal() const noexcept { return nyLocal_; }
    Index nzLocal() const noexcept { return nzLocal_; }
    Int   haloWidth() const noexcept { return haloWidth_; }

    // Total padded dimensions (including halo)
    Index nxTotal() const noexcept { return nxLocal_ + 2 * haloWidth_; }
    Index nyTotal() const noexcept { return nyLocal_ + 2 * haloWidth_; }
    Index nzTotal() const noexcept { return nzLocal_ + 2 * haloWidth_; }
    Index totalCells() const noexcept { return nxTotal() * nyTotal() * nzTotal(); }

    // Interior index range (excluding halo)
    Index iBegin() const noexcept { return haloWidth_; }
    Index iEnd() const noexcept { return nxLocal_ + haloWidth_; }
    Index jBegin() const noexcept { return haloWidth_; }
    Index jEnd() const noexcept { return nyLocal_ + haloWidth_; }
    Index kBegin() const noexcept { return haloWidth_; }
    Index kEnd() const noexcept { return nzLocal_ + haloWidth_; }

    // Flat index for padded array
    Index index(Index i, Index j, Index k = 0) const noexcept {
        return (k * nyTotal() + j) * nxTotal() + i;
    }

    // Data arrays (SoA layout)
    AlignedBuffer<Real>& u() noexcept { return u_; }
    const AlignedBuffer<Real>& u() const noexcept { return u_; }
    AlignedBuffer<Real>& uNext() noexcept { return uNext_; }
    const AlignedBuffer<Real>& uNext() const noexcept { return uNext_; }
    AlignedBuffer<Real>& rhs() noexcept { return rhs_; }
    const AlignedBuffer<Real>& rhs() const noexcept { return rhs_; }

    // Neighbor ranks in Cartesian topology
    int neighborLeft() const noexcept { return neighborLeft_; }
    int neighborRight() const noexcept { return neighborRight_; }
    int neighborDown() const noexcept { return neighborDown_; }
    int neighborUp() const noexcept { return neighborUp_; }
    int neighborBack() const noexcept { return neighborBack_; }
    int neighborFront() const noexcept { return neighborFront_; }

    void setNeighbors(int left, int right, int down, int up, int back = MPI_PROC_NULL, int front = MPI_PROC_NULL) noexcept;

    MPI_Comm comm() const noexcept { return comm_; }

    /**
     * @brief Access a cell value (from u array).
     */
    Real& at(Index i, Index j, Index k = 0) noexcept { return u_[index(i, j, k)]; }
    Real at(Index i, Index j, Index k = 0) const noexcept { return u_[index(i, j, k)]; }

    /**
     * @brief Initialize all arrays to zero.
     */
    void zeroInitialize() noexcept;

    /**
     * @brief Apply Dirichlet boundary conditions to halo cells.
     * For simplicity, sets halo to fixedValue.
     */
    void applyDirichletBC(Real fixedValue) noexcept;

    /**
     * @brief Apply Neumann (zero-flux) BC to halo cells.
     */
    void applyNeumannBC() noexcept;

private:
    Index nxLocal_ = 0;
    Index nyLocal_ = 0;
    Index nzLocal_ = 0;
    Int haloWidth_ = 1;

    AlignedBuffer<Real> u_;
    AlignedBuffer<Real> uNext_;
    AlignedBuffer<Real> rhs_;

    int neighborLeft_   = MPI_PROC_NULL;
    int neighborRight_  = MPI_PROC_NULL;
    int neighborDown_   = MPI_PROC_NULL;
    int neighborUp_     = MPI_PROC_NULL;
    int neighborBack_   = MPI_PROC_NULL;
    int neighborFront_  = MPI_PROC_NULL;

    MPI_Comm comm_ = MPI_COMM_NULL;
};

} // namespace hypo
