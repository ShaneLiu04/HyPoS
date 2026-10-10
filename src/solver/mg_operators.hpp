#pragma once

#include "core/types.hpp"
#include "grid/subgrid.hpp"

namespace hypo {

/**
 * @brief Coarse-cell ownership: coarse cell I belongs to the rank owning fine
 * point 2I (vertex-coincident convention, AR007 design §4.1/§4.2).
 * Ranges [ceil(offset/2), ceil((offset+nLocal)/2)) tile [0, nGlobal/2)
 * seam-free for any uniform or manual rectangular layout; count may be 0
 * (empty piece, e.g. nLocal == 1 with odd offset).
 */
Index coarseRangeBegin(Index offset, Index nLocal) noexcept;
Index coarseRangeCount(Index offset, Index nLocal) noexcept;

/**
 * @brief Full-weighting restriction of a fine residual field into the local
 * coarse rectangle [coarseX0, coarseX0+cnx) x [coarseY0, coarseY0+cny),
 * packed row-major. Pure kernel (no MPI).
 *
 * Caller contract (design §4.2): fineResidual halos are synchronized — face
 * exchange + exchangeCorners + applyPhysicalBoundary(fineResidual, 0.0) —
 * so every 9-point stencil read lands on a valid cell (interior, neighbor
 * halo, or physical-outside zero).
 */
void restrictResidual(const Subgrid& fine, const Real* fineResidual,
                      Index coarseX0, Index coarseY0, Index cnx, Index cny,
                      Real* coarsePacked) noexcept;

/**
 * @brief Vertex-coincident bilinear prolongation: fine even points copy the
 * coarse parent (weight 1), fine odd points average the two coarse
 * parents (1/2 + 1/2); coarse cells outside the coarse grid read as zero
 * (homogeneous correction boundary, design D3). Accumulates into
 * fine.u() over the fine interior. `coarse` holds the correction in its
 * u() interior with offsets (0, 0) (replicated full coarse grid).
 */
void prolongateCorrection(Subgrid& fine, const Subgrid& coarse) noexcept;

/**
 * @brief Exchange the four diagonal corner halo cells of a 2D field buffer
 * with the diagonal cart neighbors (1 Real per corner, Sendrecv). Physical
 * sides (PROC_NULL face neighbor) are skipped: those corners are already
 * zero from the Dirichlet band fill, and no diagonal rank exists there.
 * MG-specific communication — the face exchangers are not touched (D9).
 */
void exchangeCorners(Subgrid& subgrid, Real* data) noexcept;

} // namespace hypo
