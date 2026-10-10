#pragma once

#include "core/types.hpp"
#include "grid/subgrid.hpp"

namespace hypo {

class HaloExchanger;

/**
 * @brief Local sum of squared true residuals over the interior points:
 *        sum (D*u - sum(neighbors) + rhs)^2 with D = 4 (2D) / 6 (3D),
 *        for the system A u = b where A = D*I - S and b = -rhs.
 *
 * Caller contract: halo cells must be refreshed (halo exchange +
 * applyPhysicalBoundary) before the call — this helper reads halo values
 * as-is and performs no communication, so the result is only meaningful
 * with a synchronized ghost layer. The global reduction is the caller's
 * responsibility, which keeps this kernel unit-testable in isolation.
 *
 * Thread-safe, no exception paths.
 */
Real trueResidualSquaredLocal(const Subgrid& subgrid) noexcept;

/**
 * @brief Global true residual ||r||_2 of the current iterate:
 *        refreshes the neighbor halos of subgrid.u() via the exchanger,
 *        evaluates trueResidualSquaredLocal, and reduces across ranks.
 *
 * Carries no profiler regions: callers wrap it with the appropriate
 * region name (residual_confirm for exit scans, residual_allreduce for
 * in-loop checks). MPI paths do not throw (AGENT_SPEC).
 */
Real globalTrueResidual(Subgrid& subgrid, HaloExchanger& exchanger) noexcept;

/**
 * @brief Field version of the true residual: writes rho = b - A*u into
 *        `residual` over the interior points,
 *        rho = sum(u_neighbors) - D*u - rhs  (b = -rhs, A = D*I - S).
 *
 * Same kernel expression as trueResidualSquaredLocal, returning the field
 * instead of its squared norm (used by multigrid restriction). Halo cells
 * are not written; the caller contract is the same — u halos must already
 * be synchronized. Synchronizing the residual field's own halos (face
 * exchange + corner exchange + physical fill) is the caller's job.
 *
 * Thread-safe, no exception paths.
 */
void residualFieldLocal(const Subgrid& subgrid, Real* residual) noexcept;

} // namespace hypo
