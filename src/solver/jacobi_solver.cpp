#include "solver/solver.hpp"
#include "solver/residual.hpp"
#include "comm/halo_exchanger.hpp"
#include "perf/profiler.hpp"
#include "utils/logger.hpp"
#include "core/exception.hpp"
#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hypo {

JacobiSolver::JacobiSolver(bool overlapComm)
    : overlapComm_(overlapComm) {
}

Real JacobiSolver::updateRegion(Subgrid& subgrid,
                                Index i0, Index i1,
                                Index j0, Index j1,
                                Index k0, Index k1) const {
    if (i0 >= i1 || j0 >= j1 || k0 >= k1) {
        return 0.0;
    }

    Real* __restrict__ u = subgrid.u().data();
    Real* __restrict__ uNext = subgrid.uNext().data();
    const Real* __restrict__ rhs = subgrid.rhs().data();
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();

    Real residual = 0.0;

    if (subgrid.nzLocal() == 1) {
        // 2D: field lives in the k=0 plane
        const Real invDenom = 1.0 / 4.0;
        #pragma omp parallel for reduction(+:residual) schedule(static)
        for (Index j = j0; j < j1; ++j) {
            #pragma omp simd reduction(+:residual)
            for (Index i = i0; i < i1; ++i) {
                const Index idx = j * nxT + i;
                const Real uNew = (u[idx - 1] + u[idx + 1] +
                                   u[idx - nxT] + u[idx + nxT] -
                                   rhs[idx]) * invDenom;
                uNext[idx] = uNew;
                const Real diff = uNew - u[idx];
                residual += diff * diff;
            }
        }
    } else {
        const Real invDenom = 1.0 / 6.0;
        #pragma omp parallel for reduction(+:residual) schedule(static)
        for (Index k = k0; k < k1; ++k) {
            for (Index j = j0; j < j1; ++j) {
                #pragma omp simd reduction(+:residual)
                for (Index i = i0; i < i1; ++i) {
                    const Index idx = (k * nyT + j) * nxT + i;
                    const Real uNew = (u[idx - 1] + u[idx + 1] +
                                       u[idx - nxT] + u[idx + nxT] +
                                       u[idx - nxT * nyT] + u[idx + nxT * nyT] -
                                       rhs[idx]) * invDenom;
                    uNext[idx] = uNew;
                    const Real diff = uNew - u[idx];
                    residual += diff * diff;
                }
            }
        }
    }

    return residual;
}

Real JacobiSolver::iterate(Subgrid& subgrid, HaloExchanger& exchanger) {
    const Int hw = subgrid.haloWidth();

    // Inner box (independent of halo data) and boundary band split.
    // Endpoints are clamped so that degenerate subdomains (nLocal <= 2*hw)
    // still get an exact, non-overlapping partition of the interior.
    const Index i0 = std::min(subgrid.iBegin() + hw, subgrid.iEnd());
    const Index i1 = std::max(subgrid.iEnd() - hw, i0);
    const Index j0 = std::min(subgrid.jBegin() + hw, subgrid.jEnd());
    const Index j1 = std::max(subgrid.jEnd() - hw, j0);

    // Active layer range: [0,1) for the 2D k=0 plane, [kBegin,kEnd) for 3D.
    const bool is2D = (subgrid.nzLocal() == 1);
    const Index kA0 = is2D ? 0 : subgrid.kBegin();
    const Index kA1 = is2D ? 1 : subgrid.kEnd();
    const Index kI0 = is2D ? 0 : std::min(subgrid.kBegin() + hw, kA1);
    const Index kI1 = is2D ? 1 : std::max(subgrid.kEnd() - hw, kI0);

#ifndef NDEBUG
    {
        // Region partition self-check: the seven regions must tile the
        // interior exactly (no overlap, no gap).
        const Index volInterior = (subgrid.iEnd() - subgrid.iBegin()) *
                                  (subgrid.jEnd() - subgrid.jBegin()) *
                                  (kA1 - kA0);
        const Index volInner = (i1 - i0) * (j1 - j0) * (kI1 - kI0);
        const Index volLeft = (i0 - subgrid.iBegin()) * (subgrid.jEnd() - subgrid.jBegin()) * (kA1 - kA0);
        const Index volRight = (subgrid.iEnd() - i1) * (subgrid.jEnd() - subgrid.jBegin()) * (kA1 - kA0);
        const Index volBottom = (i1 - i0) * (j0 - subgrid.jBegin()) * (kA1 - kA0);
        const Index volTop = (i1 - i0) * (subgrid.jEnd() - j1) * (kA1 - kA0);
        const Index volBack = (i1 - i0) * (j1 - j0) * (kI0 - kA0);
        const Index volFront = (i1 - i0) * (j1 - j0) * (kA1 - kI1);
        HYPOS_ASSERT(volInner + volLeft + volRight + volBottom + volTop +
                         volBack + volFront == volInterior);
    }
#endif

    // Residual accumulation order is fixed (inner -> L -> R -> D -> U -> B -> F)
    // and identical for both overlap modes, so residuals are bit-identical.
    Real residual = 0.0;

    if (overlapComm_) {
        {
            HYPOS_PROFILE("halo_exchange");
            exchanger.beginExchange(subgrid);
        }
        {
            HYPOS_PROFILE("stencil_interior");
            residual += updateRegion(subgrid, i0, i1, j0, j1, kI0, kI1);
        }
        {
            HYPOS_PROFILE("halo_wait");
            exchanger.endExchange(subgrid);
        }
    } else {
        {
            HYPOS_PROFILE("halo_exchange");
            exchanger.beginExchange(subgrid);
            exchanger.endExchange(subgrid);
        }
        {
            HYPOS_PROFILE("stencil_interior");
            residual += updateRegion(subgrid, i0, i1, j0, j1, kI0, kI1);
        }
    }

    {
        HYPOS_PROFILE("stencil_boundary");
        residual += updateRegion(subgrid, subgrid.iBegin(), i0, subgrid.jBegin(), subgrid.jEnd(), kA0, kA1);
        residual += updateRegion(subgrid, i1, subgrid.iEnd(), subgrid.jBegin(), subgrid.jEnd(), kA0, kA1);
        residual += updateRegion(subgrid, i0, i1, subgrid.jBegin(), j0, kA0, kA1);
        residual += updateRegion(subgrid, i0, i1, j1, subgrid.jEnd(), kA0, kA1);
        residual += updateRegion(subgrid, i0, i1, j0, j1, kA0, kI0);
        residual += updateRegion(subgrid, i0, i1, j0, j1, kI1, kA1);
    }

    subgrid.swapU();
    subgrid.applyPhysicalBoundary();

    // Exact identity r = D * diff (see class doc): the converted value is
    // the true residual of the state BEFORE this update.
    const Real denom = is2D ? 4.0 : 6.0;
    return denom * std::sqrt(residual);
}

Index JacobiSolver::solve(Subgrid& subgrid,
                          HaloExchanger& exchanger,
                          Index maxIter,
                          Real tolerance) {
    Timer totalTimer;
    totalTimer.start();

    Index completed = 0;
    for (; completed < maxIter; ) {
        HYPOS_PROFILE("jacobi_iteration");
        Real localResidual = iterate(subgrid, exchanger);
        ++completed;

        // Convergence is only evaluated every residualCheckInterval_
        // iterations: each skipped check saves one global reduction.
        if (completed % residualCheckInterval_ == 0) {
            Real globalResidual = 0.0;
            {
                HYPOS_PROFILE("residual_allreduce");
                // localResidual is already the converted true residual of
                // the previous iterate (D * ||diff||), so this reduction
                // directly yields the criterion value.
                Real localSquared = localResidual * localResidual;
                MPI_Allreduce(&localSquared, &globalResidual, 1, MPI_DOUBLE, MPI_SUM, subgrid.comm());
                globalResidual = std::sqrt(globalResidual);
            }
            lastResidual_ = globalResidual;
            notifyProgress(completed);

            if (globalResidual < tolerance) {
                HYPOS_INFO("Converged at iteration " << completed << ", residual = " << globalResidual);
                break;
            }

            if (completed % 500 == 0) {
                HYPOS_INFO("Iteration " << completed << ", residual = " << globalResidual);
            }
        } else {
            notifyProgress(completed);
        }
    }

    // Exit confirmation scan (converged or maxIter): refresh the neighbor
    // halos of the final iterate and report its exact true residual, so
    // lastResidual() is always the residual of the state we return. Skipped
    // when no iteration ran, keeping the "0 before any iteration" contract.
    if (completed > 0) {
        HYPOS_PROFILE("residual_confirm");
        lastResidual_ = globalTrueResidual(subgrid, exchanger);
    }

    totalTimer.stop();
    HYPOS_INFO("Solver finished in " << completed << " iterations, time = " << totalTimer.elapsedSeconds() << " s");

    return completed;
}

} // namespace hypo
