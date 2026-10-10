#include "solver/solver.hpp"
#include "solver/residual.hpp"
#include "comm/halo_exchanger.hpp"
#include "perf/profiler.hpp"
#include "utils/logger.hpp"
#include <cmath>

namespace hypo {

Real RedBlackGSSolver::sweep(Subgrid& subgrid, int parity) const {
    const Index iB = subgrid.iBegin();
    const Index iE = subgrid.iEnd();
    const Index jB = subgrid.jBegin();
    const Index jE = subgrid.jEnd();
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    const long long hw = static_cast<long long>(subgrid.haloWidth());
    Real* u = subgrid.u().data();
    const Real* rhs = subgrid.rhs().data();

    Real residual = 0.0;

    if (subgrid.nzLocal() == 1) {
        const long long tilt = static_cast<long long>(subgrid.offsetX()) +
                               static_cast<long long>(subgrid.offsetY()) - 2LL * hw;
        const Real invDenom = 1.0 / 4.0;
        #pragma omp parallel for reduction(+:residual) schedule(static)
        for (Index j = jB; j < jE; ++j) {
            const long long row = static_cast<long long>(j) + tilt;
            Index i = iB + ((((static_cast<long long>(iB) + row) & 1LL) == parity) ? 0 : 1);
            for (; i < iE; i += 2) {
                const Index idx = j * nxT + i;
                const Real uNew = (u[idx - 1] + u[idx + 1] +
                                   u[idx - nxT] + u[idx + nxT] -
                                   rhs[idx]) * invDenom;
                const Real diff = uNew - u[idx];
                u[idx] = uNew;
                residual += diff * diff;
            }
        }
    } else {
        const long long tilt = static_cast<long long>(subgrid.offsetX()) +
                               static_cast<long long>(subgrid.offsetY()) +
                               static_cast<long long>(subgrid.offsetZ()) - 3LL * hw;
        const Real invDenom = 1.0 / 6.0;
        #pragma omp parallel for reduction(+:residual) schedule(static)
        for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
            for (Index j = jB; j < jE; ++j) {
                const long long row = static_cast<long long>(j) + static_cast<long long>(k) + tilt;
                Index i = iB + ((((static_cast<long long>(iB) + row) & 1LL) == parity) ? 0 : 1);
                for (; i < iE; i += 2) {
                    const Index idx = (k * nyT + j) * nxT + i;
                    const Real uNew = (u[idx - 1] + u[idx + 1] +
                                       u[idx - nxT] + u[idx + nxT] +
                                       u[idx - nxT * nyT] + u[idx + nxT * nyT] -
                                       rhs[idx]) * invDenom;
                    const Real diff = uNew - u[idx];
                    u[idx] = uNew;
                    residual += diff * diff;
                }
            }
        }
    }

    return residual;
}

Real RedBlackGSSolver::iterateCore(Subgrid& subgrid, HaloExchanger& exchanger) const {
    Real updates = 0.0;

    {
        HYPOS_PROFILE("halo_exchange");
        exchanger.exchange(subgrid);
    }
    {
        HYPOS_PROFILE("stencil_interior");
        updates += sweep(subgrid, 0);
    }
    {
        HYPOS_PROFILE("halo_exchange");
        exchanger.exchange(subgrid);
    }
    {
        HYPOS_PROFILE("stencil_boundary");
        updates += sweep(subgrid, 1);
    }

    subgrid.applyPhysicalBoundary();

    return updates;
}

Real RedBlackGSSolver::iterate(Subgrid& subgrid, HaloExchanger& exchanger) {
    // Single-call contract: always report the true residual of the state
    // produced by this iteration (scan = exchange + kernel + reduction).
    iterateCore(subgrid, exchanger);
    HYPOS_PROFILE("residual_allreduce");
    return globalTrueResidual(subgrid, exchanger);
}

void RedBlackGSSolver::smooth(Subgrid& subgrid,
                              HaloExchanger& exchanger,
                              Index sweeps) const {
    for (Index k = 0; k < sweeps; ++k) {
        iterateCore(subgrid, exchanger);
    }
}

Index RedBlackGSSolver::solve(Subgrid& subgrid,
                              HaloExchanger& exchanger,
                              Index maxIter,
                              Real tolerance) {
    Timer totalTimer;
    totalTimer.start();

    Index completed = 0;
    for (; completed < maxIter; ) {
        HYPOS_PROFILE("rbgs_iteration");
        iterateCore(subgrid, exchanger);
        ++completed;

        // In-place updates have no diff-to-residual relation, so the true
        // residual needs a dedicated scan; skipping it on non-check
        // iterations is exactly what the interval buys.
        if (completed % residualCheckInterval_ == 0) {
            Real globalResidual = 0.0;
            {
                HYPOS_PROFILE("residual_allreduce");
                globalResidual = globalTrueResidual(subgrid, exchanger);
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

    // Exit confirmation scan (converged or maxIter): lastResidual() is the
    // true residual of the final iterate. Skipped when no iteration ran,
    // keeping the "0 before any iteration" contract.
    if (completed > 0) {
        HYPOS_PROFILE("residual_confirm");
        lastResidual_ = globalTrueResidual(subgrid, exchanger);
    }

    totalTimer.stop();
    HYPOS_INFO("RedBlackGS finished in " << completed << " iterations, time = "
                                         << totalTimer.elapsedSeconds() << " s");

    return completed;
}

} // namespace hypo
