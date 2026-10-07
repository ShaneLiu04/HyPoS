#include "solver/solver.hpp"
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

Real RedBlackGSSolver::iterate(Subgrid& subgrid, HaloExchanger& exchanger) {
    Real residual = 0.0;

    {
        HYPOS_PROFILE("halo_exchange");
        exchanger.exchange(subgrid);
    }
    {
        HYPOS_PROFILE("stencil_interior");
        residual += sweep(subgrid, 0);
    }
    {
        HYPOS_PROFILE("halo_exchange");
        exchanger.exchange(subgrid);
    }
    {
        HYPOS_PROFILE("stencil_boundary");
        residual += sweep(subgrid, 1);
    }

    subgrid.applyPhysicalBoundary();

    return std::sqrt(residual);
}

Index RedBlackGSSolver::solve(Subgrid& subgrid,
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

        Real globalResidual = 0.0;
        {
            HYPOS_PROFILE("residual_allreduce");
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
    }

    totalTimer.stop();
    HYPOS_INFO("RedBlackGS finished in " << completed << " iterations, time = "
                                         << totalTimer.elapsedSeconds() << " s");

    return completed;
}

} // namespace hypo
