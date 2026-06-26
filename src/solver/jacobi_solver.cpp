#include "solver/solver.hpp"
#include "comm/halo_exchanger.hpp"
#include "perf/profiler.hpp"
#include "utils/logger.hpp"
#include "core/exception.hpp"
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hypo {

Real JacobiSolver::jacobiUpdate(Subgrid& subgrid) const {
    const Index nx = subgrid.nxLocal();
    const Index ny = subgrid.nyLocal();
    const Index nz = subgrid.nzLocal();
    const Int hw = subgrid.haloWidth();

    Real* u = subgrid.u().data();
    Real* uNext = subgrid.uNext().data();
    const Real* rhs = subgrid.rhs().data();

    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    const Index nzT = subgrid.nzTotal();

    Real residual = 0.0;

    if (nz == 1) {
        // 2D Jacobi
        const Real invDenom = 1.0 / 4.0;
        #pragma omp parallel for reduction(+:residual) schedule(static)
        for (Index j = 0; j < ny; ++j) {
            const Index gj = j + hw;
            #pragma omp simd reduction(+:residual)
            for (Index i = 0; i < nx; ++i) {
                const Index gi = i + hw;
                const Index idx = gj * nxT + gi;
                const Index idxL = idx - 1;
                const Index idxR = idx + 1;
                const Index idxD = idx - nxT;
                const Index idxU = idx + nxT;

                Real uNew = (u[idxL] + u[idxR] + u[idxD] + u[idxU] - rhs[idx]) * invDenom;
                uNext[idx] = uNew;
                Real diff = uNew - u[idx];
                residual += diff * diff;
            }
        }
    } else {
        // 3D Jacobi
        const Real invDenom = 1.0 / 6.0;
        #pragma omp parallel for reduction(+:residual) schedule(static)
        for (Index k = 0; k < nz; ++k) {
            const Index gk = k + hw;
            for (Index j = 0; j < ny; ++j) {
                const Index gj = j + hw;
                #pragma omp simd reduction(+:residual)
                for (Index i = 0; i < nx; ++i) {
                    const Index gi = i + hw;
                    const Index idx = (gk * nyT + gj) * nxT + gi;
                    const Index idxL = idx - 1;
                    const Index idxR = idx + 1;
                    const Index idxD = idx - nxT;
                    const Index idxU = idx + nxT;
                    const Index idxB = idx - nxT * nyT;
                    const Index idxF = idx + nxT * nyT;

                    Real uNew = (u[idxL] + u[idxR] + u[idxD] + u[idxU] + u[idxB] + u[idxF] - rhs[idx]) * invDenom;
                    uNext[idx] = uNew;
                    Real diff = uNew - u[idx];
                    residual += diff * diff;
                }
            }
        }
    }

    return std::sqrt(residual);
}

Real JacobiSolver::iterate(Subgrid& subgrid, HaloExchanger& exchanger) {
    // Exchange halos before update
    exchanger.exchange(subgrid);

    Real residual = jacobiUpdate(subgrid);

    // Swap u and uNext
    subgrid.swapU();

    return residual;
}

Index JacobiSolver::solve(Subgrid& subgrid,
                          HaloExchanger& exchanger,
                          Index maxIter,
                          Real tolerance) {
    Real residual = 0.0;
    Index iter = 0;
    Timer totalTimer;
    totalTimer.start();

    for (iter = 0; iter < maxIter; ++iter) {
        {
            HYPOS_PROFILE("jacobi_iteration");
            residual = iterate(subgrid, exchanger);
        }

        // Normalize residual by grid size
        Real globalResidual = 0.0;
        Real localResidual = residual * residual;
        MPI_Allreduce(&localResidual, &globalResidual, 1, MPI_DOUBLE, MPI_SUM, subgrid.comm());
        globalResidual = std::sqrt(globalResidual);

        if (globalResidual < tolerance) {
            if (iter % 100 == 0 || iter < 10) {
                HYPOS_INFO("Converged at iteration " << iter << ", residual = " << globalResidual);
            }
            break;
        }

        if (iter % 500 == 0) {
            HYPOS_INFO("Iteration " << iter << ", residual = " << globalResidual);
        }
    }

    totalTimer.stop();
    HYPOS_INFO("Solver finished in " << iter << " iterations, time = " << totalTimer.elapsedSeconds() << " s");

    return iter;
}

} // namespace hypo
