#include "solver/solver.hpp"
#include "comm/halo_exchanger.hpp"
#include "perf/profiler.hpp"
#include "utils/logger.hpp"
#include "core/exception.hpp"
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hypo {

JacobiSolver::JacobiSolver(bool overlapComm)
    : overlapComm_(overlapComm) {
}

void JacobiSolver::updateRegion(Subgrid& subgrid,
                                Index i0, Index i1,
                                Index j0, Index j1,
                                Index k0, Index k1) const {
    if (i0 >= i1 || j0 >= j1 || k0 >= k1) {
        return;
    }

    Real* u = subgrid.u().data();
    Real* uNext = subgrid.uNext().data();
    const Real* rhs = subgrid.rhs().data();
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();

    if (subgrid.nzLocal() == 1) {
        // 2D: field lives in the k=0 plane
        const Real invDenom = 1.0 / 4.0;
        #pragma omp parallel for schedule(static)
        for (Index j = j0; j < j1; ++j) {
            #pragma omp simd
            for (Index i = i0; i < i1; ++i) {
                const Index idx = j * nxT + i;
                uNext[idx] = (u[idx - 1] + u[idx + 1] +
                              u[idx - nxT] + u[idx + nxT] -
                              rhs[idx]) * invDenom;
            }
        }
    } else {
        const Real invDenom = 1.0 / 6.0;
        #pragma omp parallel for schedule(static)
        for (Index k = k0; k < k1; ++k) {
            for (Index j = j0; j < j1; ++j) {
                #pragma omp simd
                for (Index i = i0; i < i1; ++i) {
                    const Index idx = (k * nyT + j) * nxT + i;
                    uNext[idx] = (u[idx - 1] + u[idx + 1] +
                                  u[idx - nxT] + u[idx + nxT] +
                                  u[idx - nxT * nyT] + u[idx + nxT * nyT] -
                                  rhs[idx]) * invDenom;
                }
            }
        }
    }
}

Real JacobiSolver::residualPass(Subgrid& subgrid) const {
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    const Real* u = subgrid.u().data();
    const Real* uNext = subgrid.uNext().data();

    Real residual = 0.0;

    if (subgrid.nzLocal() == 1) {
        #pragma omp parallel for reduction(+:residual) schedule(static)
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            #pragma omp simd reduction(+:residual)
            for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                const Index idx = j * nxT + i;
                const Real diff = uNext[idx] - u[idx];
                residual += diff * diff;
            }
        }
    } else {
        #pragma omp parallel for reduction(+:residual) schedule(static)
        for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
            for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                #pragma omp simd reduction(+:residual)
                for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                    const Index idx = (k * nyT + j) * nxT + i;
                    const Real diff = uNext[idx] - u[idx];
                    residual += diff * diff;
                }
            }
        }
    }

    return std::sqrt(residual);
}

Real JacobiSolver::iterate(Subgrid& subgrid, HaloExchanger& exchanger) {
    const Int hw = subgrid.haloWidth();

    // Inner box (independent of halo data) and boundary band split.
    const Index i0 = subgrid.iBegin() + hw;
    const Index i1 = subgrid.iEnd() - hw;
    const Index j0 = subgrid.jBegin() + hw;
    const Index j1 = subgrid.jEnd() - hw;

    // Active layer range: [0,1) for the 2D k=0 plane, [kBegin,kEnd) for 3D.
    const bool is2D = (subgrid.nzLocal() == 1);
    const Index kA0 = is2D ? 0 : subgrid.kBegin();
    const Index kA1 = is2D ? 1 : subgrid.kEnd();
    const Index kI0 = is2D ? 0 : (subgrid.kBegin() + hw);
    const Index kI1 = is2D ? 1 : (subgrid.kEnd() - hw);

    if (overlapComm_) {
        {
            HYPOS_PROFILE("halo_exchange");
            exchanger.beginExchange(subgrid);
        }
        {
            HYPOS_PROFILE("stencil_interior");
            updateRegion(subgrid, i0, i1, j0, j1, kI0, kI1);
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
            updateRegion(subgrid, i0, i1, j0, j1, kI0, kI1);
        }
    }

    // Boundary band: interior minus inner box, as six disjoint slabs.
    {
        HYPOS_PROFILE("stencil_boundary");
        updateRegion(subgrid, subgrid.iBegin(), i0, subgrid.jBegin(), subgrid.jEnd(), kA0, kA1);
        updateRegion(subgrid, i1, subgrid.iEnd(), subgrid.jBegin(), subgrid.jEnd(), kA0, kA1);
        updateRegion(subgrid, i0, i1, subgrid.jBegin(), j0, kA0, kA1);
        updateRegion(subgrid, i0, i1, j1, subgrid.jEnd(), kA0, kA1);
        updateRegion(subgrid, i0, i1, j0, j1, kA0, kI0);
        updateRegion(subgrid, i0, i1, j0, j1, kI1, kA1);
    }

    Real residual = 0.0;
    {
        HYPOS_PROFILE("residual_allreduce");
        residual = residualPass(subgrid);
    }

    subgrid.swapU();

    return residual;
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

        Real globalResidual = 0.0;
        {
            HYPOS_PROFILE("residual_allreduce");
            Real localSquared = localResidual * localResidual;
            MPI_Allreduce(&localSquared, &globalResidual, 1, MPI_DOUBLE, MPI_SUM, subgrid.comm());
            globalResidual = std::sqrt(globalResidual);
        }
        lastResidual_ = globalResidual;

        if (globalResidual < tolerance) {
            HYPOS_INFO("Converged at iteration " << completed << ", residual = " << globalResidual);
            break;
        }

        if (completed % 500 == 0) {
            HYPOS_INFO("Iteration " << completed << ", residual = " << globalResidual);
        }
    }

    totalTimer.stop();
    HYPOS_INFO("Solver finished in " << completed << " iterations, time = " << totalTimer.elapsedSeconds() << " s");

    return completed;
}

} // namespace hypo
