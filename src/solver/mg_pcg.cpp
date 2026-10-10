#include "solver/mg_hierarchy.hpp"
#include "solver/cg_kernels.hpp"
#include "solver/residual.hpp"
#include "comm/halo_exchanger.hpp"
#include "core/exception.hpp"
#include "perf/profiler.hpp"
#include "perf/timer.hpp"
#include "utils/logger.hpp"

#include <cmath>

namespace hypo {

// PCG state (buffers exclude nothing: totalCells covers halos so the
// memcpy-based marshalling in applyPreconditioner sees whole fields).
Index MGPreconditionedCGSolver::solve(Subgrid& subgrid,
                                      HaloExchanger& exchanger,
                                      Index maxIter,
                                      Real tolerance) {
    Timer totalTimer;
    totalTimer.start();

    if (!ensureInitialized(subgrid, exchanger)) {
        return 0;
    }

    const Index total = subgrid.totalCells();
    if (r_.size() != total) {
        r_.allocate(total);
    }
    if (z_.size() != total) {
        z_.allocate(total);
    }
    if (p_.size() != total) {
        p_.allocate(total);
    }
    if (ap_.size() != total) {
        ap_.allocate(total);
    }
    #pragma omp parallel for schedule(static)
    for (Index idx = 0; idx < total; ++idx) {
        r_[idx] = 0.0;
        z_[idx] = 0.0;
        p_[idx] = 0.0;
        ap_[idx] = 0.0;
    }

    // System M u = b with M = 4I - S (2D only, hierarchy scope) and
    // b = -rhs. Initial guess x = 0 (set up by the driver), so r = b.
    // First PCG step: z0 = M_approx^-1 r0, rho0 = (z0, r0), p0 = z0.
    {
        const Index nxT = subgrid.nxTotal();
        const Real* rhs = subgrid.rhs().data();
        #pragma omp parallel for schedule(static)
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            #pragma omp simd
            for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                const Index idx = j * nxT + i;
                r_[idx] = -rhs[idx];
            }
        }
    }
    subgrid.applyPhysicalBoundary();

    hierarchy_.applyPreconditioner(r_.data(), z_.data());
    Real rho = cgDotGlobal(subgrid, z_.data(), r_.data());
    {
        const Index nxT = subgrid.nxTotal();
        const Real* z = z_.data();
        Real* pp = p_.data();
        #pragma omp parallel for schedule(static)
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            #pragma omp simd
            for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                const Index idx = j * nxT + i;
                pp[idx] = z[idx];
            }
        }
    }

    Real residual = std::sqrt(cgDotGlobal(subgrid, r_.data(), r_.data()));
    lastResidual_ = residual;
    stateReady_ = true;

    Index completed = 0;
    if (residual >= tolerance) {
        for (; completed < maxIter;) {
            HYPOS_PROFILE("mgcg_iteration");
            cgMatvec(subgrid, exchanger, p_.data(), ap_.data());
            const Real pap = cgDotGlobal(subgrid, p_.data(), ap_.data());
            if (!(pap > 0.0)) {
                HYPOS_WARN("mgcg breakdown: p^T M p <= 0 at iteration "
                           << completed);
                break;
            }
            const Real alpha = rho / pap;
            cgAxpyInterior(subgrid, alpha, p_.data(), subgrid.u().data());
            cgAxpyInterior(subgrid, -alpha, ap_.data(), r_.data());

            const Real rr = cgDotGlobal(subgrid, r_.data(), r_.data());
            residual = std::sqrt(rr);
            lastResidual_ = residual;
            ++completed;
            notifyProgress(completed);

            if (residual < tolerance) {
                HYPOS_INFO("Converged at iteration " << completed
                                                     << ", residual = "
                                                     << residual);
                break;
            }
            if (completed % 500 == 0) {
                HYPOS_INFO("Iteration " << completed << ", residual = "
                                        << residual);
            }

            hierarchy_.applyPreconditioner(r_.data(), z_.data());
            const Real rhoNew = cgDotGlobal(subgrid, z_.data(), r_.data());
            const Real beta = rhoNew / rho; // Fletcher-Reeves (design D5)
            rho = rhoNew;
            cgUpdatePInterior(subgrid, z_.data(), p_.data(), beta);
        }
    }

    subgrid.applyPhysicalBoundary();
    // G3 convention: lastResidual reports the TRUE residual — the
    // recursive ||r|| is exact in exact arithmetic but drifts by roundoff
    // once the iterate converges; confirm against b - M u.
    if (completed > 0) {
        HYPOS_PROFILE("residual_confirm");
        lastResidual_ = globalTrueResidual(subgrid, exchanger);
    }
    rho_ = rho; // iterate() continues from this state

    totalTimer.stop();
    HYPOS_INFO("MGPreconditionedCG finished in " << completed
                << " iterations (preconditioner W-cycles included), residual = "
                << lastResidual_ << ", time = "
                << totalTimer.elapsedSeconds() << " s");

    return completed;
}

Real MGPreconditionedCGSolver::iterate(Subgrid& subgrid,
                                       HaloExchanger& exchanger) {
    if (!stateReady_) {
        HYPOS_WARN("MGPreconditionedCGSolver::iterate called before solve(); "
                   "returning 0");
        return 0.0;
    }
    if (!(rho_ > 0.0)) {
        return lastResidual_;
    }

    cgMatvec(subgrid, exchanger, p_.data(), ap_.data());
    const Real pap = cgDotGlobal(subgrid, p_.data(), ap_.data());
    if (!(pap > 0.0)) {
        HYPOS_WARN("mgcg breakdown: p^T M p <= 0");
        return lastResidual_;
    }
    const Real alpha = rho_ / pap;
    cgAxpyInterior(subgrid, alpha, p_.data(), subgrid.u().data());
    cgAxpyInterior(subgrid, -alpha, ap_.data(), r_.data());

    hierarchy_.applyPreconditioner(r_.data(), z_.data());
    const Real rhoNew = cgDotGlobal(subgrid, z_.data(), r_.data());
    const Real beta = rhoNew / rho_;
    rho_ = rhoNew;
    cgUpdatePInterior(subgrid, z_.data(), p_.data(), beta);

    lastResidual_ = std::sqrt(cgDotGlobal(subgrid, r_.data(), r_.data()));
    return lastResidual_;
}

bool MGPreconditionedCGSolver::ensureInitialized(Subgrid& subgrid,
                                                 HaloExchanger& exchanger) {
    (void)exchanger;
    if (valid_ && subgrid.nxLocal() == nxLocal_ &&
        subgrid.nyLocal() == nyLocal_) {
        return true;
    }
    valid_ = false;
    if (!hierarchy_.initialize(subgrid)) {
        return false;
    }
    nxLocal_ = subgrid.nxLocal();
    nyLocal_ = subgrid.nyLocal();
    valid_ = true;
    return true;
}

} // namespace hypo
