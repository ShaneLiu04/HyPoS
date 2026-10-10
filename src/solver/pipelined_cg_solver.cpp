// AR009 (design §4.1/§4.2): the pipelined CG solver — CG-1 single-step
// pipeline (Chronopoulos–Gear 1989) with the nu-recurrence and ONE packed
// non-blocking Iallreduce per iteration. See the class doc in solver.hpp
// for the algorithm contract; the bodies here follow cg_solver.cpp's
// structure and conventions (same init layout, same log lines, iterate()
// mirrors one solve() iteration step without notifyProgress).

#include "solver/solver.hpp"
#include "solver/cg_kernels.hpp"
#include "solver/residual.hpp"
#include "comm/halo_exchanger.hpp"
#include "perf/profiler.hpp"
#include "utils/logger.hpp"
#include <cmath>

namespace hypo {

Index PipelinedCGSolver::solve(Subgrid& subgrid,
                               HaloExchanger& exchanger,
                               Index maxIter,
                               Real tolerance) {
    Timer totalTimer;
    totalTimer.start();

    const Index total = subgrid.totalCells();
    if (r_.size() != total) {
        r_.allocate(total);
    }
    if (p_.size() != total) {
        p_.allocate(total);
    }
    if (q_.size() != total) {
        q_.allocate(total);
    }
    if (ar_.size() != total) {
        ar_.allocate(total);
    }
    #pragma omp parallel for schedule(static)
    for (Index idx = 0; idx < total; ++idx) {
        r_[idx] = 0.0;
        p_[idx] = 0.0;
        q_[idx] = 0.0;
        ar_[idx] = 0.0;
    }

    // System M u = b with M = 4I - S (2D) / 6I - S (3D) and b = -rhs.
    // Initial guess x = 0 (set up by the driver), so r = b, p = r.
    {
        const Index nxT = subgrid.nxTotal();
        const Index nyT = subgrid.nyTotal();
        const Real* rhs = subgrid.rhs().data();
        if (subgrid.nzLocal() == 1) {
            #pragma omp parallel for schedule(static)
            for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                #pragma omp simd
                for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                    const Index idx = j * nxT + i;
                    r_[idx] = -rhs[idx];
                    p_[idx] = r_[idx];
                }
            }
        } else {
            #pragma omp parallel for schedule(static)
            for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
                for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                    #pragma omp simd
                    for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                        const Index idx = (k * nyT + j) * nxT + i;
                        r_[idx] = -rhs[idx];
                        p_[idx] = r_[idx];
                    }
                }
            }
        }
    }
    subgrid.applyPhysicalBoundary();

    // Ar0 = A r0 (the one matvec of the init phase); q0 = Ar0. q's halo
    // cells are rewritten by the exchange inside every later cgMatvec, so
    // copying the whole buffer is safe.
    cgMatvec(subgrid, exchanger, r_.data(), ar_.data());
    #pragma omp parallel for schedule(static)
    for (Index idx = 0; idx < total; ++idx) {
        q_[idx] = ar_[idx];
    }

    // Init fused reduction: [rho0, m0] = [r0·r0, r0·Ar0], issued and
    // completed immediately (loop entry invariant: no pending request).
    Real local[2] = {cgDotLocal(subgrid, r_.data(), r_.data()),
                     cgDotLocal(subgrid, r_.data(), ar_.data())};
    Real global[2] = {0.0, 0.0};
    MPI_Request request;
    cgIallreduce2Start(subgrid, local, global, &request);
    cgIallreduce2Wait(&request);
    rho_ = global[0];
    m_ = global[1];
    nu_ = m_;  // nu0 = p0·q0 = r0·Ar0 = m0
    lastResidual_ = std::sqrt(rho_);
    stateReady_ = true;

    Index completed = 0;
    if (lastResidual_ >= tolerance) {
        while (completed < maxIter) {
            HYPOS_PROFILE("pcg_iteration");
            // ① Convergence BEFORE ② breakdown: on the convergence
            // iteration the nu-recurrence m − beta^2*nu cancels
            // catastrophically and would false-positive (design §4.1).
            if (std::sqrt(rho_) < tolerance) {
                HYPOS_INFO("Converged at iteration " << completed
                           << ", residual = " << std::sqrt(rho_));
                break;
            }
            if (!(nu_ > 0.0)) {
                HYPOS_WARN("PipelinedCG breakdown: nu = p^T M p <= 0 at iteration "
                           << completed);
                break;
            }
            // ③ alpha via the nu-recurrence value (p·Ap != r·Ar for m >= 1).
            const Real alpha = rho_ / nu_;
            // ④ x/r updates (x is hosted in subgrid.u()).
            cgAxpyInterior(subgrid, alpha, p_.data(), subgrid.u().data());
            cgAxpyInterior(subgrid, -alpha, q_.data(), r_.data());
            // ⑤ The single matvec of the iteration.
            cgMatvec(subgrid, exchanger, r_.data(), ar_.data());
            // ⑥ Fused issue+Wait for [rho+, m+] — no true overlap (D1);
            // the value is a halved sync count, not hidden latency.
            local[0] = cgDotLocal(subgrid, r_.data(), r_.data());
            local[1] = cgDotLocal(subgrid, r_.data(), ar_.data());
            cgIallreduce2Start(subgrid, local, global, &request);
            cgIallreduce2Wait(&request);
            // ⑦ Direction updates + nu/rho recurrence state.
            const Real rhoNew = global[0];
            const Real mNew = global[1];
            const Real beta = rhoNew / rho_;
            cgUpdatePInterior(subgrid, r_.data(), p_.data(), beta);   // p = r + beta*p
            cgUpdatePInterior(subgrid, ar_.data(), q_.data(), beta);  // q = Ar + beta*q
            nu_ = mNew - beta * beta * nu_;
            rho_ = rhoNew;
            m_ = mNew;
            // ⑧ Fresh residual + progress (no one-step lag — gate fix 1).
            lastResidual_ = std::sqrt(rhoNew);
            ++completed;
            notifyProgress(completed);
            if (completed % 500 == 0) {
                HYPOS_INFO("Iteration " << completed
                           << ", residual = " << lastResidual_);
            }
        }
    }

    subgrid.applyPhysicalBoundary();
    // G3: exit review against the true residual so lastResidual() is the
    // reviewed value, not the raw recurrence (mgcg convention).
    {
        HYPOS_PROFILE("residual_confirm");
        lastResidual_ = globalTrueResidual(subgrid, exchanger);
    }
    totalTimer.stop();
    HYPOS_INFO("PipelinedCG finished in " << completed << " iterations, time = "
                                       << totalTimer.elapsedSeconds() << " s");

    return completed;
}

Real PipelinedCGSolver::iterate(Subgrid& subgrid, HaloExchanger& exchanger) {
    if (!stateReady_) {
        HYPOS_WARN("PipelinedCGSolver::iterate called before solve(); returning 0");
        return 0.0;
    }
    if (!(rho_ > 0.0)) {
        return lastResidual_;
    }
    if (!(nu_ > 0.0)) {
        HYPOS_WARN("PipelinedCG breakdown: nu = p^T M p <= 0");
        return lastResidual_;
    }

    // One full iteration body (steps ③–⑧ of design §4.1): entry has no
    // pending request and leaves none (issue+Wait inside).
    const Real alpha = rho_ / nu_;
    cgAxpyInterior(subgrid, alpha, p_.data(), subgrid.u().data());
    cgAxpyInterior(subgrid, -alpha, q_.data(), r_.data());
    cgMatvec(subgrid, exchanger, r_.data(), ar_.data());
    Real local[2] = {cgDotLocal(subgrid, r_.data(), r_.data()),
                     cgDotLocal(subgrid, r_.data(), ar_.data())};
    Real global[2] = {0.0, 0.0};
    MPI_Request request;
    cgIallreduce2Start(subgrid, local, global, &request);
    cgIallreduce2Wait(&request);
    const Real rhoNew = global[0];
    const Real mNew = global[1];
    const Real beta = rhoNew / rho_;
    cgUpdatePInterior(subgrid, r_.data(), p_.data(), beta);
    cgUpdatePInterior(subgrid, ar_.data(), q_.data(), beta);
    nu_ = mNew - beta * beta * nu_;
    rho_ = rhoNew;
    m_ = mNew;
    lastResidual_ = std::sqrt(rhoNew);
    return lastResidual_;
}

} // namespace hypo
