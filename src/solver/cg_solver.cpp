#include "solver/solver.hpp"
#include "solver/cg_kernels.hpp"
#include "comm/halo_exchanger.hpp"
#include "perf/profiler.hpp"
#include "utils/logger.hpp"
#include <cmath>

namespace hypo {

// FP4 (AR008 design D5): kernel bodies moved verbatim from the CGSolver
// members below; the members are now thin forwarders so CGSolver call
// sites stay bit-identical (U6 golden guard).

Real cgDotLocal(const Subgrid& subgrid, const Real* a, const Real* b) {
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    Real local = 0.0;

    if (subgrid.nzLocal() == 1) {
        #pragma omp parallel for reduction(+:local) schedule(static)
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            #pragma omp simd reduction(+:local)
            for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                const Index idx = j * nxT + i;
                local += a[idx] * b[idx];
            }
        }
    } else {
        #pragma omp parallel for reduction(+:local) schedule(static)
        for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
            for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                #pragma omp simd reduction(+:local)
                for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                    const Index idx = (k * nyT + j) * nxT + i;
                    local += a[idx] * b[idx];
                }
            }
        }
    }
    return local;
}

Real cgBlockingAllreduce(const Subgrid& subgrid, Real local) {
    HYPOS_PROFILE("cg_blocking_allreduce");
    Real global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, subgrid.comm());
    return global;
}

Real cgDotGlobal(const Subgrid& subgrid, const Real* a, const Real* b) {
    // FP2 (AR009 design D2): split into cgDotLocal + cgBlockingAllreduce;
    // the composition is bit-identical to the pre-split body (U6 golden).
    return cgBlockingAllreduce(subgrid, cgDotLocal(subgrid, a, b));
}

void cgIallreduce2Start(const Subgrid& subgrid, const Real local[2], Real global[2], MPI_Request* req) {
    // FP3 (AR009 design §4.3): one packed count=2 Iallreduce behind the
    // "pcg_iallreduce" region. The recvbuf is a caller-provided buffer —
    // MPI forbids aliasing send/recv, and an MPI_Request cannot carry the
    // pointer across to the Wait side, so the design signature gains an
    // explicit `global` out-param (implementation-level deviation, logged
    // in tasks.md T002; the algorithm/timing/region contract is unchanged).
    HYPOS_PROFILE("pcg_iallreduce");
    MPI_Iallreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, subgrid.comm(), req);
}

void cgIallreduce2Wait(MPI_Request* req) {
    // The global sums are already in the caller's buffer once the request
    // completes (issue-then-immediately-Wait; no true overlap — D1).
    MPI_Wait(req, MPI_STATUS_IGNORE);
}

void cgMatvec(Subgrid& subgrid, HaloExchanger& exchanger, Real* p, Real* ap) {
    {
        HYPOS_PROFILE("halo_exchange");
        exchanger.exchange(subgrid, p);
        subgrid.applyPhysicalBoundary(p);
    }

    {
        HYPOS_PROFILE("stencil_interior");
        const Index nxT = subgrid.nxTotal();
        const Index nyT = subgrid.nyTotal();

        if (subgrid.nzLocal() == 1) {
            #pragma omp parallel for schedule(static)
            for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                #pragma omp simd
                for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                    const Index idx = j * nxT + i;
                    ap[idx] = 4.0 * p[idx] -
                              (p[idx - 1] + p[idx + 1] + p[idx - nxT] + p[idx + nxT]);
                }
            }
        } else {
            #pragma omp parallel for schedule(static)
            for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
                for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                    #pragma omp simd
                    for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                        const Index idx = (k * nyT + j) * nxT + i;
                        ap[idx] = 6.0 * p[idx] -
                                  (p[idx - 1] + p[idx + 1] +
                                   p[idx - nxT] + p[idx + nxT] +
                                   p[idx - nxT * nyT] + p[idx + nxT * nyT]);
                    }
                }
            }
        }
    }
}

void cgAxpyInterior(const Subgrid& subgrid, Real alpha, const Real* x, Real* y) {
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();

    if (subgrid.nzLocal() == 1) {
        #pragma omp parallel for schedule(static)
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            #pragma omp simd
            for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                const Index idx = j * nxT + i;
                y[idx] += alpha * x[idx];
            }
        }
    } else {
        #pragma omp parallel for schedule(static)
        for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
            for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                #pragma omp simd
                for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                    const Index idx = (k * nyT + j) * nxT + i;
                    y[idx] += alpha * x[idx];
                }
            }
        }
    }
}

void cgUpdatePInterior(const Subgrid& subgrid, const Real* v, Real* p, Real beta) {
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    const Real* rp = v;
    Real* pp = p;

    if (subgrid.nzLocal() == 1) {
        #pragma omp parallel for schedule(static)
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            #pragma omp simd
            for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                const Index idx = j * nxT + i;
                pp[idx] = rp[idx] + beta * pp[idx];
            }
        }
    } else {
        #pragma omp parallel for schedule(static)
        for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
            for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                #pragma omp simd
                for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                    const Index idx = (k * nyT + j) * nxT + i;
                    pp[idx] = rp[idx] + beta * pp[idx];
                }
            }
        }
    }
}

Real CGSolver::dotGlobal(const Subgrid& subgrid, const Real* a, const Real* b) const {
    return cgDotGlobal(subgrid, a, b);
}

void CGSolver::matvec(Subgrid& subgrid, HaloExchanger& exchanger, Real* p, Real* ap) {
    cgMatvec(subgrid, exchanger, p, ap);
}

void CGSolver::axpyInterior(Subgrid& subgrid, Real alpha, const Real* x, Real* y) const {
    cgAxpyInterior(subgrid, alpha, x, y);
}

void CGSolver::updatePInterior(Subgrid& subgrid, Real beta) {
    cgUpdatePInterior(subgrid, r_.data(), p_.data(), beta);
}

Index CGSolver::solve(Subgrid& subgrid,
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
    if (ap_.size() != total) {
        ap_.allocate(total);
    }
    #pragma omp parallel for schedule(static)
    for (Index idx = 0; idx < total; ++idx) {
        r_[idx] = 0.0;
        p_[idx] = 0.0;
        ap_[idx] = 0.0;
    }

    // System M u = b with M = 4I - S (2D) / 6I - S (3D) and b = -rhs.
    // Initial guess x = 0 (set up by the driver), so r = b.
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

    rho_ = dotGlobal(subgrid, r_.data(), r_.data());
    Real residual = std::sqrt(rho_);
    lastResidual_ = residual;
    stateReady_ = true;

    Index completed = 0;
    if (residual >= tolerance) {
        for (; completed < maxIter; ) {
            HYPOS_PROFILE("cg_iteration");
            matvec(subgrid, exchanger, p_.data(), ap_.data());
            const Real pap = dotGlobal(subgrid, p_.data(), ap_.data());
            if (!(pap > 0.0)) {
                HYPOS_WARN("CG breakdown: p^T M p <= 0 at iteration " << completed);
                break;
            }
            const Real alpha = rho_ / pap;
            axpyInterior(subgrid, alpha, p_.data(), subgrid.u().data());
            axpyInterior(subgrid, -alpha, ap_.data(), r_.data());

            const Real rhoNew = dotGlobal(subgrid, r_.data(), r_.data());
            residual = std::sqrt(rhoNew);
            lastResidual_ = residual;
            ++completed;
            notifyProgress(completed);

            if (residual < tolerance) {
                HYPOS_INFO("Converged at iteration " << completed << ", residual = " << residual);
                break;
            }
            if (completed % 500 == 0) {
                HYPOS_INFO("Iteration " << completed << ", residual = " << residual);
            }

            const Real beta = rhoNew / rho_;
            rho_ = rhoNew;
            updatePInterior(subgrid, beta);
        }
    }

    subgrid.applyPhysicalBoundary();
    totalTimer.stop();
    HYPOS_INFO("CG finished in " << completed << " iterations, time = "
                                 << totalTimer.elapsedSeconds() << " s");

    return completed;
}

Real CGSolver::iterate(Subgrid& subgrid, HaloExchanger& exchanger) {
    if (!stateReady_) {
        HYPOS_WARN("CGSolver::iterate called before solve(); returning 0");
        return 0.0;
    }
    if (!(rho_ > 0.0)) {
        return lastResidual_;
    }

    matvec(subgrid, exchanger, p_.data(), ap_.data());
    const Real pap = dotGlobal(subgrid, p_.data(), ap_.data());
    if (!(pap > 0.0)) {
        HYPOS_WARN("CG breakdown: p^T M p <= 0");
        return lastResidual_;
    }
    const Real alpha = rho_ / pap;
    axpyInterior(subgrid, alpha, p_.data(), subgrid.u().data());
    axpyInterior(subgrid, -alpha, ap_.data(), r_.data());

    const Real rhoNew = dotGlobal(subgrid, r_.data(), r_.data());
    lastResidual_ = std::sqrt(rhoNew);
    const Real beta = rhoNew / rho_;
    rho_ = rhoNew;
    updatePInterior(subgrid, beta);
    return lastResidual_;
}

} // namespace hypo
