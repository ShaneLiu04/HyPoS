#pragma once

#include "core/types.hpp"
#include "grid/subgrid.hpp"
#include "comm/halo_exchanger.hpp"
#include <mpi.h>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace hypo {

class HaloExchanger;

/**
 * @brief Abstract interface for Poisson equation solvers.
 */
class PoissonSolver {
public:
    virtual ~PoissonSolver() = default;

    /**
     * @brief Callback invoked after each completed iteration (1-based count).
     * Used by the driver for intermediate solution saves. Must not throw.
     */
    using ProgressCallback = std::function<void(Index iteration)>;

    void setProgressCallback(ProgressCallback callback) {
        progressCallback_ = std::move(callback);
    }

    /**
     * @brief Set how often (in iterations) the convergence criterion is
     * evaluated via a global reduction. Values < 1 and wrapped-around
     * negatives (size_t two's complement) are clamped to 1 (every
     * iteration); intervals beyond any sane iteration budget are treated
     * the same. Only stationary-iteration solvers (Jacobi, Red-Black GS,
     * and the mg2 outer cycle) honor this — CG's beta recurrence needs a
     * reduction every iteration and ignores the setting.
     */
    void setResidualCheckInterval(Index interval) noexcept {
        constexpr Index kMaxResidualCheckInterval = Index(1) << 30;
        residualCheckInterval_ =
            (interval < Index(1) || interval > kMaxResidualCheckInterval)
                ? Index(1)
                : interval;
    }

    /**
     * @brief Solve the Poisson equation on the given subdomain.
     * @param subgrid Local subdomain with current u, uNext, and rhs.
     * @param exchanger Halo exchange handler for boundary synchronization.
     * @param maxIter Maximum number of iterations.
     * @param tolerance Convergence criterion on L2 residual.
     * @return Number of iterations actually performed.
     */
    virtual Index solve(Subgrid& subgrid,
                        HaloExchanger& exchanger,
                        Index maxIter,
                        Real tolerance) = 0;

    /**
     * @brief Perform one iteration (for external convergence monitoring).
     * @return L2 residual after this iteration.
     */
    virtual Real iterate(Subgrid& subgrid, HaloExchanger& exchanger) = 0;

    virtual std::string name() const = 0;

    /**
     * @brief Global L2 residual after the most recent iteration.
     * Returns 0 before any iteration has been performed.
     */
    virtual Real lastResidual() const noexcept { return 0.0; }

protected:
    void notifyProgress(Index iteration) const {
        if (progressCallback_) {
            progressCallback_(iteration);
        }
    }

    /**
     * @brief Convergence-check interval in iterations (see
     * setResidualCheckInterval); 1 = check after every iteration.
     */
    Index residualCheckInterval_ = 1;

private:
    ProgressCallback progressCallback_;
};

/**
 * @brief Classic Jacobi solver with OpenMP acceleration.
 * Supports optional communication-computation overlap: the interior box is
 * updated while halo messages are in flight, followed by the boundary band.
 *
 * Residual semantics (AR004): the fused update accumulates sum(diff^2) with
 * diff = uNew - u, which satisfies the exact algebraic identity
 * r = D * diff for the system A u = b (A = D*I - S, b = -rhs). iterate()
 * therefore returns denom * ||diff|| — the true residual of the state
 * BEFORE the update. solve() evaluates the convergence criterion on this
 * converted value and, at loop exit (converged or maxIter), runs one
 * confirmation scan so lastResidual() is the true residual of the final
 * iterate.
 */
class JacobiSolver : public PoissonSolver {
public:
    explicit JacobiSolver(bool overlapComm = false);

    Index solve(Subgrid& subgrid,
                HaloExchanger& exchanger,
                Index maxIter,
                Real tolerance) override;

    Real iterate(Subgrid& subgrid, HaloExchanger& exchanger) override;

    std::string name() const override { return "jacobi"; }

    bool overlapEnabled() const noexcept { return overlapComm_; }

    Real lastResidual() const noexcept override { return lastResidual_; }

private:
    /**
     * @brief Update the given interior region, accumulating the local L2
     * residual contribution (sum of squared updates) into the return value.
     * The accumulation order across regions is fixed by iterate() so that
     * overlap on/off produce bit-identical residuals.
     * @param k0,k1 Active layer range: [0,1) for 2D, [kBegin,kEnd) for 3D.
     */
    Real updateRegion(Subgrid& subgrid,
                      Index i0, Index i1,
                      Index j0, Index j1,
                      Index k0, Index k1) const;

    bool overlapComm_ = false;
    Real lastResidual_ = 0.0;
};

/**
 * @brief Red-Black Gauss-Seidel solver.
 * Two colored half-sweeps per iteration; a halo exchange after each sweep
 * shares the freshly updated color. Global coloring uses the subgrid offsets,
 * so the same color assignment holds across ranks.
 *
 * Residual semantics (AR004): in-place updates leave no algebraic relation
 * between the update diff and the residual, so the true residual
 * ||A u - b|| is obtained from a dedicated scan (exchange + kernel +
 * reduction). iterate() always scans; solve() scans every
 * residualCheckInterval_ iterations and, at loop exit, runs one
 * confirmation scan so lastResidual() is the true residual of the final
 * iterate.
 */
class RedBlackGSSolver : public PoissonSolver {
public:
    Index solve(Subgrid& subgrid,
                HaloExchanger& exchanger,
                Index maxIter,
                Real tolerance) override;

    Real iterate(Subgrid& subgrid, HaloExchanger& exchanger) override;

    std::string name() const override { return "red_black_gs"; }

    Real lastResidual() const noexcept override { return lastResidual_; }

    /**
     * @brief Lightweight smoothing entry for multigrid drivers: `sweeps`
     * full RBGS iterations (two colored half-sweeps with halo exchanges
     * plus physical boundaries), no residual scan. Composes the exact
     * iterate() core, so smooth(n) and n x iterate() produce bit-identical
     * u (existing paths untouched).
     */
    void smooth(Subgrid& subgrid, HaloExchanger& exchanger, Index sweeps) const;

private:
    /**
     * @brief In-place sweep over one color; returns the local sum of squared
     * updates (residual contribution).
     * @param parity 0 for red cells, 1 for black cells (global parity).
     */
    Real sweep(Subgrid& subgrid, int parity) const;

    /**
     * @brief One full iteration without the residual scan: two colored
     * half-sweeps with halo exchanges, then applyPhysicalBoundary.
     * @return Local sum of squared updates (diagnostic only; not a
     *         convergence criterion — see class doc).
     */
    Real iterateCore(Subgrid& subgrid, HaloExchanger& exchanger) const;

    Real lastResidual_ = 0.0;
};

/**
 * @brief Conjugate gradient solver (unpreconditioned) for the SPD system
 * M u = -rhs, where M = 4I - S (2D) / 6I - S (3D) is the discrete negative
 * Laplacian. Each iteration performs one matvec (with a halo exchange).
 */
class CGSolver : public PoissonSolver {
public:
    Index solve(Subgrid& subgrid,
                HaloExchanger& exchanger,
                Index maxIter,
                Real tolerance) override;

    Real iterate(Subgrid& subgrid, HaloExchanger& exchanger) override;

    std::string name() const override { return "cg"; }

    Real lastResidual() const noexcept override { return lastResidual_; }

private:
    /**
     * @brief Apply the operator: ap = M * p (exchange + physical halos first).
     */
    void matvec(Subgrid& subgrid, HaloExchanger& exchanger, Real* p, Real* ap);

    /**
     * @brief Global sum of products over the interior: allreduce(sum(a*b)).
     */
    Real dotGlobal(const Subgrid& subgrid, const Real* a, const Real* b) const;

    void axpyInterior(Subgrid& subgrid, Real alpha, const Real* x, Real* y) const;

    /**
     * @brief Direction update on the interior: p = r + beta * p.
     * Element-wise, so the thread partition does not affect any value.
     */
    void updatePInterior(Subgrid& subgrid, Real beta);

    AlignedBuffer<Real> r_;
    AlignedBuffer<Real> p_;
    AlignedBuffer<Real> ap_;
    Real rho_ = 0.0;
    Real lastResidual_ = 0.0;
    bool stateReady_ = false;
};

/**
 * @brief Pipelined conjugate gradient solver (AR009, CG-1 single-step
 * pipeline a la Chronopoulos–Gear 1989) for the same SPD system M u = -rhs.
 *
 * Per iteration: ONE packed non-blocking MPI_Iallreduce carrying
 * [rho, r·Ar] (issue + immediate Wait — no true computation overlap, D1)
 * and ZERO blocking MPI_Allreduce (vs cg's two); q = A p is maintained by
 * recurrence (one matvec per iteration on r). alpha uses the nu-recurrence
 * nu_m = m_m − beta_{m−1}^2 * nu_{m−1} (p·Ap != r·Ar for m >= 1, so r·Ar
 * alone is NOT a valid denominator — design D1 identity). Loop entry
 * always has no pending request (issue-then-Wait), so every exit path is
 * free of dangling requests by construction. Convergence is judged on
 * sqrt(rho) BEFORE the breakdown check (on the convergence iteration the
 * nu-recurrence cancels catastrophically and would false-positive); at
 * loop exit the true residual is re-evaluated via globalTrueResidual (G3)
 * so lastResidual() is the reviewed value.
 */
class PipelinedCGSolver : public PoissonSolver {
public:
    Index solve(Subgrid& subgrid,
                HaloExchanger& exchanger,
                Index maxIter,
                Real tolerance) override;

    Real iterate(Subgrid& subgrid, HaloExchanger& exchanger) override;

    std::string name() const override { return "pcg"; }

    Real lastResidual() const noexcept override { return lastResidual_; }

private:
    AlignedBuffer<Real> r_;
    AlignedBuffer<Real> p_;
    AlignedBuffer<Real> q_;
    AlignedBuffer<Real> ar_;
    Real rho_ = 0.0;
    Real m_ = 0.0;
    Real nu_ = 0.0;
    Real lastResidual_ = 0.0;
    bool stateReady_ = false;
};

/**
 * @brief Two-level geometric multigrid correction scheme (AR007, GUIDE B1a):
 * RBGS pre-smoothing -> residual restriction (full weighting) -> replicated
 * coarse CG exact solve on MPI_COMM_SELF -> vertex-coincident prolongation
 * of the correction -> RBGS post-smoothing, repeated until the global true
 * residual meets the tolerance. One "iteration" = one outer cycle (the
 * solver's maxIter/return/progress/residual-check-interval are all in
 * cycle units).
 *
 * Scope guards (HYPOS_ERROR + zero iterations + u untouched): 2D only
 * (nz == 1), even global nx/ny, Dirichlet BC, coarse grid >= 4x4.
 */
class TwoLevelMGSolver : public PoissonSolver {
public:
    Index solve(Subgrid& subgrid,
                HaloExchanger& exchanger,
                Index maxIter,
                Real tolerance) override;

    Real iterate(Subgrid& subgrid, HaloExchanger& exchanger) override;

    std::string name() const override { return "mg2"; }

    Real lastResidual() const noexcept override { return lastResidual_; }

    static constexpr Index kPreSmoothSweeps = 2;
    static constexpr Index kPostSmoothSweeps = 2;

private:
    bool ensureInitialized(Subgrid& subgrid, HaloExchanger& exchanger);
    void cycleOnce(Subgrid& subgrid, HaloExchanger& exchanger);

    RedBlackGSSolver smoother_;
    CGSolver coarseSolver_;
    std::unique_ptr<Subgrid> coarse_;
    std::unique_ptr<PointToPointExchanger> coarseExchanger_;
    AlignedBuffer<Real> fineResidual_;
    AlignedBuffer<Real> coarsePack_;
    AlignedBuffer<Real> coarseAll_;
    std::vector<Index> rectX0_;
    std::vector<Index> rectY0_;
    std::vector<Index> rectNX_;
    std::vector<Index> rectNY_;
    std::vector<int> counts_;
    std::vector<int> displs_;
    Index myBx_ = 0;
    Index myBy_ = 0;
    Index myCnx_ = 0;
    Index myCny_ = 0;
    Index nxH_ = 0;
    Index nyH_ = 0;
    Index nxLocal_ = 0;
    Index nyLocal_ = 0;
    Index coarseIters_ = 0;
    Real coarseTolerance_ = 1e-10;
    Real lastResidual_ = 0.0;
    bool valid_ = false;
};

} // namespace hypo
