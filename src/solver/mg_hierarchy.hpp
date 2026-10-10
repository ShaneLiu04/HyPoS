#pragma once

#include "solver/solver.hpp"

#include <memory>
#include <vector>

namespace hypo {

/**
 * @brief Multigrid coarsening chain + recursive V-cycle engine (AR008,
 * GUIDE B1b).
 *
 * The hierarchy owns every level as an independent Subgrid (design D9):
 * level 0 is an owned "shadow" mirror of the driver layout (same local
 * size, offsets, neighbors, and communicator), levels >= 1 are replicated
 * full-grid COMM_SELF subgrids. Owning the shadow matters because every
 * reused kernel (RBGS smoothing, residualFieldLocal, restriction,
 * prolongation) reads and writes the Subgrid's own u()/rhs() fields — the
 * mgcg preconditioner step has no host fields on the driver subgrid (they
 * are occupied by the PCG state), so both drivers marshal data in and out
 * of the shadow instead.
 *
 * Scale compensation (design D3): restricting level l gives
 * rhs_{l+1}() = -4 * R(rho_l) — the factor 4 = (H_l/H_{l+1})^2 compensates
 * the rediscretization eigenvalue mismatch (PERFORMANCE §14), the minus
 * sign is the repo-wide kernel convention (kernels solve A x = -rhs()).
 *
 * Cycle shape (design D11, W-cycle): the compensation above is exact only
 * for smooth modes; stray mid-frequency content (prolongation/interpolation
 * residue) gets over-compensated and a plain V-cycle recursively amplifies
 * it (measured: [1,1] factor 2.01/cycle on 256^2, 2.69 on 512^2). Doing
 * TWO coarse corrections per level (the second one on the refreshed
 * residual) restores the uniform convergence the variational Galerkin
 * scheme would give (measured: 1e-4/cycle, level-count independent).
 *
 * Root solve (design D4): CG with the RELATIVE tolerance
 * 1e-12 * ||b_root||, which removes the absolute-tolerance starvation mg2
 * exhibits below ~1e-8 (AR007 leftover Minor).
 */
class MGHierarchy {
public:
    /**
     * @brief Build the coarsening chain for the driver layout mirrored by
     * `fine`. Idempotent for an unchanged local size.
     *
     * Scope guards (HYPOS_ERROR + false): 2D only, Dirichlet BC, even
     * global nx/ny, chain root >= 4 per dim (narrow grids like 8x20
     * bottom out at 2x5 and are rejected), and the level-0 rectangles
     * must tile the first coarse grid.
     */
    bool initialize(const Subgrid& fine);

    /**
     * @brief Recursive W-cycle on the hierarchy's own fields: on level l
     * approximately solve A_l u_l = b_l (kCycleDepth coarse corrections
     * per level, design D11). Level 0's u halos must be coherent on entry
     * (solveCycle/applyPreconditioner guarantee this).
     */
    void vcycle(Index level) noexcept;

    /**
     * @brief One cycle for the mgv driver: `b` is the driver rhs() field
     * content (rhs-space, copied verbatim into the shadow), `u` is the
     * current iterate marshalled in and out of the shadow.
     */
    void solveCycle(const Real* b, Real* u) noexcept;

    /**
     * @brief Preconditioner step for the mgcg driver: z ~= A^-1 r.
     * `r` is the PCG recurrence residual (f-space); the kernel convention
     * A x = -rhs() means the shadow stores -r (design N1 sign contract).
     * The shadow u() is zeroed here (fresh initial guess every call).
     */
    void applyPreconditioner(const Real* r, Real* z) noexcept;

    Index levels() const noexcept { return static_cast<Index>(dims_.size() / 2); }

    /** @brief Global dims of the chain root (last level). */
    void rootDims(Index& nx, Index& ny) const noexcept;

    /** @brief Flat per-level dims: (nx, ny) per level, level 0 first. */
    const std::vector<Index>& levelDims() const noexcept { return dims_; }

    /** @brief Accumulated coarse-root CG iterations since the last reset. */
    Index coarseCgIterations() const noexcept { return coarseIters_; }

    void resetCoarseCgIterations() noexcept { coarseIters_ = 0; }

    static constexpr Index kPreSmoothSweeps = 2;
    static constexpr Index kPostSmoothSweeps = 2;
    /** @brief Coarse-correction visits per level: 2 = W-cycle (D11). */
    static constexpr Index kCycleDepth = 2;

private:
    struct Level {
        std::unique_ptr<Subgrid> grid;
        std::unique_ptr<PointToPointExchanger> ex;
        Index nxH = 0; // global dims at this level
        Index nyH = 0;
    };

    std::vector<Level> levels_;
    std::vector<Index> dims_; // flat (nx, ny) per level

    // Per-level residual field. AlignedBuffer does not zero on allocate
    // (posix_memalign), so buffers are filled with 0 in initialize; the
    // per-cycle applyPhysicalBoundary(rho, 0) keeps the physical halo band
    // coherent (design M2 — restriction reads that band at boundaries).
    std::vector<AlignedBuffer<Real>> rho_;
    // Packed restriction buffer for the l >= 1 pure-kernel segments
    // (packs_[l] receives the restriction of level l; index 0 unused —
    // the level 0 -> 1 segment goes through the Allgatherv tables).
    std::vector<AlignedBuffer<Real>> packs_;
    AlignedBuffer<Real> coarsePack_; // level 0 -> 1 local rectangle
    AlignedBuffer<Real> coarseAll_;  // level 0 -> 1 assembled full grid

    // Replicated-assembly tables for the level 0 -> 1 restriction
    // (mg2's ensureInitialized pattern).
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

    RedBlackGSSolver smoother_;
    CGSolver coarseSolver_;
    Index coarseIters_ = 0;
    Index nxLocal_ = 0;
    Index nyLocal_ = 0;
    bool ready_ = false;
};

/**
 * @brief Multigrid V-cycle driver (AR008, GUIDE B1b, --solver mgv):
 * repeated VCycles on the driver problem until the global true residual
 * meets the tolerance. One "iteration" = one outer cycle (maxIter/return/
 * progress/residual-check-interval are all in cycle units, mirroring mg2).
 *
 * Scope guards (HYPOS_ERROR + zero iterations + u untouched): 2D only,
 * even global nx/ny, Dirichlet BC, chain root >= 4 per dim.
 */
class VCycleMGSolver : public PoissonSolver {
public:
    Index solve(Subgrid& subgrid,
                HaloExchanger& exchanger,
                Index maxIter,
                Real tolerance) override;

    Real iterate(Subgrid& subgrid, HaloExchanger& exchanger) override;

    std::string name() const override { return "mgv"; }

    Real lastResidual() const noexcept override { return lastResidual_; }

private:
    bool ensureInitialized(Subgrid& subgrid, HaloExchanger& exchanger);

    MGHierarchy hierarchy_;
    Index nxLocal_ = 0;
    Index nyLocal_ = 0;
    Real lastResidual_ = 0.0;
    bool valid_ = false;
};

/**
 * @brief MG-preconditioned CG (AR008, GUIDE B1b, --solver mgcg):
 * PCG with one V-cycle per iteration as the preconditioner
 * (z = V-cycle(r), Fletcher-Reeves beta), reusing the CGSolver kernels
 * via the FP4 free-function extraction. See mg_pcg.cpp.
 *
 * Scope guards: identical to VCycleMGSolver.
 */
class MGPreconditionedCGSolver : public PoissonSolver {
public:
    Index solve(Subgrid& subgrid,
                HaloExchanger& exchanger,
                Index maxIter,
                Real tolerance) override;

    Real iterate(Subgrid& subgrid, HaloExchanger& exchanger) override;

    std::string name() const override { return "mgcg"; }

    Real lastResidual() const noexcept override { return lastResidual_; }

private:
    bool ensureInitialized(Subgrid& subgrid, HaloExchanger& exchanger);

    MGHierarchy hierarchy_;
    Index nxLocal_ = 0;
    Index nyLocal_ = 0;
    Real lastResidual_ = 0.0;
    bool valid_ = false;

    // PCG state (iterate() continues from the solve() state).
    AlignedBuffer<Real> r_;
    AlignedBuffer<Real> z_;
    AlignedBuffer<Real> p_;
    AlignedBuffer<Real> ap_;
    Real rho_ = 0.0;
    bool stateReady_ = false;
};

} // namespace hypo
