#include "solver/mg_hierarchy.hpp"
#include "solver/mg_operators.hpp"
#include "solver/residual.hpp"
#include "comm/halo_exchanger.hpp"
#include "core/exception.hpp"
#include "perf/profiler.hpp"
#include "perf/timer.hpp"
#include "utils/logger.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <vector>

namespace hypo {

bool MGHierarchy::initialize(const Subgrid& fine) {
    if (ready_ && fine.nxLocal() == nxLocal_ && fine.nyLocal() == nyLocal_) {
        return true;
    }
    ready_ = false;

    // Scope guards: the whole chain is the mg2 contract (2D Dirichlet on
    // even global dims); both drivers and the standalone hierarchy tests
    // enter through this single gate.
    if (fine.nzLocal() != 1) {
        HYPOS_ERROR("mgv: 3D grids are not supported (B1b scope is 2D-only)");
        return false;
    }
    if (fine.boundaryCondition() != BoundaryCondition::Dirichlet) {
        HYPOS_ERROR("mgv: only Dirichlet BC is supported (the rediscretized "
                    "Neumann coarse operator is singular)");
        return false;
    }

    unsigned long long localEnd[3] = {
        static_cast<unsigned long long>(fine.offsetX() + fine.nxLocal()),
        static_cast<unsigned long long>(fine.offsetY() + fine.nyLocal()),
        static_cast<unsigned long long>(fine.offsetZ() + fine.nzLocal())};
    unsigned long long globalEnd[3] = {0, 0, 0};
    MPI_Allreduce(localEnd, globalEnd, 3, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                  fine.comm());
    const Index nxGlobal = static_cast<Index>(globalEnd[0]);
    const Index nyGlobal = static_cast<Index>(globalEnd[1]);
    if (globalEnd[2] > 1) {
        HYPOS_ERROR("mgv: 3D grids are not supported (B1b scope is 2D-only)");
        return false;
    }
    if (nxGlobal % 2 != 0 || nyGlobal % 2 != 0) {
        HYPOS_ERROR("mgv: even global nx/ny required for coarsening, got "
                    << nxGlobal << " x " << nyGlobal);
        return false;
    }

    // Design D2: coarsen while both dims are even and some dim exceeds 8.
    // An odd mid-chain dim (34 -> 17) truncates naturally — the 17^2 root
    // is CG-solvable and not an error.
    std::vector<Index> chain;
    Index nx = nxGlobal;
    Index ny = nyGlobal;
    chain.push_back(nx);
    chain.push_back(ny);
    while (nx % 2 == 0 && ny % 2 == 0 && (nx > 8 || ny > 8)) {
        nx /= 2;
        ny /= 2;
        chain.push_back(nx);
        chain.push_back(ny);
    }
    if (nx < 4 || ny < 4) {
        // Narrow grids (8x20 -> 4x10 -> 2x5) bottom out below the 4x4
        // minimum; the chain generator refuses even though every driver
        // pre-check passed (design D2 backstop, U7).
        HYPOS_ERROR("mgv: coarsening chain roots at " << nx << " x " << ny
                    << ", below the 4x4 minimum");
        return false;
    }

    levels_.clear();
    dims_ = std::move(chain);
    rho_.clear();
    packs_.clear();
    rectX0_.clear();
    rectY0_.clear();
    rectNX_.clear();
    rectNY_.clear();
    counts_.clear();
    displs_.clear();

    // Design D9: level 0 is an OWNED mirror of the driver layout — the
    // reused kernels bind to Subgrid::u()/rhs(), so the hierarchy must
    // host its own working fields (the mgcg preconditioner has none on
    // the driver subgrid). Same offsets/neighbors keep the RBGS global
    // coloring and halo exchanges identical to the driver's.
    {
        Level l0;
        l0.nxH = nxGlobal;
        l0.nyH = nyGlobal;
        l0.grid = std::make_unique<Subgrid>(fine.nxLocal(), fine.nyLocal(), 1,
                                            fine.haloWidth(), fine.comm());
        l0.grid->setOffsets(fine.offsetX(), fine.offsetY(), 0);
        l0.grid->setNeighbors(fine.neighborLeft(), fine.neighborRight(),
                              fine.neighborDown(), fine.neighborUp());
        l0.grid->setBoundaryCondition(fine.boundaryCondition());
        l0.ex = std::make_unique<PointToPointExchanger>();
        l0.ex->initialize(*l0.grid);
        levels_.push_back(std::move(l0));
    }

    // Levels >= 1: replicated full-grid COMM_SELF subgrids (design D6).
    const Index nLevels = levels();
    for (Index l = 1; l < nLevels; ++l) {
        Level lv;
        lv.nxH = dims_[2 * l];
        lv.nyH = dims_[2 * l + 1];
        lv.grid = std::make_unique<Subgrid>(lv.nxH, lv.nyH, 1, 1, MPI_COMM_SELF);
        lv.grid->setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
                              MPI_PROC_NULL);
        lv.ex = std::make_unique<PointToPointExchanger>();
        lv.ex->initialize(*lv.grid);
        levels_.push_back(std::move(lv));
    }

    rho_.resize(nLevels);
    packs_.resize(nLevels);
    for (Index l = 0; l < nLevels; ++l) {
        rho_[l].allocate(levels_[l].grid->totalCells());
        rho_[l].fill(0.0);
        if (l >= 1 && l + 1 < nLevels) {
            packs_[l].allocate(levels_[l + 1].nxH * levels_[l + 1].nyH);
        }
    }

    // Level 0 -> 1 replicated assembly tables (mg2's ensureInitialized
    // pattern): each rank restricts its own coarse rectangle, the pieces
    // are Allgatherv'ed into the full level-1 rhs.
    {
        const Index nx1 = dims_[2];
        const Index ny1 = dims_[3];
        int size = 0;
        MPI_Comm_size(fine.comm(), &size);
        myBx_ = coarseRangeBegin(fine.offsetX(), fine.nxLocal());
        myBy_ = coarseRangeBegin(fine.offsetY(), fine.nyLocal());
        myCnx_ = coarseRangeCount(fine.offsetX(), fine.nxLocal());
        myCny_ = coarseRangeCount(fine.offsetY(), fine.nyLocal());

        std::vector<unsigned long long> meta(4 * static_cast<size_t>(size), 0);
        const unsigned long long mine[4] = {static_cast<unsigned long long>(myBx_),
                                            static_cast<unsigned long long>(myBy_),
                                            static_cast<unsigned long long>(myCnx_),
                                            static_cast<unsigned long long>(myCny_)};
        MPI_Allgather(mine, 4, MPI_UNSIGNED_LONG_LONG, meta.data(), 4,
                      MPI_UNSIGNED_LONG_LONG, fine.comm());

        rectX0_.resize(size);
        rectY0_.resize(size);
        rectNX_.resize(size);
        rectNY_.resize(size);
        counts_.resize(size);
        displs_.resize(size);
        Index total = 0;
        for (int r = 0; r < size; ++r) {
            rectX0_[r] = static_cast<Index>(meta[4 * r + 0]);
            rectY0_[r] = static_cast<Index>(meta[4 * r + 1]);
            rectNX_[r] = static_cast<Index>(meta[4 * r + 2]);
            rectNY_[r] = static_cast<Index>(meta[4 * r + 3]);
            const Index cnt = rectNX_[r] * rectNY_[r];
            counts_[r] = static_cast<int>(cnt);
            displs_[r] = static_cast<int>(total);
            total += cnt;
        }
        // Same tiling defense as mg2's fifth guard: a manual layout whose
        // coarse rectangles do not tile the coarse grid would assemble
        // silently wrong data.
        if (total != nx1 * ny1) {
            HYPOS_ERROR("mgv: subdomain layout does not tile the coarse grid ("
                        << total << " vs " << nx1 * ny1 << " cells)");
            return false;
        }

        coarsePack_.allocate(std::max<Index>(1, myCnx_ * myCny_));
        coarseAll_.allocate(nx1 * ny1);
    }

    nxLocal_ = fine.nxLocal();
    nyLocal_ = fine.nyLocal();
    ready_ = true;
    return true;
}

void MGHierarchy::rootDims(Index& nx, Index& ny) const noexcept {
    nx = dims_[dims_.size() - 2];
    ny = dims_.back();
}

void MGHierarchy::vcycle(Index level) noexcept {
    const Index L = static_cast<Index>(levels_.size()) - 1;
    Subgrid& sg = *levels_[level].grid;

    if (level == L) {
        // Design D4: relative-tolerance root solve. The caller (parent
        // level's restriction segment) zeroed u, which is CG's required
        // zero initial guess. b_root = -rhs(); skipping the solve on a
        // zero rhs avoids the CG breakdown on an already-converged system.
        Real bSq = 0.0;
        const Real* crhs = sg.rhs().data();
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                const Real b = -crhs[sg.index(i, j)];
                bSq += b * b;
            }
        }
        if (bSq > 0.0) {
            HYPOS_PROFILE("mgv_coarse_solve");
            coarseIters_ += coarseSolver_.solve(
                sg, *levels_[level].ex, levels_[level].nxH * levels_[level].nyH,
                1e-12 * std::sqrt(bSq));
        }
        return;
    }

    smoother_.smooth(sg, *levels_[level].ex, kPreSmoothSweeps);

    // Design D11 (W-cycle): kCycleDepth coarse corrections per level. The
    // second pass re-restricts the residual of the already-corrected
    // iterate, which cleans up the stray-mode over-compensation a single
    // V-shaped pass leaves behind (see the class comment for the measured
    // divergence data).
    Subgrid& csg = *levels_[level + 1].grid;
    for (Index g = 0; g < kCycleDepth; ++g) {
        if (level == 0) {
            // AR007 deviation-4 contract + prolongation staleness: the
            // residual reads u's neighbor halos, which the pre-smoother's
            // black sweep (first pass) and every prolongation (later
            // passes) leave stale.
            HYPOS_PROFILE("halo_exchange");
            levels_[0].ex->exchange(sg);
        }

        residualFieldLocal(sg, rho_[level].data());
        if (level == 0) {
            HYPOS_PROFILE("halo_exchange");
            levels_[0].ex->exchange(sg, rho_[level].data());
            exchangeCorners(sg, rho_[level].data());
        }
        // Design M2: the residual's physical halo band must be zeroed on
        // EVERY level — AlignedBuffer does not zero on allocate and the
        // restriction kernel reads the physical band at domain boundaries.
        sg.applyPhysicalBoundary(rho_[level].data(), 0.0);

        {
            HYPOS_PROFILE("mgv_restrict");
            if (level == 0) {
                // Replicated assembly: local rectangle restriction +
                // Allgatherv.
                restrictResidual(sg, rho_[0].data(), myBx_, myBy_, myCnx_,
                                 myCny_, coarsePack_.data());
                MPI_Allgatherv(coarsePack_.data(),
                               static_cast<int>(myCnx_ * myCny_), MPI_DOUBLE,
                               coarseAll_.data(), counts_.data(),
                               displs_.data(), MPI_DOUBLE, sg.comm());
                csg.zeroInitialize();
                const Index cnxT = csg.nxTotal();
                const Int chw = csg.haloWidth();
                Real* crhs = csg.rhs().data();
                const int size = static_cast<int>(rectX0_.size());
                for (int r = 0; r < size; ++r) {
                    const Real* src = coarseAll_.data() + displs_[r];
                    for (Index j = 0; j < rectNY_[r]; ++j) {
                        for (Index i = 0; i < rectNX_[r]; ++i) {
                            crhs[(rectY0_[r] + j + chw) * cnxT +
                                 rectX0_[r] + i + chw] =
                                -4.0 * src[j * rectNX_[r] + i];
                        }
                    }
                }
            } else {
                // Pure-kernel full-grid restriction (every rank already
                // owns the whole coarse level; no communication).
                restrictResidual(sg, rho_[level].data(), 0, 0,
                                 levels_[level + 1].nxH,
                                 levels_[level + 1].nyH,
                                 packs_[level].data());
                csg.zeroInitialize();
                const Index cnxT = csg.nxTotal();
                const Int chw = csg.haloWidth();
                const Index cnx = levels_[level + 1].nxH;
                const Index cny = levels_[level + 1].nyH;
                Real* crhs = csg.rhs().data();
                const Real* src = packs_[level].data();
                for (Index j = 0; j < cny; ++j) {
                    for (Index i = 0; i < cnx; ++i) {
                        crhs[(j + chw) * cnxT + i + chw] =
                            -4.0 * src[j * cnx + i];
                    }
                }
            }
        }

        vcycle(level + 1);

        {
            HYPOS_PROFILE("mgv_prolongate");
            prolongateCorrection(sg, csg);
        }
    }

    smoother_.smooth(sg, *levels_[level].ex, kPostSmoothSweeps);
}

void MGHierarchy::solveCycle(const Real* b, Real* u) noexcept {
    Subgrid& s0 = *levels_[0].grid;
    const Index total = s0.totalCells();
    Real* srhs = s0.rhs().data();
    Real* su = s0.u().data();
    // Whole-field copies including halos: the driver's halo state is part
    // of the iterate (the marshalled-back shadow halos are coherent after
    // each cycle's smoothing + physical boundaries).
    std::memcpy(srhs, b, static_cast<size_t>(total) * sizeof(Real));
    std::memcpy(su, u, static_cast<size_t>(total) * sizeof(Real));
    vcycle(0);
    std::memcpy(u, su, static_cast<size_t>(total) * sizeof(Real));
}

void MGHierarchy::applyPreconditioner(const Real* r, Real* z) noexcept {
    Subgrid& s0 = *levels_[0].grid;
    const Index total = s0.totalCells();
    Real* srhs = s0.rhs().data();
    Real* su = s0.u().data();
    // Design N1 sign contract: kernels solve A x = -rhs(), while r is the
    // PCG recurrence residual (f-space), so the shadow rhs() stores -r.
    // Element-wise negation is exact in IEEE arithmetic, which is what
    // makes the z == -solveCycle(b) identity in the unit test bit-exact.
    for (Index idx = 0; idx < total; ++idx) {
        srhs[idx] = -r[idx];
        su[idx] = 0.0;
    }
    vcycle(0);
    std::memcpy(z, su, static_cast<size_t>(total) * sizeof(Real));
}

// ---------------------------------------------------------------------------
// VCycleMGSolver
// ---------------------------------------------------------------------------

bool VCycleMGSolver::ensureInitialized(Subgrid& subgrid,
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

Index VCycleMGSolver::solve(Subgrid& subgrid,
                            HaloExchanger& exchanger,
                            Index maxIter,
                            Real tolerance) {
    Timer totalTimer;
    totalTimer.start();

    if (!ensureInitialized(subgrid, exchanger)) {
        return 0;
    }
    hierarchy_.resetCoarseCgIterations();

    Index cycle = 0;
    while (cycle < maxIter) {
        HYPOS_PROFILE("mgv_cycle");
        hierarchy_.solveCycle(subgrid.rhs().data(), subgrid.u().data());
        ++cycle;

        if (cycle % residualCheckInterval_ == 0) {
            Real globalResidual = 0.0;
            {
                HYPOS_PROFILE("residual_allreduce");
                globalResidual = globalTrueResidual(subgrid, exchanger);
            }
            lastResidual_ = globalResidual;
            notifyProgress(cycle);

            if (globalResidual < tolerance) {
                HYPOS_INFO("Converged at cycle " << cycle << ", residual = "
                                                << globalResidual);
                break;
            }
            if (cycle % 500 == 0) {
                HYPOS_INFO("Cycle " << cycle << ", residual = "
                                    << globalResidual);
            }
        } else {
            notifyProgress(cycle);
        }
    }

    if (cycle > 0) {
        HYPOS_PROFILE("residual_confirm");
        lastResidual_ = globalTrueResidual(subgrid, exchanger);
    }

    totalTimer.stop();
    // Smoothing accounting: every non-root level runs kPre + kPost sweeps
    // per cycle (the root is an exact CG solve), so the sweep total spans
    // the whole chain, not just the fine level.
    const Index nonRootLevels = hierarchy_.levels() - 1;
    const Index sweepsPerCycle =
        nonRootLevels * (MGHierarchy::kPreSmoothSweeps + MGHierarchy::kPostSmoothSweeps);
    HYPOS_INFO("VCycleMG finished in " << cycle << " cycles ("
               << cycle * sweepsPerCycle << " smoothing sweeps, "
               << hierarchy_.coarseCgIterations()
               << " coarse CG iterations), time = "
               << totalTimer.elapsedSeconds() << " s");

    return cycle;
}

Real VCycleMGSolver::iterate(Subgrid& subgrid, HaloExchanger& exchanger) {
    if (!ensureInitialized(subgrid, exchanger)) {
        return 0.0;
    }
    {
        HYPOS_PROFILE("mgv_cycle");
        hierarchy_.solveCycle(subgrid.rhs().data(), subgrid.u().data());
    }
    HYPOS_PROFILE("residual_allreduce");
    lastResidual_ = globalTrueResidual(subgrid, exchanger);
    return lastResidual_;
}

} // namespace hypo
