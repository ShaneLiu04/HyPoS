#include "solver/mg_operators.hpp"
#include "solver/solver.hpp"
#include "solver/residual.hpp"
#include "comm/halo_exchanger.hpp"
#include "core/exception.hpp"
#include "perf/profiler.hpp"
#include "perf/timer.hpp"
#include "utils/logger.hpp"

#include <algorithm>
#include <memory>
#include <vector>

namespace hypo {

Index coarseRangeBegin(Index offset, Index nLocal) noexcept {
    (void)nLocal;
    return (offset + 1) / 2;
}

Index coarseRangeCount(Index offset, Index nLocal) noexcept {
    const Index begin = coarseRangeBegin(offset, nLocal);
    const Index end = coarseRangeBegin(offset + nLocal, nLocal);
    return end - begin;
}

void restrictResidual(const Subgrid& fine, const Real* fineResidual,
                      Index coarseX0, Index coarseY0, Index cnx, Index cny,
                      Real* coarsePacked) noexcept {
    HYPOS_ASSERT(fine.nzLocal() == 1);
    const Index nxT = fine.nxTotal();
    const Index iB = fine.iBegin();
    const Index jB = fine.jBegin();
    const long long ox = static_cast<long long>(fine.offsetX());
    const long long oy = static_cast<long long>(fine.offsetY());
    const long long hw = static_cast<long long>(fine.haloWidth());

    #pragma omp parallel for schedule(static)
    for (Index j = 0; j < cny; ++j) {
        const long long cj = 2LL * static_cast<long long>(coarseY0 + j);
        for (Index i = 0; i < cnx; ++i) {
            const long long ci = 2LL * static_cast<long long>(coarseX0 + i);
            Real acc = 0.0;
            for (int dj = -1; dj <= 1; ++dj) {
                const long long fj = cj + dj - oy + hw;
                for (int di = -1; di <= 1; ++di) {
                    const long long fi = ci + di - ox + hw;
                    const Real w = (di == 0 ? 2.0 : 1.0) * (dj == 0 ? 2.0 : 1.0);
                    acc += w * fineResidual[fj * static_cast<long long>(nxT) + fi];
                }
            }
            coarsePacked[j * cnx + i] = acc / 16.0;
        }
    }
}

void prolongateCorrection(Subgrid& fine, const Subgrid& coarse) noexcept {
    HYPOS_ASSERT(fine.nzLocal() == 1);
    HYPOS_ASSERT(coarse.nzLocal() == 1);
    const Index nxT = fine.nxTotal();
    const Index cnxT = coarse.nxTotal();
    const long long hw = static_cast<long long>(fine.haloWidth());
    const long long chw = static_cast<long long>(coarse.haloWidth());
    const long long nxH = static_cast<long long>(coarse.nxLocal());
    const long long nyH = static_cast<long long>(coarse.nyLocal());
    const long long ox = static_cast<long long>(fine.offsetX());
    const long long oy = static_cast<long long>(fine.offsetY());
    const Real* ec = coarse.u().data();
    Real* u = fine.u().data();

    auto coarseAt = [&](long long I, long long J) -> Real {
        if (I < 0 || J < 0 || I >= nxH || J >= nyH) {
            return 0.0;
        }
        return ec[(J + chw) * static_cast<long long>(cnxT) + (I + chw)];
    };

    #pragma omp parallel for schedule(static)
    for (Index j = fine.jBegin(); j < fine.jEnd(); ++j) {
        const long long gj = oy + static_cast<long long>(j) - hw;
        const long long J0 = gj / 2;
        const Real wy0 = (gj % 2 == 0) ? 1.0 : 0.5;
        const Real wy1 = (gj % 2 == 0) ? 0.0 : 0.5;
        const long long J1 = (gj % 2 == 0) ? J0 : J0 + 1;
        for (Index i = fine.iBegin(); i < fine.iEnd(); ++i) {
            const long long gi = ox + static_cast<long long>(i) - hw;
            const long long I0 = gi / 2;
            const Real wx0 = (gi % 2 == 0) ? 1.0 : 0.5;
            const Real wx1 = (gi % 2 == 0) ? 0.0 : 0.5;
            const long long I1 = (gi % 2 == 0) ? I0 : I0 + 1;
            const Real v = wx0 * wy0 * coarseAt(I0, J0) +
                           wx0 * wy1 * coarseAt(I0, J1) +
                           wx1 * wy0 * coarseAt(I1, J0) +
                           wx1 * wy1 * coarseAt(I1, J1);
            u[static_cast<long long>(j) * static_cast<long long>(nxT) +
              static_cast<long long>(i)] += v;
        }
    }
}

void exchangeCorners(Subgrid& subgrid, Real* data) noexcept {
    if (subgrid.nzLocal() != 1) {
        return;
    }
    MPI_Comm comm = subgrid.comm();
    if (comm == MPI_COMM_NULL) {
        return;
    }
    int topology = MPI_UNDEFINED;
    if (MPI_Topo_test(comm, &topology) != MPI_SUCCESS || topology != MPI_CART) {
        return;
    }
    int ndims = 0;
    if (MPI_Cartdim_get(comm, &ndims) != MPI_SUCCESS || ndims < 2) {
        return;
    }
    if (ndims > 3) {
        return;
    }
    int dims[3] = {0, 0, 0};
    int periods[3] = {0, 0, 0};
    int coords[3] = {0, 0, 0};
    if (MPI_Cart_get(comm, ndims, dims, periods, coords) != MPI_SUCCESS) {
        return;
    }

    const bool hasX[2] = {subgrid.neighborLeft() != MPI_PROC_NULL,
                          subgrid.neighborRight() != MPI_PROC_NULL};
    const bool hasY[2] = {subgrid.neighborDown() != MPI_PROC_NULL,
                          subgrid.neighborUp() != MPI_PROC_NULL};
    if (!hasX[0] && !hasX[1] && !hasY[0] && !hasY[1]) {
        return;
    }

    constexpr int kCornerTag = 42;
    const Index nxT = subgrid.nxTotal();
    const Index iB = subgrid.iBegin();
    const Index iE = subgrid.iEnd();
    const Index jB = subgrid.jBegin();
    const Index jE = subgrid.jEnd();

    const int kDir[4][2] = {{-1, -1}, {1, -1}, {-1, 1}, {1, 1}};
    for (int d = 0; d < 4; ++d) {
        const int dx = kDir[d][0];
        const int dy = kDir[d][1];
        const bool hasFaceX = dx < 0 ? hasX[0] : hasX[1];
        const bool hasFaceY = dy < 0 ? hasY[0] : hasY[1];
        if (!hasFaceX || !hasFaceY) {
            continue;
        }
        int diagCoords[3] = {coords[0] + dx, coords[1] + dy, 0};
        for (int a = 2; a < ndims; ++a) {
            diagCoords[a] = coords[a];
        }
        int diagRank = MPI_PROC_NULL;
        if (MPI_Cart_rank(comm, diagCoords, &diagRank) != MPI_SUCCESS) {
            continue;
        }
        const Index si = dx < 0 ? iB : iE - 1;
        const Index sj = dy < 0 ? jB : jE - 1;
        const Index ri = dx < 0 ? iB - 1 : iE;
        const Index rj = dy < 0 ? jB - 1 : jE;
        Real sendVal = data[sj * nxT + si];
        Real recvVal = 0.0;
        MPI_Sendrecv(&sendVal, 1, MPI_DOUBLE, diagRank, kCornerTag,
                     &recvVal, 1, MPI_DOUBLE, diagRank, kCornerTag,
                     comm, MPI_STATUS_IGNORE);
        data[rj * nxT + ri] = recvVal;
    }
}

// ---------------------------------------------------------------------------
// TwoLevelMGSolver
// ---------------------------------------------------------------------------

bool TwoLevelMGSolver::ensureInitialized(Subgrid& subgrid, HaloExchanger& exchanger) {
    (void)exchanger;
    if (valid_ && subgrid.nxLocal() == nxLocal_ && subgrid.nyLocal() == nyLocal_) {
        return true;
    }
    valid_ = false;

    if (subgrid.nzLocal() != 1) {
        HYPOS_ERROR("mg2: 3D grids are not supported (B1a scope is 2D-only)");
        return false;
    }
    if (subgrid.boundaryCondition() != BoundaryCondition::Dirichlet) {
        HYPOS_ERROR("mg2: only Dirichlet BC is supported (the rediscretized "
                    "Neumann coarse operator is singular)");
        return false;
    }

    unsigned long long localEnd[3] = {
        static_cast<unsigned long long>(subgrid.offsetX() + subgrid.nxLocal()),
        static_cast<unsigned long long>(subgrid.offsetY() + subgrid.nyLocal()),
        static_cast<unsigned long long>(subgrid.offsetZ() + subgrid.nzLocal())};
    unsigned long long globalEnd[3] = {0, 0, 0};
    MPI_Allreduce(localEnd, globalEnd, 3, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                  subgrid.comm());
    const unsigned long long nxGlobal = globalEnd[0];
    const unsigned long long nyGlobal = globalEnd[1];
    const unsigned long long nzGlobal = globalEnd[2];

    if (nzGlobal > 1) {
        HYPOS_ERROR("mg2: 3D grids are not supported (B1a scope is 2D-only)");
        return false;
    }
    if (nxGlobal % 2 != 0 || nyGlobal % 2 != 0) {
        HYPOS_ERROR("mg2: even global nx/ny required for coarsening, got "
                    << nxGlobal << " x " << nyGlobal);
        return false;
    }

    nxH_ = static_cast<Index>(nxGlobal / 2);
    nyH_ = static_cast<Index>(nyGlobal / 2);
    if (nxH_ < 4 || nyH_ < 4) {
        HYPOS_ERROR("mg2: coarse grid " << nxH_ << " x " << nyH_
                    << " is below the 4x4 minimum");
        return false;
    }

    int size = 0;
    MPI_Comm_size(subgrid.comm(), &size);
    myBx_ = coarseRangeBegin(subgrid.offsetX(), subgrid.nxLocal());
    myBy_ = coarseRangeBegin(subgrid.offsetY(), subgrid.nyLocal());
    myCnx_ = coarseRangeCount(subgrid.offsetX(), subgrid.nxLocal());
    myCny_ = coarseRangeCount(subgrid.offsetY(), subgrid.nyLocal());

    std::vector<unsigned long long> meta(4 * static_cast<size_t>(size), 0);
    const unsigned long long mine[4] = {static_cast<unsigned long long>(myBx_),
                                        static_cast<unsigned long long>(myBy_),
                                        static_cast<unsigned long long>(myCnx_),
                                        static_cast<unsigned long long>(myCny_)};
    MPI_Allgather(mine, 4, MPI_UNSIGNED_LONG_LONG, meta.data(), 4,
                  MPI_UNSIGNED_LONG_LONG, subgrid.comm());

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
    // Fifth defense beyond the §4.4 "four defense rules": a manual layout
    // whose coarse counts do not tile the coarse grid would make the
    // Allgatherv assembly silently wrong; reject it instead. Recorded in
    // tasks.md (AR007 deviations).
    if (total != nxH_ * nyH_) {
        HYPOS_ERROR("mg2: subdomain layout does not tile the coarse grid ("
                    << total << " vs " << nxH_ * nyH_ << " cells)");
        return false;
    }

    coarsePack_.allocate(std::max<Index>(1, myCnx_ * myCny_));
    coarseAll_.allocate(nxH_ * nyH_);
    fineResidual_.allocate(subgrid.totalCells());
    coarse_ = std::make_unique<Subgrid>(nxH_, nyH_, 1, 1, MPI_COMM_SELF);
    coarse_->setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    coarseExchanger_ = std::make_unique<PointToPointExchanger>();
    coarseExchanger_->initialize(*coarse_);

    nxLocal_ = subgrid.nxLocal();
    nyLocal_ = subgrid.nyLocal();
    valid_ = true;
    return true;
}

// Design deviation from the §4.4 pseudocode: the tolerance lives in the
// coarseTolerance_ member (set once per solve() call) instead of a
// cycleOnce parameter — cycleOnce is also driven by iterate(), which has
// no tolerance argument. Recorded in tasks.md (AR007 deviations).
void TwoLevelMGSolver::cycleOnce(Subgrid& subgrid, HaloExchanger& exchanger) {
    smoother_.smooth(subgrid, exchanger, kPreSmoothSweeps);
    {
        // The smoother's black sweep leaves the neighbor halos of u one
        // half-sweep stale (its last halo exchange follows the red sweep);
        // residualFieldLocal reads those halos, so refresh them here.
        HYPOS_PROFILE("halo_exchange");
        exchanger.exchange(subgrid);
    }

    residualFieldLocal(subgrid, fineResidual_.data());
    {
        HYPOS_PROFILE("halo_exchange");
        exchanger.exchange(subgrid, fineResidual_.data());
    }
    exchangeCorners(subgrid, fineResidual_.data());
    subgrid.applyPhysicalBoundary(fineResidual_.data(), 0.0);

    {
        HYPOS_PROFILE("mg2_restrict");
        restrictResidual(subgrid, fineResidual_.data(), myBx_, myBy_, myCnx_,
                         myCny_, coarsePack_.data());
    }
    MPI_Allgatherv(coarsePack_.data(), static_cast<int>(myCnx_ * myCny_),
                   MPI_DOUBLE, coarseAll_.data(), counts_.data(),
                   displs_.data(), MPI_DOUBLE, subgrid.comm());

    coarse_->zeroInitialize();
    {
        // Coarse equation (design D2, rediscretized 5-point stencil):
        //   ((4 - S)_H / 4) e_H = rho_H, i.e. (4 - S)_H e_H = 4 * rho_H.
        // The factor 4 = (H/h)^2 compensates the discretization-scale
        // mismatch: the fine system is the h^2-scaled Laplacian (4 - S)_h
        // with lowest-mode eigenvalue lambda_f ~ 2 k^2 h^2, while the same
        // stencil pattern rediscretized at H = 2h gives lambda_H = 4 lambda_f.
        // Without the compensation the coarse solve returns e_H ~ e/4
        // (probe-verified 0.258), i.e. every V-cycle under-corrects by 4x
        // and the cycle degenerates to plain smoothing (cycles ~ O(n^2)).
        const Index cnxT = coarse_->nxTotal();
        const Int chw = coarse_->haloWidth();
        Real* crhs = coarse_->rhs().data();
        const int size = static_cast<int>(rectX0_.size());
        for (int r = 0; r < size; ++r) {
            const Real* src = coarseAll_.data() + displs_[r];
            for (Index j = 0; j < rectNY_[r]; ++j) {
                for (Index i = 0; i < rectNX_[r]; ++i) {
                    crhs[(rectY0_[r] + j + chw) * cnxT + rectX0_[r] + i + chw] =
                        -4.0 * src[j * rectNX_[r] + i];
                }
            }
        }
    }

    {
        HYPOS_PROFILE("mg2_coarse_solve");
        coarseIters_ += coarseSolver_.solve(*coarse_, *coarseExchanger_,
                                            nxH_ * nyH_, coarseTolerance_);
    }
    {
        HYPOS_PROFILE("mg2_prolongate");
        prolongateCorrection(subgrid, *coarse_);
    }

    smoother_.smooth(subgrid, exchanger, kPostSmoothSweeps);
}

Index TwoLevelMGSolver::solve(Subgrid& subgrid,
                              HaloExchanger& exchanger,
                              Index maxIter,
                              Real tolerance) {
    Timer totalTimer;
    totalTimer.start();

    coarseIters_ = 0;
    if (!ensureInitialized(subgrid, exchanger)) {
        return 0;
    }
    // Design deviation from D2 ("coarse tol = outer tolerance"): with
    // tolerance == 0 the coarse CG would iterate to machine precision and
    // can break down (pAp == 0) on a converged system; the floor keeps the
    // coarse solve well-posed while preserving the exact-tolerance contract
    // for any tolerance >= 1e-14. Recorded in tasks.md (AR007 deviations).
    coarseTolerance_ = std::max(tolerance, 1e-14);

    Index cycle = 0;
    while (cycle < maxIter) {
        HYPOS_PROFILE("mg2_cycle");
        cycleOnce(subgrid, exchanger);
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
                HYPOS_INFO("Cycle " << cycle << ", residual = " << globalResidual);
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
    HYPOS_INFO("TwoLevelMG finished in " << cycle << " cycles ("
               << cycle * (kPreSmoothSweeps + kPostSmoothSweeps)
               << " smoothing sweeps, " << coarseIters_
               << " coarse CG iterations), time = "
               << totalTimer.elapsedSeconds() << " s");

    return cycle;
}

Real TwoLevelMGSolver::iterate(Subgrid& subgrid, HaloExchanger& exchanger) {
    if (!ensureInitialized(subgrid, exchanger)) {
        return 0.0;
    }
    {
        HYPOS_PROFILE("mg2_cycle");
        cycleOnce(subgrid, exchanger);
    }
    HYPOS_PROFILE("residual_allreduce");
    lastResidual_ = globalTrueResidual(subgrid, exchanger);
    return lastResidual_;
}

} // namespace hypo
