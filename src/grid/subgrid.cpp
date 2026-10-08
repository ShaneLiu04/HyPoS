#include "grid/subgrid.hpp"
#include "core/exception.hpp"
#include <algorithm>
#include <cmath>

namespace hypo {

Subgrid::Subgrid(Index nxLocal, Index nyLocal, Index nzLocal, Int haloWidth, MPI_Comm comm)
    : nxLocal_(nxLocal)
    , nyLocal_(nyLocal)
    , nzLocal_(nzLocal)
    , haloWidth_(haloWidth)
    , comm_(comm)
{
    Index total = totalCells();
    u_.allocate(total);
    uNext_.allocate(total);
    rhs_.allocate(total);
    zeroInitialize();
}

void Subgrid::setNeighbors(int left, int right, int down, int up, int back, int front) noexcept {
    neighborLeft_  = left;
    neighborRight_ = right;
    neighborDown_  = down;
    neighborUp_    = up;
    neighborBack_  = back;
    neighborFront_ = front;
}

void Subgrid::zeroInitialize() noexcept {
    // AR004 A1 first-touch: zero the ENTIRE buffer of all three fields in a
    // single pass, parallel over the outer dimension with schedule(static)
    // to match the stencil loops' partitioning. A 2D subgrid still
    // allocates nzTotal() z-planes (halo on both sides), so the 2D branch
    // keeps an inner k loop; parallelizing over j there preserves thread
    // utilization (nzTotal is tiny for 2D).
    Real* u = u_.data();
    Real* un = uNext_.data();
    Real* r = rhs_.data();
    const Index nxT = nxTotal();
    const Index nyT = nyTotal();
    const Index nzT = nzTotal();

    if (nzLocal_ == 1) {
        #pragma omp parallel for schedule(static)
        for (Index j = 0; j < nyT; ++j) {
            for (Index k = 0; k < nzT; ++k) {
                #pragma omp simd
                for (Index i = 0; i < nxT; ++i) {
                    const Index idx = (k * nyT + j) * nxT + i;
                    u[idx] = 0.0;
                    un[idx] = 0.0;
                    r[idx] = 0.0;
                }
            }
        }
    } else {
        #pragma omp parallel for schedule(static)
        for (Index k = 0; k < nzT; ++k) {
            for (Index j = 0; j < nyT; ++j) {
                #pragma omp simd
                for (Index i = 0; i < nxT; ++i) {
                    const Index idx = (k * nyT + j) * nxT + i;
                    u[idx] = 0.0;
                    un[idx] = 0.0;
                    r[idx] = 0.0;
                }
            }
        }
    }
}

void Subgrid::applyDirichletBC(Real fixedValue) noexcept {
    // Left halo (i = 0..halo-1)
    for (Index k = 0; k < nzTotal(); ++k) {
        for (Index j = 0; j < nyTotal(); ++j) {
            for (Index i = 0; i < static_cast<Index>(haloWidth_); ++i) {
                u_[index(i, j, k)] = fixedValue;
            }
        }
    }
    // Right halo (i = nxLocal+halo..nxTotal-1)
    for (Index k = 0; k < nzTotal(); ++k) {
        for (Index j = 0; j < nyTotal(); ++j) {
            for (Index i = nxLocal_ + haloWidth_; i < nxTotal(); ++i) {
                u_[index(i, j, k)] = fixedValue;
            }
        }
    }
    // Bottom halo (j = 0..halo-1)
    for (Index k = 0; k < nzTotal(); ++k) {
        for (Index j = 0; j < static_cast<Index>(haloWidth_); ++j) {
            for (Index i = 0; i < nxTotal(); ++i) {
                u_[index(i, j, k)] = fixedValue;
            }
        }
    }
    // Top halo (j = nyLocal+halo..nyTotal-1)
    for (Index k = 0; k < nzTotal(); ++k) {
        for (Index j = nyLocal_ + haloWidth_; j < nyTotal(); ++j) {
            for (Index i = 0; i < nxTotal(); ++i) {
                u_[index(i, j, k)] = fixedValue;
            }
        }
    }
    // Back halo (k = 0..halo-1) for 3D
    if (nzLocal_ > 1) {
        for (Index k = 0; k < static_cast<Index>(haloWidth_); ++k) {
            for (Index j = 0; j < nyTotal(); ++j) {
                for (Index i = 0; i < nxTotal(); ++i) {
                    u_[index(i, j, k)] = fixedValue;
                }
            }
        }
    }
    // Front halo (k = nzLocal+halo..nzTotal-1) for 3D
    if (nzLocal_ > 1) {
        for (Index k = nzLocal_ + haloWidth_; k < nzTotal(); ++k) {
            for (Index j = 0; j < nyTotal(); ++j) {
                for (Index i = 0; i < nxTotal(); ++i) {
                    u_[index(i, j, k)] = fixedValue;
                }
            }
        }
    }
}

void Subgrid::applyPhysicalBoundary(Real dirichletValue) noexcept {
    applyPhysicalBoundary(u_.data(), dirichletValue);
}

void Subgrid::applyPhysicalBoundary(Real* data, Real dirichletValue) noexcept {
    const bool neumann = (boundaryCondition_ == BoundaryCondition::Neumann);

    auto clampInto = [](Index v, Index lo, Index hi) -> Index {
        if (v < lo) {
            return lo;
        }
        if (v > hi) {
            return hi;
        }
        return v;
    };

    const Index kSourceLo = kBegin();
    const Index kSourceHi = (nzLocal_ == 1) ? 0 : (kEnd() - 1);

    // Left / right faces (i halo bands), full transverse band.
    for (Index k = 0; k < nzTotal(); ++k) {
        const Index kSrc = (nzLocal_ == 1) ? 0 : clampInto(k, kBegin(), kEnd() - 1);
        for (Index j = 0; j < nyTotal(); ++j) {
            const Index jSrc = clampInto(j, jBegin(), jEnd() - 1);
            if (neighborLeft_ == MPI_PROC_NULL) {
                for (Index i = 0; i < static_cast<Index>(haloWidth_); ++i) {
                    data[index(i, j, k)] = neumann ? data[index(iBegin(), jSrc, kSrc)] : dirichletValue;
                }
            }
            if (neighborRight_ == MPI_PROC_NULL) {
                for (Index i = iEnd(); i < nxTotal(); ++i) {
                    data[index(i, j, k)] = neumann ? data[index(iEnd() - 1, jSrc, kSrc)] : dirichletValue;
                }
            }
        }
    }

    // Bottom / top faces (j halo bands), full x band.
    for (Index k = 0; k < nzTotal(); ++k) {
        const Index kSrc = (nzLocal_ == 1) ? 0 : clampInto(k, kBegin(), kEnd() - 1);
        for (Index i = 0; i < nxTotal(); ++i) {
            const Index iSrc = clampInto(i, iBegin(), iEnd() - 1);
            if (neighborDown_ == MPI_PROC_NULL) {
                for (Index j = 0; j < static_cast<Index>(haloWidth_); ++j) {
                    data[index(i, j, k)] = neumann ? data[index(iSrc, jBegin(), kSrc)] : dirichletValue;
                }
            }
            if (neighborUp_ == MPI_PROC_NULL) {
                for (Index j = jEnd(); j < nyTotal(); ++j) {
                    data[index(i, j, k)] = neumann ? data[index(iSrc, jEnd() - 1, kSrc)] : dirichletValue;
                }
            }
        }
    }

    // Back / front faces (k halo bands, 3D only; 2D has no z boundary).
    if (nzLocal_ > 1) {
        for (Index j = 0; j < nyTotal(); ++j) {
            const Index jSrc = clampInto(j, jBegin(), jEnd() - 1);
            for (Index i = 0; i < nxTotal(); ++i) {
                const Index iSrc = clampInto(i, iBegin(), iEnd() - 1);
                if (neighborBack_ == MPI_PROC_NULL) {
                    for (Index k = 0; k < static_cast<Index>(haloWidth_); ++k) {
                        data[index(i, j, k)] = neumann ? data[index(iSrc, jSrc, kSourceLo)] : dirichletValue;
                    }
                }
                if (neighborFront_ == MPI_PROC_NULL) {
                    for (Index k = kEnd(); k < nzTotal(); ++k) {
                        data[index(i, j, k)] = neumann ? data[index(iSrc, jSrc, kSourceHi)] : dirichletValue;
                    }
                }
            }
        }
    }
}

} // namespace hypo
