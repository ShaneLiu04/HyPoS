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
    u_.fill(0.0);
    uNext_.fill(0.0);
    rhs_.fill(0.0);
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

void Subgrid::applyNeumannBC() noexcept {
    // Left: copy from first interior column
    for (Index k = 0; k < nzTotal(); ++k) {
        for (Index j = 0; j < nyTotal(); ++j) {
            for (Index i = 0; i < static_cast<Index>(haloWidth_); ++i) {
                u_[index(i, j, k)] = u_[index(iBegin(), j, k)];
            }
        }
    }
    // Right: copy from last interior column
    for (Index k = 0; k < nzTotal(); ++k) {
        for (Index j = 0; j < nyTotal(); ++j) {
            for (Index i = nxLocal_ + haloWidth_; i < nxTotal(); ++i) {
                u_[index(i, j, k)] = u_[index(iEnd() - 1, j, k)];
            }
        }
    }
    // Bottom: copy from first interior row
    for (Index k = 0; k < nzTotal(); ++k) {
        for (Index j = 0; j < static_cast<Index>(haloWidth_); ++j) {
            for (Index i = 0; i < nxTotal(); ++i) {
                u_[index(i, j, k)] = u_[index(i, jBegin(), k)];
            }
        }
    }
    // Top: copy from last interior row
    for (Index k = 0; k < nzTotal(); ++k) {
        for (Index j = nyLocal_ + haloWidth_; j < nyTotal(); ++j) {
            for (Index i = 0; i < nxTotal(); ++i) {
                u_[index(i, j, k)] = u_[index(i, jEnd() - 1, k)];
            }
        }
    }
    // Back/front for 3D
    if (nzLocal_ > 1) {
        for (Index k = 0; k < static_cast<Index>(haloWidth_); ++k) {
            for (Index j = 0; j < nyTotal(); ++j) {
                for (Index i = 0; i < nxTotal(); ++i) {
                    u_[index(i, j, k)] = u_[index(i, j, kBegin())];
                }
            }
        }
        for (Index k = nzLocal_ + haloWidth_; k < nzTotal(); ++k) {
            for (Index j = 0; j < nyTotal(); ++j) {
                for (Index i = 0; i < nxTotal(); ++i) {
                    u_[index(i, j, k)] = u_[index(i, j, kEnd() - 1)];
                }
            }
        }
    }
}

} // namespace hypo
