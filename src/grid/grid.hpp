#pragma once

#include "core/types.hpp"

namespace hypo {

/**
 * @brief Global 2D/3D Cartesian grid descriptor.
 */
struct Grid {
    Index nx = 1024;  ///< Global cells in x
    Index ny = 1024;  ///< Global cells in y
    Index nz = 1;     ///< Global cells in z (1 for 2D)
    Real  dx = 1.0;   ///< Cell size in x
    Real  dy = 1.0;   ///< Cell size in y
    Real  dz = 1.0;   ///< Cell size in z
    Int   haloWidth = 1;

    bool is2D() const noexcept { return nz == 1; }
    bool is3D() const noexcept { return nz > 1; }

    Index totalCells() const noexcept { return nx * ny * nz; }
};

} // namespace hypo
