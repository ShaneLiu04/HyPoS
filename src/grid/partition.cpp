#include "grid/partition.hpp"
#include "core/exception.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace hypo {

UniformPartition::UniformPartition(int dimsX, int dimsY, int dimsZ)
    : dimsX_(dimsX), dimsY_(dimsY), dimsZ_(dimsZ) {
}

SubgridInfo UniformPartition::partition(const Grid& grid, MPI_Comm comm, int rank) const {
    int size = 0;
    MPI_Comm_size(comm, &size);

    int dims[3] = {dimsX_, dimsY_, dimsZ_};
    int periods[3] = {0, 0, 0};

    if (dims[0] == 0 && dims[1] == 0 && dims[2] == 0) {
        // Auto-detect optimal dimensions
        if (grid.is2D()) {
            dims[2] = 1;
            MPI_Dims_create(size, 2, dims);
        } else {
            MPI_Dims_create(size, 3, dims);
        }
    }

    // Store dims
    const_cast<UniformPartition*>(this)->dimsX_ = dims[0];
    const_cast<UniformPartition*>(this)->dimsY_ = dims[1];
    const_cast<UniformPartition*>(this)->dimsZ_ = dims[2];

    // Create Cartesian communicator
    MPI_Comm cartComm;
    int reorder = 1;
    MPI_Cart_create(comm, 3, dims, periods, reorder, &cartComm);

    // Get coordinates in Cartesian grid
    int coords[3];
    MPI_Cart_coords(cartComm, rank, 3, coords);

    // Compute local sizes with remainder distribution
    auto divide = [](Index total, int procs, int coord) -> Index {
        Index base = total / procs;
        Index rem = total % procs;
        return base + (coord < static_cast<int>(rem) ? 1 : 0);
    };

    SubgridInfo info;
    info.nxLocal = divide(grid.nx, dims[0], coords[0]);
    info.nyLocal = divide(grid.ny, dims[1], coords[1]);
    info.nzLocal = divide(grid.nz, dims[2], coords[2]);
    info.rank = rank;

    // Compute offsets
    info.offsetX = 0;
    for (int i = 0; i < coords[0]; ++i) {
        info.offsetX += divide(grid.nx, dims[0], i);
    }
    info.offsetY = 0;
    for (int j = 0; j < coords[1]; ++j) {
        info.offsetY += divide(grid.ny, dims[1], j);
    }
    info.offsetZ = 0;
    for (int k = 0; k < coords[2]; ++k) {
        info.offsetZ += divide(grid.nz, dims[2], k);
    }

    MPI_Comm_free(&cartComm);
    return info;
}

void UniformPartition::getTopologyDims(int& dimsX, int& dimsY, int& dimsZ) const {
    dimsX = dimsX_;
    dimsY = dimsY_;
    dimsZ = dimsZ_;
}

} // namespace hypo
