#pragma once

#include "core/types.hpp"
#include "grid/grid.hpp"
#include "grid/subgrid.hpp"
#include <string>
#include <vector>

namespace hypo {

/**
 * @brief Extent of one rank's interior domain within the global grid
 * (0-based global point indices, used by parallel file formats).
 */
struct PieceExtent {
    Index offsetX = 0, offsetY = 0, offsetZ = 0;
    Index nxLocal = 0, nyLocal = 0, nzLocal = 0;
};

/**
 * @brief Abstract interface for I/O backends.
 */
class IOBackend {
public:
    virtual ~IOBackend() = default;

    /**
     * @brief Write the current solution field to file.
     * @param subgrid Local subdomain data (offsets define its global position).
     * @param filename Output filename (may be rank-qualified).
     * @param step Iteration number (for time-series output).
     */
    virtual void write(const Subgrid& subgrid, const std::string& filename, int step = 0) = 0;

    /**
     * @brief Optional: write a parallel-format index referencing all pieces.
     * Only the root rank should call this; default implementation is a no-op.
     * @param grid Global grid description.
     * @param baseName Common base path (pieces are "<baseName>_r<i>.<ext>").
     * @param step Iteration number.
     * @param pieces Per-rank extents (ordered by rank).
     */
    virtual void writeParallelIndex(const Grid& grid,
                                    const std::string& baseName,
                                    int step,
                                    const std::vector<PieceExtent>& pieces) {
        (void)grid;
        (void)baseName;
        (void)step;
        (void)pieces;
    }

    virtual std::string format() const = 0;
};

/**
 * @brief Simple binary output (raw float64 arrays).
 */
class BinaryIOBackend : public IOBackend {
public:
    void write(const Subgrid& subgrid, const std::string& filename, int step = 0) override;
    std::string format() const override { return "binary"; }
};

/**
 * @brief VTK XML output for ParaView visualization.
 * Each rank writes its own `.vti` (ImageData) piece carrying the global
 * origin/extent; rank 0 can additionally write a `.pvti` index.
 */
class VTKIOBackend : public IOBackend {
public:
    VTKIOBackend(Real dx, Real dy, Real dz = 1.0);

    void write(const Subgrid& subgrid, const std::string& filename, int step = 0) override;
    void writeParallelIndex(const Grid& grid,
                            const std::string& baseName,
                            int step,
                            const std::vector<PieceExtent>& pieces) override;
    std::string format() const override { return "vtk"; }

private:
    Real dx_ = 1.0, dy_ = 1.0, dz_ = 1.0;
};

} // namespace hypo
