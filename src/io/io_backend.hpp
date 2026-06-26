#pragma once

#include "core/types.hpp"
#include "grid/subgrid.hpp"
#include <string>

namespace hypo {

/**
 * @brief Abstract interface for I/O backends.
 */
class IOBackend {
public:
    virtual ~IOBackend() = default;

    /**
     * @brief Write the current solution field to file.
     * @param subgrid Local subdomain data.
     * @param filename Output filename (may be rank-qualified).
     * @param step Iteration number (for time-series output).
     */
    virtual void write(const Subgrid& subgrid, const std::string& filename, int step = 0) = 0;

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
 * @brief VTK legacy format output for ParaView visualization.
 * Writes structured grid points with scalar data.
 */
class VTKIOBackend : public IOBackend {
public:
    VTKIOBackend(Real dx, Real dy, Real dz = 1.0);

    void write(const Subgrid& subgrid, const std::string& filename, int step = 0) override;
    std::string format() const override { return "vtk"; }

private:
    Real dx_ = 1.0, dy_ = 1.0, dz_ = 1.0;
};

} // namespace hypo
