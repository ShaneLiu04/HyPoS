#include "io/io_backend.hpp"
#include "utils/logger.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>

namespace hypo {

// ============================================================================
// BinaryIOBackend
// ============================================================================

void BinaryIOBackend::write(const Subgrid& subgrid, const std::string& filename, int step) {
    std::string fname = filename;
    if (fname.find(".") == std::string::npos) {
        fname += ".bin";
    }

    std::ofstream ofs(fname, std::ios::binary);
    if (!ofs) {
        HYPOS_WARN("Cannot open binary output file: " << fname);
        return;
    }

    // Write header: nx, ny, nz, halo, offsetX, offsetY, offsetZ as Index (8 bytes each)
    Index nx = subgrid.nxLocal();
    Index ny = subgrid.nyLocal();
    Index nz = subgrid.nzLocal();
    Index hw = subgrid.haloWidth();
    Index offsetX = 0, offsetY = 0, offsetZ = 0;

    ofs.write(reinterpret_cast<const char*>(&nx), sizeof(Index));
    ofs.write(reinterpret_cast<const char*>(&ny), sizeof(Index));
    ofs.write(reinterpret_cast<const char*>(&nz), sizeof(Index));
    ofs.write(reinterpret_cast<const char*>(&hw), sizeof(Index));
    ofs.write(reinterpret_cast<const char*>(&offsetX), sizeof(Index));
    ofs.write(reinterpret_cast<const char*>(&offsetY), sizeof(Index));
    ofs.write(reinterpret_cast<const char*>(&offsetZ), sizeof(Index));

    // Write interior data only
    const Real* u = subgrid.u().data();
    for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                Real val = u[subgrid.index(i, j, k)];
                ofs.write(reinterpret_cast<const char*>(&val), sizeof(Real));
            }
        }
    }

    HYPOS_INFO("Wrote binary output to " << fname);
}

// ============================================================================
// VTKIOBackend
// ============================================================================

VTKIOBackend::VTKIOBackend(Real dx, Real dy, Real dz)
    : dx_(dx), dy_(dy), dz_(dz) {
}

void VTKIOBackend::write(const Subgrid& subgrid, const std::string& filename, int step) {
    std::string fname = filename;
    if (fname.find(".") == std::string::npos) {
        fname += "_" + std::to_string(step) + ".vtk";
    }

    std::ofstream ofs(fname);
    if (!ofs) {
        HYPOS_WARN("Cannot open VTK output file: " << fname);
        return;
    }

    Index nx = subgrid.nxLocal();
    Index ny = subgrid.nyLocal();
    Index nz = subgrid.nzLocal();

    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "HyPoS Output\n";
    ofs << "ASCII\n";
    ofs << "DATASET STRUCTURED_POINTS\n";
    ofs << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    ofs << "ORIGIN 0 0 0\n";
    ofs << "SPACING " << dx_ << " " << dy_ << " " << dz_ << "\n";
    ofs << "POINT_DATA " << (nx * ny * nz) << "\n";
    ofs << "SCALARS u double 1\n";
    ofs << "LOOKUP_TABLE default\n";

    const Real* u = subgrid.u().data();
    for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                ofs << std::scientific << std::setprecision(12)
                    << u[subgrid.index(i, j, k)] << "\n";
            }
        }
    }

    HYPOS_INFO("Wrote VTK output to " << fname);
}

} // namespace hypo
