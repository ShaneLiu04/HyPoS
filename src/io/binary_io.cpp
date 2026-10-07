#include "io/io_backend.hpp"
#include "utils/logger.hpp"
#include <fstream>

namespace hypo {

// ============================================================================
// BinaryIOBackend
// ============================================================================

void BinaryIOBackend::write(const Subgrid& subgrid, const std::string& filename, int step) {
    (void)step;
    std::string fname = filename;
    const std::string suffix = ".bin";
    if (fname.size() < suffix.size() ||
        fname.compare(fname.size() - suffix.size(), suffix.size(), suffix) != 0) {
        fname += suffix;
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
    Index offsetX = subgrid.offsetX();
    Index offsetY = subgrid.offsetY();
    Index offsetZ = subgrid.offsetZ();

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

} // namespace hypo
