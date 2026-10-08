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

    // Write interior data as row-major (x fastest) memory image. Interior
    // cells along x are contiguous (row-major layout, no padding), so each
    // (k, j) row is emitted with a single block write instead of one call
    // per element. File byte layout is identical to the elementwise loop.
    const Real* u = subgrid.u().data();
    const std::streamsize rowBytes = static_cast<std::streamsize>(nx) * sizeof(Real);
    for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            ofs.write(reinterpret_cast<const char*>(&u[subgrid.index(subgrid.iBegin(), j, k)]),
                      rowBytes);
        }
    }

    HYPOS_INFO("Wrote binary output to " << fname);
}

} // namespace hypo
