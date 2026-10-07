#include "io/io_backend.hpp"
#include "utils/logger.hpp"
#include <fstream>
#include <iomanip>
#include <string>

namespace hypo {

// ============================================================================
// VTKIOBackend (XML ImageData pieces + PImageData index)
// ============================================================================

VTKIOBackend::VTKIOBackend(Real dx, Real dy, Real dz)
    : dx_(dx), dy_(dy), dz_(dz) {
}

void VTKIOBackend::write(const Subgrid& subgrid, const std::string& filename, int step) {
    (void)step;
    std::string fname = filename;
    if (fname.find(".vti") == std::string::npos) {
        fname += ".vti";
    }

    std::ofstream ofs(fname);
    if (!ofs) {
        HYPOS_WARN("Cannot open VTK output file: " << fname);
        return;
    }

    const Index nx = subgrid.nxLocal();
    const Index ny = subgrid.nyLocal();
    const Index nz = subgrid.nzLocal();
    const Index x0 = subgrid.offsetX();
    const Index x1 = x0 + nx - 1;
    const Index y0 = subgrid.offsetY();
    const Index y1 = y0 + ny - 1;
    const Index z0 = subgrid.offsetZ();
    const Index z1 = z0 + nz - 1;

    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    ofs << "  <ImageData WholeExtent=\"" << x0 << " " << x1 << " " << y0 << " " << y1
        << " " << z0 << " " << z1 << "\" Origin=\""
        << static_cast<Real>(x0) * dx_ << " " << static_cast<Real>(y0) * dy_ << " "
        << static_cast<Real>(z0) * dz_ << "\" Spacing=\""
        << dx_ << " " << dy_ << " " << dz_ << "\">\n";
    ofs << "    <Piece Extent=\"" << x0 << " " << x1 << " " << y0 << " " << y1
        << " " << z0 << " " << z1 << "\">\n";
    ofs << "      <PointData Scalars=\"u\">\n";
    ofs << "        <DataArray type=\"Float64\" Name=\"u\" format=\"ascii\">\n";

    const Real* u = subgrid.u().data();
    ofs << std::scientific << std::setprecision(12);
    for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                ofs << u[subgrid.index(i, j, k)] << " ";
            }
        }
    }
    ofs << "\n";
    ofs << "        </DataArray>\n";
    ofs << "      </PointData>\n";
    ofs << "    </Piece>\n";
    ofs << "  </ImageData>\n";
    ofs << "</VTKFile>\n";

    HYPOS_INFO("Wrote VTK piece to " << fname);
}

void VTKIOBackend::writeParallelIndex(const Grid& grid,
                                      const std::string& baseName,
                                      int step,
                                      const std::vector<PieceExtent>& pieces) {
    (void)step;
    if (pieces.empty()) {
        HYPOS_WARN("VTK parallel index skipped: empty piece list");
        return;
    }
    for (const PieceExtent& p : pieces) {
        if (p.nxLocal == 0 || p.nyLocal == 0 || p.nzLocal == 0 ||
            p.offsetX + p.nxLocal > grid.nx ||
            p.offsetY + p.nyLocal > grid.ny ||
            p.offsetZ + p.nzLocal > grid.nz) {
            HYPOS_WARN("VTK parallel index skipped: piece extents exceed the global grid");
            return;
        }
    }

    const std::size_t slash = baseName.find_last_of("/\\");
    const std::string leaf = (slash == std::string::npos) ? baseName : baseName.substr(slash + 1);

    std::ofstream ofs(baseName + ".pvti");
    if (!ofs) {
        HYPOS_WARN("Cannot open VTK parallel index file: " << baseName << ".pvti");
        return;
    }

    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"PImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    ofs << "  <PImageData WholeExtent=\"0 " << (grid.nx - 1) << " 0 " << (grid.ny - 1)
        << " 0 " << (grid.nz - 1) << "\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\""
        << dx_ << " " << dy_ << " " << dz_ << "\">\n";
    ofs << "    <PPointData Scalars=\"u\">\n";
    ofs << "      <PDataArray type=\"Float64\" Name=\"u\"/>\n";
    ofs << "    </PPointData>\n";

    for (std::size_t r = 0; r < pieces.size(); ++r) {
        const PieceExtent& p = pieces[r];
        ofs << "    <Piece Extent=\"" << p.offsetX << " " << (p.offsetX + p.nxLocal - 1)
            << " " << p.offsetY << " " << (p.offsetY + p.nyLocal - 1)
            << " " << p.offsetZ << " " << (p.offsetZ + p.nzLocal - 1)
            << "\" Source=\"" << leaf << "_r" << r << ".vti\"/>\n";
    }
    ofs << "  </PImageData>\n";
    ofs << "</VTKFile>\n";

    HYPOS_INFO("Wrote VTK parallel index to " << baseName << ".pvti");
}

} // namespace hypo
