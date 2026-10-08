#include "io/io_backend.hpp"
#include "utils/logger.hpp"
#include <mpi.h>

#include <cstdint>
#include <cstring>
#include <string>

namespace hypo {

namespace {

// File layout (design.md AR003 §4.1.2), little-endian, fixed widths:
//   [0..4)   magic "HYPS"
//   [4..8)   version (uint32)
//   [8..16)  globalNx (uint64)      [16..24) globalNy   [24..32) globalNz
//   [32..40) dx (double)            [40..48) dy         [48..56) dz
//   [56..60) bcType (uint32: 0=dirichlet, 1=neumann)
//   [60..64) reserved (uint32, zero)
//   [64..72) dataOffset (uint64, = 72)
constexpr char kMagic[4] = {'H', 'Y', 'P', 'S'};
constexpr std::uint32_t kVersion = 1;
constexpr std::uint64_t kHeaderBytes = 72;

void putU32(unsigned char* p, std::uint32_t v) {
    std::memcpy(p, &v, sizeof(v));
}

void putU64(unsigned char* p, std::uint64_t v) {
    std::memcpy(p, &v, sizeof(v));
}

void putF64(unsigned char* p, double v) {
    std::memcpy(p, &v, sizeof(v));
}

} // namespace

// ============================================================================
// MPIIOBinaryBackend
// ============================================================================

MPIIOBinaryBackend::MPIIOBinaryBackend(const Grid& grid)
    : grid_(grid) {}

void MPIIOBinaryBackend::write(const Subgrid& subgrid, const std::string& filename, int step) {
    (void)step;

    std::string fname = filename;
    const std::string suffix = ".bin";
    if (fname.size() < suffix.size() ||
        fname.compare(fname.size() - suffix.size(), suffix.size(), suffix) != 0) {
        fname += suffix;
    }

    // Collective open: with the default file error handler (MPI_ERRORS_RETURN)
    // a failure such as a missing directory is reported on ALL ranks, so the
    // early return below keeps every rank on the same collective path.
    MPI_File fh = MPI_FILE_NULL;
    int err = MPI_File_open(subgrid.comm(), fname.c_str(),
                            MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    if (err != MPI_SUCCESS) {
        if (!warned_) {
            HYPOS_WARN("MPI-IO: cannot open output file: " << fname);
            warned_ = true;
        }
        return;
    }

    // Exact final size; truncates any stale longer file left by earlier runs
    // (MPI_MODE_CREATE does not imply truncate).
    const MPI_Offset dataBytes =
        static_cast<MPI_Offset>(grid_.totalCells()) * static_cast<MPI_Offset>(sizeof(Real));
    const MPI_Offset fileSize = static_cast<MPI_Offset>(kHeaderBytes) + dataBytes;
    err = MPI_File_set_size(fh, fileSize);
    if (err != MPI_SUCCESS && !warned_) {
        HYPOS_WARN("MPI-IO: cannot set file size for: " << fname);
        warned_ = true;
    }

    // File-side view: this rank's interior box inside the global grid. Grid
    // dims are bounded by CLI config (<< 2^31), so the int casts are safe.
    // NOTE: MPI_ORDER_C expects dimensions in C declaration order, i.e.
    // slowest-varying first — for our row-major x-fastest layout that is
    // {z, y, x}.
    int gsize[3] = {static_cast<int>(grid_.nz), static_cast<int>(grid_.ny),
                    static_cast<int>(grid_.nx)};
    int lsize[3] = {static_cast<int>(subgrid.nzLocal()), static_cast<int>(subgrid.nyLocal()),
                    static_cast<int>(subgrid.nxLocal())};
    int lstart[3] = {static_cast<int>(subgrid.offsetZ()), static_cast<int>(subgrid.offsetY()),
                     static_cast<int>(subgrid.offsetX())};

    MPI_Datatype fileType = MPI_DATATYPE_NULL;
    MPI_Type_create_subarray(3, gsize, lsize, lstart, MPI_ORDER_C, MPI_DOUBLE, &fileType);
    MPI_Type_commit(&fileType);

    // Memory-side view: the same interior box inside the halo-padded buffer,
    // so ghost cells never reach the file. Same C-order dimension rule.
    int msize[3] = {static_cast<int>(subgrid.nzTotal()), static_cast<int>(subgrid.nyTotal()),
                    static_cast<int>(subgrid.nxTotal())};
    int mstart[3] = {static_cast<int>(subgrid.haloWidth()), static_cast<int>(subgrid.haloWidth()),
                     static_cast<int>(subgrid.haloWidth())};

    MPI_Datatype memType = MPI_DATATYPE_NULL;
    MPI_Type_create_subarray(3, msize, lsize, mstart, MPI_ORDER_C, MPI_DOUBLE, &memType);
    MPI_Type_commit(&memType);

    // Header is written by rank 0 at absolute byte 0 BEFORE any view is
    // installed: MPI_File_write_at offsets are view-relative (in etypes), so
    // with the default view (whole file, etype MPI_BYTE) offset 0 is the
    // file start. write_at is local, so ranks != 0 proceed straight to the
    // collective set_view below — no ordering hazard.
    int rank = 0;
    MPI_Comm_rank(subgrid.comm(), &rank);
    if (rank == 0) {
        unsigned char header[kHeaderBytes];
        std::memset(header, 0, sizeof(header));
        std::memcpy(header, kMagic, sizeof(kMagic));
        putU32(header + 4, kVersion);
        putU64(header + 8, static_cast<std::uint64_t>(grid_.nx));
        putU64(header + 16, static_cast<std::uint64_t>(grid_.ny));
        putU64(header + 24, static_cast<std::uint64_t>(grid_.nz));
        putF64(header + 32, grid_.dx);
        putF64(header + 40, grid_.dy);
        putF64(header + 48, grid_.dz);
        putU32(header + 56, subgrid.boundaryCondition() == BoundaryCondition::Neumann ? 1u : 0u);
        // header + 60: reserved, stays zero
        putU64(header + 64, kHeaderBytes);

        err = MPI_File_write_at(fh, 0, header, static_cast<int>(kHeaderBytes), MPI_BYTE,
                                MPI_STATUS_IGNORE);
        if (err != MPI_SUCCESS && !warned_) {
            HYPOS_WARN("MPI-IO: header write failed for: " << fname);
            warned_ = true;
        }
    }

    MPI_File_set_view(fh, static_cast<MPI_Offset>(kHeaderBytes), MPI_DOUBLE, fileType,
                      "native", MPI_INFO_NULL);

    err = MPI_File_write_all(fh, subgrid.u().data(), 1, memType, MPI_STATUS_IGNORE);
    if (err != MPI_SUCCESS && !warned_) {
        HYPOS_WARN("MPI-IO: collective data write failed for: " << fname);
        warned_ = true;
    }

    MPI_Type_free(&fileType);
    MPI_Type_free(&memType);
    MPI_File_close(&fh);
}

} // namespace hypo
