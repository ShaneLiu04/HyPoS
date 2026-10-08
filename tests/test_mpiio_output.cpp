#include <gtest/gtest.h>
#include <mpi.h>

#include "grid/grid.hpp"
#include "grid/subgrid.hpp"
#include "io/io_backend.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

using namespace hypo;

namespace {

// File layout constants (design.md §4.1.2), defined independently here:
// tests decode the layout by offset instead of including implementation
// internals.
constexpr std::size_t kHeaderBytes = 72;
constexpr uint32_t kVersion = 1;

struct FileHeader {
    char magic[4] = {0, 0, 0, 0};
    uint32_t version = 0;
    uint64_t globalNx = 0;
    uint64_t globalNy = 0;
    uint64_t globalNz = 0;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    uint32_t bcType = 0;
    uint32_t reserved = 1; // init to non-zero so stale values cannot pass
    uint64_t dataOffset = 0;
};

uint32_t decodeU32(const unsigned char* p) {
    uint32_t v = 0;
    std::memcpy(&v, p, sizeof(v));
    return v;
}

uint64_t decodeU64(const unsigned char* p) {
    uint64_t v = 0;
    std::memcpy(&v, p, sizeof(v));
    return v;
}

double decodeF64(const unsigned char* p) {
    double v = 0.0;
    std::memcpy(&v, p, sizeof(v));
    return v;
}

bool readFileBytes(const std::string& path, std::vector<unsigned char>& out) {
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) return false;
    f.seekg(0, std::ios::end);
    const std::streamoff len = f.tellg();
    if (len < 0) return false;
    out.resize(static_cast<std::size_t>(len));
    if (len > 0) {
        f.seekg(0, std::ios::beg);
        f.read(reinterpret_cast<char*>(out.data()), len);
        if (!f.good() && !f.eof()) return false;
    }
    return true;
}

// Decode the fixed 72-byte header from an already-validated buffer.
FileHeader parseHeader(const std::vector<unsigned char>& buf) {
    FileHeader h;
    std::memcpy(h.magic, buf.data(), 4);
    h.version = decodeU32(buf.data() + 4);
    h.globalNx = decodeU64(buf.data() + 8);
    h.globalNy = decodeU64(buf.data() + 16);
    h.globalNz = decodeU64(buf.data() + 24);
    h.dx = decodeF64(buf.data() + 32);
    h.dy = decodeF64(buf.data() + 40);
    h.dz = decodeF64(buf.data() + 48);
    h.bcType = decodeU32(buf.data() + 56);
    h.reserved = decodeU32(buf.data() + 60);
    h.dataOffset = decodeU64(buf.data() + 64);
    return h;
}

// Compare every double in the data area against the reference function.
// All reference values are exactly representable integers, so a plain
// bitwise == comparison is exact.
bool dataAreaMatches(const std::vector<unsigned char>& buf,
                     const FileHeader& h,
                     const std::function<double(uint64_t, uint64_t, uint64_t)>& ref) {
    if (buf.size() < kHeaderBytes) return false;
    if (h.dataOffset != kHeaderBytes) return false;
    const uint64_t n = h.globalNx * h.globalNy * h.globalNz;
    if (buf.size() != kHeaderBytes + static_cast<std::size_t>(n) * sizeof(double)) return false;
    for (uint64_t k = 0; k < h.globalNz; ++k) {
        for (uint64_t j = 0; j < h.globalNy; ++j) {
            for (uint64_t i = 0; i < h.globalNx; ++i) {
                const uint64_t idx = i + j * h.globalNx + k * h.globalNx * h.globalNy;
                const double actual = decodeF64(buf.data() + kHeaderBytes +
                                                static_cast<std::size_t>(idx) * sizeof(double));
                if (actual != ref(i, j, k)) return false;
            }
        }
    }
    return true;
}

// Reference fields: f(gI,gJ[,gK]) built from exact-integer doubles.
double ref2D(uint64_t gI, uint64_t gJ, uint64_t) {
    return static_cast<double>(gI) + 1000.0 * static_cast<double>(gJ);
}

double ref3D(uint64_t gI, uint64_t gJ, uint64_t gK) {
    return static_cast<double>(gI) + 1000.0 * static_cast<double>(gJ) +
           1000000.0 * static_cast<double>(gK);
}

bool headerMatches(const FileHeader& h,
                   uint64_t gnx, uint64_t gny, uint64_t gnz,
                   double dx, double dy, double dz,
                   uint32_t bcType) {
    return h.magic[0] == 'H' && h.magic[1] == 'Y' && h.magic[2] == 'P' && h.magic[3] == 'S' &&
           h.version == kVersion &&
           h.globalNx == gnx && h.globalNy == gny && h.globalNz == gnz &&
           h.dx == dx && h.dy == dy && h.dz == dz &&
           h.bcType == bcType && h.reserved == 0 &&
           h.dataOffset == kHeaderBytes;
}

void expectHeaderMatches(const FileHeader& h,
                         uint64_t gnx, uint64_t gny, uint64_t gnz,
                         double dx, double dy, double dz,
                         uint32_t bcType) {
    EXPECT_EQ(0, std::memcmp(h.magic, "HYPS", 4));
    EXPECT_EQ(kVersion, h.version);
    EXPECT_EQ(gnx, h.globalNx);
    EXPECT_EQ(gny, h.globalNy);
    EXPECT_EQ(gnz, h.globalNz);
    EXPECT_EQ(dx, h.dx);
    EXPECT_EQ(dy, h.dy);
    EXPECT_EQ(dz, h.dz);
    EXPECT_EQ(bcType, h.bcType);
    EXPECT_EQ(0u, h.reserved);
    EXPECT_EQ(static_cast<uint64_t>(kHeaderBytes), h.dataOffset);
}

// Fill interior cells of u with the 2D reference field given the subdomain
// global offsets.
void fill2D(Subgrid& sg, Index offsetX, Index offsetY) {
    sg.zeroInitialize();
    Real* u = sg.u().data();
    const Index halo = sg.haloWidth();
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const Index gI = offsetX + i - halo;
            const Index gJ = offsetY + j - halo;
            u[sg.index(i, j, sg.kBegin())] = static_cast<Real>(gI) + 1000.0 * static_cast<Real>(gJ);
        }
    }
}

void fill3D(Subgrid& sg, Index offsetX, Index offsetY, Index offsetZ) {
    sg.zeroInitialize();
    Real* u = sg.u().data();
    const Index halo = sg.haloWidth();
    for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                const Index gI = offsetX + i - halo;
                const Index gJ = offsetY + j - halo;
                const Index gK = offsetZ + k - halo;
                u[sg.index(i, j, k)] = static_cast<Real>(gI) +
                                       1000.0 * static_cast<Real>(gJ) +
                                       1000000.0 * static_cast<Real>(gK);
            }
        }
    }
}

bool fileExists(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    return f.is_open();
}

} // namespace

// ---------------------------------------------------------------------------
// S-1 (§6.1/§6.3): np=1 final-solution output, 16x16x1. Header fields and all
// 256 doubles must match the reference bit-for-bit.
// ---------------------------------------------------------------------------
TEST(MpiIoOutputTest, SingleRankHeaderAndDataMatchReference) {
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 1) {
        GTEST_SKIP() << "SingleRankHeaderAndDataMatchReference requires exactly 1 process";
    }

    const Index nx = 16, ny = 16, nz = 1;
    const double dx = 0.125, dy = 0.25, dz = 1.0;
    const std::string base = "./mpiio_test_np1";
    const std::string path = base + ".bin"; // backend appends ".bin" if missing

    Subgrid sg(nx, ny, nz, 1, MPI_COMM_SELF);
    sg.setOffsets(0, 0, 0);
    sg.setBoundaryCondition(BoundaryCondition::Dirichlet);
    fill2D(sg, 0, 0);

    Grid grid;
    grid.nx = nx;
    grid.ny = ny;
    grid.nz = nz;
    grid.dx = dx;
    grid.dy = dy;
    grid.dz = dz;

    MPIIOBinaryBackend backend(grid);
    EXPECT_NO_THROW(backend.write(sg, base));

    std::vector<unsigned char> buf;
    ASSERT_TRUE(readFileBytes(path, buf)) << "output file " << path << " must exist";
    ASSERT_EQ(static_cast<std::size_t>(kHeaderBytes) +
                  static_cast<std::size_t>(nx) * ny * nz * sizeof(double),
              buf.size());

    const FileHeader h = parseHeader(buf);
    expectHeaderMatches(h, nx, ny, nz, dx, dy, dz, 0u /* dirichlet */);

    // Every data value must equal ref2D bit-for-bit (exact integers).
    EXPECT_TRUE(dataAreaMatches(buf, h, &ref2D))
        << "data area must match f(gI,gJ)=gI+1000*gJ for all 256 cells";

    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// S-2 (§6.1/§6.3): np=4 (2x2 cart), 32x32 local -> 64x64 global. Rank 0 reads
// the file back and checks the globally assembled field; the verdict flag is
// broadcast so every rank asserts it.
// ---------------------------------------------------------------------------
TEST(MpiIoOutputTest, FourRanksGlobalAssemblyMatchesReference) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "FourRanksGlobalAssemblyMatchesReference requires exactly 4 processes";
    }

    int dims[2] = {2, 2};
    int periods[2] = {0, 0};
    MPI_Comm cart = MPI_COMM_NULL;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart);

    int coords[2] = {0, 0};
    MPI_Cart_coords(cart, rank, 2, coords);
    const Index offsetX = static_cast<Index>(coords[0]) * 32;
    const Index offsetY = static_cast<Index>(coords[1]) * 32;

    Subgrid sg(32, 32, 1, 1, cart);
    sg.setOffsets(offsetX, offsetY, 0);
    sg.setBoundaryCondition(BoundaryCondition::Dirichlet);
    fill2D(sg, offsetX, offsetY);

    Grid grid;
    grid.nx = 64;
    grid.ny = 64;
    grid.nz = 1;
    grid.dx = 0.5;
    grid.dy = 0.75;
    grid.dz = 1.0;

    MPIIOBinaryBackend backend(grid);
    // Collective contract: every rank writes with the same filename.
    EXPECT_NO_THROW(backend.write(sg, "./mpiio_test_np4"));
    const std::string path = "./mpiio_test_np4.bin";

    MPI_Barrier(cart);

    int flag = 1;
    if (rank == 0) {
        std::vector<unsigned char> buf;
        if (!readFileBytes(path, buf)) {
            ADD_FAILURE() << "output file " << path << " must exist";
            flag = 0;
        } else if (buf.size() != kHeaderBytes +
                                     static_cast<std::size_t>(64) * 64 * 1 * sizeof(double)) {
            ADD_FAILURE() << "unexpected file size " << buf.size();
            flag = 0;
        } else {
            const FileHeader h = parseHeader(buf);
            expectHeaderMatches(h, 64, 64, 1, 0.5, 0.75, 1.0, 0u /* dirichlet */);
            if (!headerMatches(h, 64, 64, 1, 0.5, 0.75, 1.0, 0u)) flag = 0;
            if (!dataAreaMatches(buf, h, &ref2D)) {
                ADD_FAILURE() << "global assembly mismatch (gap/overlap/misplacement)";
                flag = 0;
            }
        }
        std::remove(path.c_str());
    }
    MPI_Bcast(&flag, 1, MPI_INT, 0, cart);
    EXPECT_TRUE(flag == 1);

    MPI_Comm_free(&cart);
}

// ---------------------------------------------------------------------------
// S-3 (§6.1/§6.3): np=4 3D layout, 16x16x8 local -> 32x32x8 global. Header
// dims, data-area size, and sampled points (rank 0 checks all cells; the
// verdict is broadcast).
// ---------------------------------------------------------------------------
TEST(MpiIoOutputTest, FourRanks3DLayout) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 4) {
        GTEST_SKIP() << "FourRanks3DLayout requires exactly 4 processes";
    }

    int dims[3] = {2, 2, 1};
    int periods[3] = {0, 0, 0};
    MPI_Comm cart = MPI_COMM_NULL;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart);

    int coords[3] = {0, 0, 0};
    MPI_Cart_coords(cart, rank, 3, coords);
    const Index offsetX = static_cast<Index>(coords[0]) * 16;
    const Index offsetY = static_cast<Index>(coords[1]) * 16;
    const Index offsetZ = 0;

    Subgrid sg(16, 16, 8, 1, cart);
    sg.setOffsets(offsetX, offsetY, offsetZ);
    // Exercise the other bcType encoding while we are here.
    sg.setBoundaryCondition(BoundaryCondition::Neumann);
    fill3D(sg, offsetX, offsetY, offsetZ);

    Grid grid;
    grid.nx = 32;
    grid.ny = 32;
    grid.nz = 8;
    grid.dx = 0.5;
    grid.dy = 0.125;
    grid.dz = 0.0625;

    MPIIOBinaryBackend backend(grid);
    EXPECT_NO_THROW(backend.write(sg, "./mpiio_test_np4_3d"));
    const std::string path = "./mpiio_test_np4_3d.bin";

    MPI_Barrier(cart);

    int flag = 1;
    if (rank == 0) {
        std::vector<unsigned char> buf;
        if (!readFileBytes(path, buf)) {
            ADD_FAILURE() << "output file " << path << " must exist";
            flag = 0;
        } else if (buf.size() != kHeaderBytes +
                                     static_cast<std::size_t>(32) * 32 * 8 * sizeof(double)) {
            ADD_FAILURE() << "unexpected file size " << buf.size();
            flag = 0;
        } else {
            const FileHeader h = parseHeader(buf);
            expectHeaderMatches(h, 32, 32, 8, 0.5, 0.125, 0.0625, 1u /* neumann */);
            if (!headerMatches(h, 32, 32, 8, 0.5, 0.125, 0.0625, 1u)) flag = 0;
            if (!dataAreaMatches(buf, h, &ref3D)) {
                ADD_FAILURE() << "3D data area mismatch";
                flag = 0;
            }
        }
        std::remove(path.c_str());
    }
    MPI_Bcast(&flag, 1, MPI_INT, 0, cart);
    EXPECT_TRUE(flag == 1);

    MPI_Comm_free(&cart);
}

// ---------------------------------------------------------------------------
// X-1 (§6.4): writing to an unwritable/nonexistent directory must not throw,
// must not kill the process, and a subsequent write to a valid path must
// still succeed. Uses MPI_COMM_SELF subgrids so the test is valid at any np.
// ---------------------------------------------------------------------------
TEST(MpiIoOutputTest, WriteToInvalidPathWarnsAndSkips) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const Index nx = 8, ny = 8, nz = 1;
    Subgrid sg(nx, ny, nz, 1, MPI_COMM_SELF);
    sg.setOffsets(0, 0, 0);
    sg.setBoundaryCondition(BoundaryCondition::Dirichlet);
    fill2D(sg, 0, 0);

    Grid grid;
    grid.nx = nx;
    grid.ny = ny;
    grid.nz = nz;

    MPIIOBinaryBackend backend(grid);
    EXPECT_NO_THROW(backend.write(sg, "/nonexistent_dir_xyz/solution"));

    // Process is alive and a follow-up write to a valid path still works.
    const std::string path = "./mpiio_recovery_r" + std::to_string(rank) + ".bin";
    EXPECT_NO_THROW(backend.write(sg, "./mpiio_recovery_r" + std::to_string(rank)));

    std::vector<unsigned char> buf;
    ASSERT_TRUE(readFileBytes(path, buf)) << "recovery write must produce " << path;
    ASSERT_EQ(static_cast<std::size_t>(kHeaderBytes) +
                  static_cast<std::size_t>(nx) * ny * nz * sizeof(double),
              buf.size());
    const FileHeader h = parseHeader(buf);
    expectHeaderMatches(h, nx, ny, nz, 1.0, 1.0, 1.0, 0u /* dirichlet */);
    EXPECT_TRUE(dataAreaMatches(buf, h, &ref2D));

    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// X-2 (§6.4): a pre-existing longer file with the same name must be truncated
// to exactly 72 + N*8 bytes with no stale tail left behind.
// ---------------------------------------------------------------------------
TEST(MpiIoOutputTest, TruncatesStaleLongerFile) {
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 1) {
        GTEST_SKIP() << "TruncatesStaleLongerFile requires exactly 1 process";
    }

    const Index nx = 12, ny = 12, nz = 1;
    const std::string path = "./mpiio_truncate_test.bin";

    // Pre-create a longer garbage file under the final product name.
    {
        std::ofstream garbage(path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(garbage.is_open());
        const std::size_t garbageSize =
            kHeaderBytes + static_cast<std::size_t>(nx) * ny * nz * sizeof(double) + 999;
        for (std::size_t i = 0; i < garbageSize; ++i) {
            garbage.put(static_cast<char>(0xA5));
        }
        ASSERT_TRUE(garbage.good());
    }

    Subgrid sg(nx, ny, nz, 1, MPI_COMM_SELF);
    sg.setOffsets(0, 0, 0);
    sg.setBoundaryCondition(BoundaryCondition::Dirichlet);
    fill2D(sg, 0, 0);

    Grid grid;
    grid.nx = nx;
    grid.ny = ny;
    grid.nz = nz;

    MPIIOBinaryBackend backend(grid);
    // Pass the name without suffix; the product is "<base>.bin".
    EXPECT_NO_THROW(backend.write(sg, "./mpiio_truncate_test"));

    std::vector<unsigned char> buf;
    ASSERT_TRUE(readFileBytes(path, buf));
    ASSERT_EQ(static_cast<std::size_t>(kHeaderBytes) +
                  static_cast<std::size_t>(nx) * ny * nz * sizeof(double),
              buf.size())
        << "stale longer file must be truncated to 72 + N*8 bytes";
    const FileHeader h = parseHeader(buf);
    expectHeaderMatches(h, nx, ny, nz, 1.0, 1.0, 1.0, 0u /* dirichlet */);
    EXPECT_TRUE(dataAreaMatches(buf, h, &ref2D));

    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// E-1 / E-1b / E-1c (§6.2): suffix handling — missing suffix is appended,
// an existing ".bin" is not doubled, and an empty name does not crash.
// ---------------------------------------------------------------------------
TEST(MpiIoOutputTest, SuffixHandling) {
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 1) {
        GTEST_SKIP() << "SuffixHandling requires exactly 1 process";
    }

    const Index nx = 4, ny = 4, nz = 1;

    Subgrid sg(nx, ny, nz, 1, MPI_COMM_SELF);
    sg.setOffsets(0, 0, 0);
    sg.setBoundaryCondition(BoundaryCondition::Dirichlet);
    fill2D(sg, 0, 0);

    Grid grid;
    grid.nx = nx;
    grid.ny = ny;
    grid.nz = nz;

    MPIIOBinaryBackend backend(grid);

    // E-1: no suffix -> "<base>.bin".
    EXPECT_NO_THROW(backend.write(sg, "./mpiio_suffix_base"));
    EXPECT_TRUE(fileExists("./mpiio_suffix_base.bin"))
        << "E-1: expected product ./mpiio_suffix_base.bin";
    if (fileExists("./mpiio_suffix_base.bin")) {
        std::vector<unsigned char> buf;
        ASSERT_TRUE(readFileBytes("./mpiio_suffix_base.bin", buf));
        ASSERT_EQ(static_cast<std::size_t>(kHeaderBytes) +
                      static_cast<std::size_t>(nx) * ny * nz * sizeof(double),
                  buf.size());
    }
    std::remove("./mpiio_suffix_base.bin");

    // E-1b: already suffixed -> no double suffix.
    EXPECT_NO_THROW(backend.write(sg, "./mpiio_suffix_explicit.bin"));
    EXPECT_TRUE(fileExists("./mpiio_suffix_explicit.bin"))
        << "E-1b: expected product ./mpiio_suffix_explicit.bin";
    EXPECT_FALSE(fileExists("./mpiio_suffix_explicit.bin.bin"))
        << "E-1b: double suffix must not appear";
    std::remove("./mpiio_suffix_explicit.bin");
    std::remove("./mpiio_suffix_explicit.bin.bin");

    // E-1c: empty filename -> smoke only, must not crash; clean up ".bin".
    EXPECT_NO_THROW(backend.write(sg, ""));
    std::remove(".bin");
}

// ---------------------------------------------------------------------------
// E-2 (§6.2): format() reports "mpibin".
// ---------------------------------------------------------------------------
TEST(MpiIoOutputTest, FormatIsMpibin) {
    Grid grid;
    grid.nx = 4;
    grid.ny = 4;
    grid.nz = 1;

    MPIIOBinaryBackend backend(grid);
    EXPECT_EQ("mpibin", backend.format());
}
