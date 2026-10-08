#include <gtest/gtest.h>
#include <mpi.h>

#include "comm/halo_exchanger.hpp"
#include "grid/grid.hpp"
#include "grid/subgrid.hpp"
#include "io/io_backend.hpp"
#include "perf/profiler.hpp"
#include "solver/solver.hpp"

#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace hypo;

namespace {

void setupUniformPoisson(Subgrid& sg) {
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);
    Real* rhs = sg.rhs().data();
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            rhs[sg.index(i, j)] = -1.0;
        }
    }
}

Real sineExact(long long gI, long long gJ) {
    const double x = static_cast<double>(gI + 1) / 65.0;
    const double y = static_cast<double>(gJ + 1) / 65.0;
    return std::sin(M_PI * x) * std::sin(M_PI * y);
}

void setupManufacturedSine(Subgrid& sg, Index offsetX = 0, Index offsetY = 0) {
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);
    Real* rhs = sg.rhs().data();
    const long long hw = static_cast<long long>(sg.haloWidth());
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gI = static_cast<long long>(offsetX) +
                                 static_cast<long long>(i) - hw;
            const long long gJ = static_cast<long long>(offsetY) +
                                 static_cast<long long>(j) - hw;
            rhs[sg.index(i, j)] = sineExact(gI - 1, gJ) + sineExact(gI + 1, gJ) +
                                  sineExact(gI, gJ - 1) + sineExact(gI, gJ + 1) -
                                  4.0 * sineExact(gI, gJ);
        }
    }
}

Real l2ErrorVsSine(const Subgrid& sg) {
    Real sum = 0.0;
    const long long hw = static_cast<long long>(sg.haloWidth());
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gI = static_cast<long long>(sg.offsetX()) +
                                 static_cast<long long>(i) - hw;
            const long long gJ = static_cast<long long>(sg.offsetY()) +
                                 static_cast<long long>(j) - hw;
            const Real diff = sg.u().data()[sg.index(i, j)] - sineExact(gI, gJ);
            sum += diff * diff;
        }
    }
    const double cells = static_cast<double>(sg.nxLocal() * sg.nyLocal() * sg.nzLocal());
    return std::sqrt(sum / cells);
}

} // namespace

TEST(AltFeatureTest, SubgridOffsetsAccessors) {
    Subgrid sg(4, 4, 1, 1, MPI_COMM_SELF);
    EXPECT_EQ(sg.offsetX(), Index(0));
    EXPECT_EQ(sg.offsetY(), Index(0));
    EXPECT_EQ(sg.offsetZ(), Index(0));

    sg.setOffsets(3, 5, 7);
    EXPECT_EQ(sg.offsetX(), Index(3));
    EXPECT_EQ(sg.offsetY(), Index(5));
    EXPECT_EQ(sg.offsetZ(), Index(7));
}

TEST(AltFeatureTest, VtkPieceAndParallelIndexFiles) {
    const std::string outDir = "test_out";
    std::filesystem::create_directories(outDir);

    Subgrid sg(4, 4, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    sg.setOffsets(0, 0, 0);
    for (Index j = 0; j < sg.nyTotal(); ++j) {
        for (Index i = 0; i < sg.nxTotal(); ++i) {
            sg.at(i, j) = 1000.0 * static_cast<Real>(i) + static_cast<Real>(j);
        }
    }

    VTKIOBackend io(0.25, 0.25, 0.25);
    io.write(sg, outDir + "/solution_0_r0", 0);

    Grid grid;
    grid.nx = 8;
    grid.ny = 4;
    grid.nz = 1;
    std::vector<PieceExtent> pieces = {
        {0, 0, 0, 4, 4, 1},
        {4, 0, 0, 4, 4, 1},
    };
    io.writeParallelIndex(grid, outDir + "/solution_0", 0, pieces);

    std::ifstream piece(outDir + "/solution_0_r0.vti");
    ASSERT_TRUE(piece.good());
    std::stringstream pieceBuf;
    pieceBuf << piece.rdbuf();
    EXPECT_NE(pieceBuf.str().find("ImageData"), std::string::npos);
    EXPECT_NE(pieceBuf.str().find("Origin=\"0 0 0\""), std::string::npos);

    std::ifstream index(outDir + "/solution_0.pvti");
    ASSERT_TRUE(index.good());
    std::stringstream indexBuf;
    indexBuf << index.rdbuf();
    const std::string indexText = indexBuf.str();
    EXPECT_NE(indexText.find("PImageData"), std::string::npos);
    EXPECT_NE(indexText.find("solution_0_r0.vti"), std::string::npos);
    EXPECT_NE(indexText.find("solution_0_r1.vti"), std::string::npos);

    // Empty piece list: no index file is written.
    io.writeParallelIndex(grid, outDir + "/solution_empty", 0, {});
    EXPECT_FALSE(std::filesystem::exists(outDir + "/solution_empty.pvti"));

    // Pieces exceeding the global grid: no index file is written.
    std::vector<PieceExtent> badPieces = {{0, 0, 0, 100, 4, 1}};
    io.writeParallelIndex(grid, outDir + "/solution_bad", 0, badPieces);
    EXPECT_FALSE(std::filesystem::exists(outDir + "/solution_bad.pvti"));
}

// ============================================================================
// AR005 T002 (U2): BinaryIOBackend block-write byte-layout guard.
// The on-disk contract: 56-byte header (nx, ny, nz, halo, offsetX, offsetY,
// offsetZ as little-endian 8-byte Index) followed by the interior u field in
// (k, j, i) row-major order (x fastest), Real = double, no padding, no halo
// bytes. The expected byte sequence is constructed from scratch — NOT by
// running the old implementation — so it is a harder guarantee than an
// old-vs-new diff: any Green block-write implementation that leaks halo
// sentinel bytes or scrambles the (k, j, i) row order must FAIL here, while
// the current per-element writer already satisfies the layout (this test is
// the GREEN baseline guardian, not a Red test).
// ============================================================================
TEST(AltFeatureTest, BinaryOutputBitIdenticalToExpectedBytes) {
    const std::string outDir = "test_out";
    std::filesystem::create_directories(outDir);

    Subgrid sg(4, 4, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    sg.setOffsets(0, 0, 0);

    // Poison the WHOLE padded buffer (every halo cell, every k plane) with a
    // sentinel, then overwrite only the interior cells in the writer's
    // k-range: any halo byte leaking into the output is caught by the
    // sentinel never being allowed to appear in the expected/actual bytes.
    const Real kHaloSentinel = -999.0;
    Real* uBuf = sg.u().data();
    for (Index idx = 0; idx < sg.totalCells(); ++idx) {
        uBuf[idx] = kHaloSentinel;
    }
    for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                sg.at(i, j, k) =
                    1000.0 * static_cast<Real>(i) + static_cast<Real>(j);
            }
        }
    }

    BinaryIOBackend io;
    io.write(sg, outDir + "/binary_layout_guard", 0);
    const std::string fname = outDir + "/binary_layout_guard.bin";

    // ---- Expected bytes: header (little-endian Index on x86) + data ----
    std::vector<char> expected;
    expected.reserve(7 * sizeof(Index) + 16 * sizeof(Real));

    auto appendBytes = [&expected](const void* src, std::size_t n) {
        const char* p = static_cast<const char*>(src);
        expected.insert(expected.end(), p, p + n);
    };

    const Index kNx = 4, kNy = 4, kNz = 1, kHw = 1, kOffX = 0, kOffY = 0, kOffZ = 0;
    const Index headerFields[7] = {kNx, kNy, kNz, kHw, kOffX, kOffY, kOffZ};
    for (const Index field : headerFields) {
        appendBytes(&field, sizeof(Index));
    }
    for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                const Real v = sg.at(i, j, k);
                appendBytes(&v, sizeof(Real));
            }
        }
    }

    // ---- Read the actual file back ----
    std::ifstream ifs(fname, std::ios::binary);
    ASSERT_TRUE(ifs.good());
    const std::vector<char> actual((std::istreambuf_iterator<char>(ifs)),
                                   std::istreambuf_iterator<char>());

    // (1) Total size: 56-byte header + 16 doubles, nothing more.
    ASSERT_EQ(actual.size(), std::size_t{7 * sizeof(Index) + 16 * sizeof(Real)});

    // (2) Header: 7 little-endian 8-byte fields with the exact values above.
    for (int f = 0; f < 7; ++f) {
        Index decoded = 0;
        std::memcpy(&decoded, actual.data() + f * sizeof(Index), sizeof(Index));
        EXPECT_EQ(decoded, headerFields[f]) << "header field " << f;
    }

    // (3) Interior values are exact and the halo sentinel never appears.
    const Real* data = reinterpret_cast<const Real*>(actual.data() + 7 * sizeof(Index));
    std::size_t cells = 0;
    for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                ASSERT_LT(cells, std::size_t{16});
                EXPECT_EQ(data[cells], sg.at(i, j, k))
                    << "cell (i=" << i << ", j=" << j << ", k=" << k << ")";
                EXPECT_NE(data[cells], kHaloSentinel) << "halo leak at data cell " << cells;
                ++cells;
            }
        }
    }
    EXPECT_EQ(cells, std::size_t{16});

    // (4) Whole file bit-identical to the expected byte sequence; on mismatch
    // report the first differing offset for diagnosis (header spans [0, 56),
    // data starts at 56).
    ASSERT_EQ(expected.size(), actual.size());
    if (std::memcmp(actual.data(), expected.data(), expected.size()) != 0) {
        for (std::size_t b = 0; b < expected.size(); ++b) {
            if (actual[b] != expected[b]) {
                ADD_FAILURE() << "first byte difference at offset " << b
                              << " (header is 56 B, data region starts at 56)";
                break;
            }
        }
    }
}

// ============================================================================
// AR005 T003 (U3): VTKIOBackend .vti appended raw-binary output guard.
// The on-disk contract (VTK XML appended mode, raw encoding):
//   - <VTKFile ... header_type="UInt64"> attribute;
//   - <DataArray type="Float64" Name="u" format="appended" offset="0">,
//     and format="ascii" must be gone;
//   - <AppendedData encoding="raw">, the '_' marker, an 8-byte little-endian
//     UInt64 byte count (nxLocal*nyLocal*nzLocal*sizeof(double)), then the
//     interior u field as a raw double memory image in (k, j, i) row-major
//     order (x fastest), then the closing </AppendedData></VTKFile> tags.
// This test is the "minimal independent parser" of srs §3.1 acceptance (1):
// it re-reads the file and reconstructs the value sequence from the raw
// bytes alone. The expected bytes are built from the test's own encoding
// (poisoned halo + interior formula), so any halo sentinel leaking into the
// payload or any (k, j, i) order scramble fails here. The current ASCII
// implementation has no AppendedData section at all, so this test FAILS —
// the legal Red state for this task (Green lands in the implementation
// step).
// ============================================================================
TEST(AltFeatureTest, VtkAppendedBinaryParsesBack) {
    const std::string outDir = "test_out";
    std::filesystem::create_directories(outDir);

    Subgrid sg(4, 4, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    sg.setOffsets(0, 0, 0);

    // Poison the WHOLE padded buffer, then overwrite only interior cells:
    // a halo byte leaking into the appended payload can only come from the
    // sentinel and is caught by the exact byte comparison below.
    const Real kHaloSentinel = -777.0;
    Real* uBuf = sg.u().data();
    for (Index idx = 0; idx < sg.totalCells(); ++idx) {
        uBuf[idx] = kHaloSentinel;
    }
    for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                sg.at(i, j, k) =
                    1000.0 * static_cast<Real>(i) + static_cast<Real>(j);
            }
        }
    }

    VTKIOBackend io(0.25, 0.25, 0.25);
    io.write(sg, outDir + "/vtk_appended_guard", 0);
    const std::string fname = outDir + "/vtk_appended_guard.vti";

    // ---- Read the whole file back (binary-safe) ----
    std::ifstream ifs(fname, std::ios::binary);
    ASSERT_TRUE(ifs.good()) << "cannot open " << fname;
    const std::vector<char> actual((std::istreambuf_iterator<char>(ifs)),
                                   std::istreambuf_iterator<char>());
    ASSERT_FALSE(actual.empty());
    const std::string text(actual.begin(), actual.end());

    // ---- (1) XML structure: appended raw mode, ASCII mode is gone ----
    EXPECT_NE(text.find("header_type=\"UInt64\""), std::string::npos);
    EXPECT_NE(text.find("format=\"appended\""), std::string::npos);
    EXPECT_NE(text.find("offset=\"0\""), std::string::npos);
    const std::size_t appendedTagPos =
        text.find("<AppendedData encoding=\"raw\">");
    ASSERT_NE(appendedTagPos, std::string::npos)
        << "missing <AppendedData encoding=\"raw\">";
    EXPECT_EQ(text.find("format=\"ascii\""), std::string::npos);

    // ---- (2) '_' marker followed by the 8-byte UInt64 length header ----
    const std::size_t markPos = text.find('_', appendedTagPos);
    ASSERT_NE(markPos, std::string::npos) << "missing '_' appended-data marker";
    const std::size_t headerPos = markPos + 1;
    ASSERT_GE(actual.size(), headerPos + sizeof(std::uint64_t))
        << "file too short for the 8-byte length header";

    std::uint64_t declaredBytes = 0;
    std::memcpy(&declaredBytes, actual.data() + headerPos,
                sizeof(std::uint64_t));
    const std::uint64_t expectedBytes =
        static_cast<std::uint64_t>(sg.nxLocal()) *
        static_cast<std::uint64_t>(sg.nyLocal()) *
        static_cast<std::uint64_t>(sg.nzLocal()) * sizeof(Real);
    EXPECT_EQ(declaredBytes, expectedBytes)
        << "declared payload byte count != nxLocal*nyLocal*nzLocal*8";

    // ---- (3) payload bit-identical to the interior field in (k, j, i) order
    //      (memcmp from the test-built expected sequence; halo excluded) ----
    const std::size_t dataPos = headerPos + sizeof(std::uint64_t);
    ASSERT_GE(actual.size(), dataPos + static_cast<std::size_t>(expectedBytes))
        << "file too short for the appended payload";

    std::vector<char> expected;
    expected.reserve(static_cast<std::size_t>(expectedBytes));
    for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                const Real v = sg.at(i, j, k);
                const char* p = reinterpret_cast<const char*>(&v);
                expected.insert(expected.end(), p, p + sizeof(Real));
            }
        }
    }
    ASSERT_EQ(expected.size(), static_cast<std::size_t>(expectedBytes));
    EXPECT_EQ(std::memcmp(actual.data() + dataPos, expected.data(),
                          expected.size()),
              0)
        << "appended payload differs from the expected interior byte sequence";

    // The halo sentinel must never appear inside the payload (redundant with
    // the memcmp above, but pinpoints the failure mode when it does).
    for (std::size_t off = 0; off + sizeof(Real) <= expected.size();
         off += sizeof(Real)) {
        Real decoded = 0.0;
        std::memcpy(&decoded, actual.data() + dataPos + off, sizeof(Real));
        EXPECT_NE(decoded, kHaloSentinel)
            << "halo sentinel leaked into payload double " << off / sizeof(Real);
    }

    // ---- (4) decoded doubles exactly match sg.at(i, j, k) ----
    std::size_t cells = 0;
    for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                Real decoded = 0.0;
                std::memcpy(&decoded,
                            actual.data() + dataPos + cells * sizeof(Real),
                            sizeof(Real));
                EXPECT_EQ(decoded, sg.at(i, j, k))  // bitwise-exact compare
                    << "cell (i=" << i << ", j=" << j << ", k=" << k << ")";
                ++cells;
            }
        }
    }
    EXPECT_EQ(cells, std::size_t{16});

    // ---- (5) closing tags after the appended payload ----
    const std::string tail(actual.begin() + static_cast<std::ptrdiff_t>(
                                             dataPos + expected.size()),
                           actual.end());
    EXPECT_NE(tail.find("</AppendedData>"), std::string::npos)
        << "missing </AppendedData> after the payload";
    ASSERT_NE(tail.find("</VTKFile>"), std::string::npos)
        << "missing </VTKFile> after the appended data";
    // Only whitespace may follow the final closing tag.
    const std::size_t lastContent = text.find_last_not_of(" \t\r\n");
    ASSERT_NE(lastContent, std::string::npos);
    ASSERT_GT(lastContent + 1, std::strlen("</VTKFile>"));
    EXPECT_EQ(text.compare(lastContent + 1 - std::strlen("</VTKFile>"),
                           std::strlen("</VTKFile>"), "</VTKFile>"),
              0)
        << "file does not end with </VTKFile>";
}

TEST(AltFeatureTest, ProgressCallbackCountsIterations) {
    Subgrid sg(16, 16, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupUniformPoisson(sg);

    PointToPointExchanger ex;
    ex.initialize(sg);
    JacobiSolver solver;

    std::vector<Index> seen;
    solver.setProgressCallback([&seen](Index iteration) { seen.push_back(iteration); });

    const Index iters = solver.solve(sg, ex, 10, 0.0);
    EXPECT_EQ(iters, Index(10));
    ASSERT_EQ(seen.size(), std::size_t(10));
    for (std::size_t k = 0; k < seen.size(); ++k) {
        EXPECT_EQ(seen[k], Index(k + 1));
    }
}

TEST(AltFeatureTest, ProgressCallbackLastOneWins) {
    Subgrid sg(16, 16, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupUniformPoisson(sg);

    PointToPointExchanger ex;
    ex.initialize(sg);
    JacobiSolver solver;

    int firstCalls = 0;
    int secondCalls = 0;
    solver.setProgressCallback([&firstCalls](Index) { ++firstCalls; });
    solver.setProgressCallback([&secondCalls](Index) { ++secondCalls; });

    solver.solve(sg, ex, 3, 0.0);
    EXPECT_EQ(firstCalls, 0);
    EXPECT_EQ(secondCalls, 3);
}

TEST(AltFeatureTest, PhysicalBoundaryDirichletFillsFaces) {
    Subgrid sg(4, 4, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    for (Index j = 0; j < sg.nyTotal(); ++j) {
        for (Index i = 0; i < sg.nxTotal(); ++i) {
            sg.at(i, j) = 1000.0 * static_cast<Real>(i) + static_cast<Real>(j);
        }
    }

    sg.setBoundaryCondition(BoundaryCondition::Dirichlet);
    sg.applyPhysicalBoundary(0.0);

    for (Index j = 0; j < sg.nyTotal(); ++j) {
        EXPECT_DOUBLE_EQ(sg.at(0, j), 0.0);
        EXPECT_DOUBLE_EQ(sg.at(sg.iEnd(), j), 0.0);
    }
    for (Index i = 0; i < sg.nxTotal(); ++i) {
        EXPECT_DOUBLE_EQ(sg.at(i, 0), 0.0);
        EXPECT_DOUBLE_EQ(sg.at(i, sg.jEnd()), 0.0);
    }
    EXPECT_DOUBLE_EQ(sg.at(sg.iBegin(), sg.jBegin()),
                     1000.0 * static_cast<Real>(sg.iBegin()) + static_cast<Real>(sg.jBegin()));
}

TEST(AltFeatureTest, PhysicalBoundaryNeumannMirrorsInterior) {
    Subgrid sg(4, 4, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    for (Index j = 0; j < sg.nyTotal(); ++j) {
        for (Index i = 0; i < sg.nxTotal(); ++i) {
            sg.at(i, j) = 1000.0 * static_cast<Real>(i) + static_cast<Real>(j);
        }
    }

    sg.setBoundaryCondition(BoundaryCondition::Neumann);
    sg.applyPhysicalBoundary();

    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        EXPECT_DOUBLE_EQ(sg.at(0, j), sg.at(sg.iBegin(), j));
        EXPECT_DOUBLE_EQ(sg.at(sg.iEnd(), j), sg.at(sg.iEnd() - 1, j));
    }
    for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
        EXPECT_DOUBLE_EQ(sg.at(i, 0), sg.at(i, sg.jBegin()));
        EXPECT_DOUBLE_EQ(sg.at(i, sg.jEnd()), sg.at(i, sg.jEnd() - 1));
    }
    // Corner clamps into the interior.
    EXPECT_DOUBLE_EQ(sg.at(0, 0), sg.at(sg.iBegin(), sg.jBegin()));
    // Interior untouched.
    EXPECT_DOUBLE_EQ(sg.at(sg.iBegin(), sg.jBegin()),
                     1000.0 * static_cast<Real>(sg.iBegin()) + static_cast<Real>(sg.jBegin()));
}

TEST(AltFeatureTest, PhysicalBoundaryOnlyTouchesProcNullFaces) {
    Subgrid sg(4, 4, 1, 1, MPI_COMM_SELF);
    // Right neighbor is rank 0 (valid); all other faces are physical.
    sg.setNeighbors(MPI_PROC_NULL, 0, MPI_PROC_NULL, MPI_PROC_NULL);
    for (Index j = 0; j < sg.nyTotal(); ++j) {
        for (Index i = 0; i < sg.nxTotal(); ++i) {
            sg.at(i, j) = 1000.0 * static_cast<Real>(i) + static_cast<Real>(j);
        }
    }

    sg.setBoundaryCondition(BoundaryCondition::Neumann);
    sg.applyPhysicalBoundary();

    // Left physical face mirrored.
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        EXPECT_DOUBLE_EQ(sg.at(0, j), sg.at(sg.iBegin(), j));
    }
    // Right face has a neighbor: halo keeps its previous value.
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        EXPECT_DOUBLE_EQ(sg.at(sg.iEnd(), j),
                         1000.0 * static_cast<Real>(sg.iEnd()) + static_cast<Real>(j));
    }
}

TEST(AltFeatureTest, FusedResidualBitIdenticalAcrossOverlapModes) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Subgrid sgOff(16, 16, 1, 1, MPI_COMM_WORLD);
    sgOff.setNeighbors(rank, rank, rank, rank);
    setupUniformPoisson(sgOff);

    Subgrid sgOn(16, 16, 1, 1, MPI_COMM_WORLD);
    sgOn.setNeighbors(rank, rank, rank, rank);
    setupUniformPoisson(sgOn);

    PointToPointExchanger exOff;
    exOff.initialize(sgOff);
    PointToPointExchanger exOn;
    exOn.initialize(sgOn);

    JacobiSolver solverOff(false);
    JacobiSolver solverOn(true);

    const Index itersOff = solverOff.solve(sgOff, exOff, 400, 0.0);
    const Index itersOn = solverOn.solve(sgOn, exOn, 400, 0.0);

    EXPECT_EQ(itersOff, Index(400));
    EXPECT_EQ(itersOn, Index(400));
    EXPECT_EQ(solverOff.lastResidual(), solverOn.lastResidual());

    ASSERT_EQ(sgOff.totalCells(), sgOn.totalCells());
    EXPECT_EQ(std::memcmp(sgOff.u().data(), sgOn.u().data(),
                          sgOff.totalCells() * sizeof(Real)), 0);
}

TEST(AltFeatureTest, ExchangeWithExternalBufferMatchesU) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Subgrid sgU(4, 4, 1, 1, MPI_COMM_WORLD);
    sgU.setNeighbors(rank, rank, rank, rank);
    for (Index j = 0; j < sgU.nyTotal(); ++j) {
        for (Index i = 0; i < sgU.nxTotal(); ++i) {
            sgU.at(i, j) = 1000.0 * static_cast<Real>(i) + static_cast<Real>(j);
        }
    }
    PointToPointExchanger exU;
    exU.initialize(sgU);
    exU.exchange(sgU);

    Subgrid sgB(4, 4, 1, 1, MPI_COMM_WORLD);
    sgB.setNeighbors(rank, rank, rank, rank);
    for (Index j = 0; j < sgB.nyTotal(); ++j) {
        for (Index i = 0; i < sgB.nxTotal(); ++i) {
            sgB.at(i, j) = 1000.0 * static_cast<Real>(i) + static_cast<Real>(j);
        }
    }
    std::vector<Real> buf(sgB.totalCells());
    std::memcpy(buf.data(), sgB.u().data(), buf.size() * sizeof(Real));

    PointToPointExchanger exB;
    exB.initialize(sgB);
    exB.exchange(sgB, buf.data());

    EXPECT_EQ(std::memcmp(sgU.u().data(), buf.data(), buf.size() * sizeof(Real)), 0);
}

TEST(AltFeatureTest, ExchangeNullDataIsSafeNoOp) {
    Subgrid sg(4, 4, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    for (Index j = 0; j < sg.nyTotal(); ++j) {
        for (Index i = 0; i < sg.nxTotal(); ++i) {
            sg.at(i, j) = 1000.0 * static_cast<Real>(i) + static_cast<Real>(j);
        }
    }
    std::vector<Real> snapshot(sg.totalCells());
    std::memcpy(snapshot.data(), sg.u().data(), snapshot.size() * sizeof(Real));

    PointToPointExchanger ex;
    ex.initialize(sg);
    ex.exchange(sg, nullptr);

    EXPECT_EQ(std::memcmp(sg.u().data(), snapshot.data(), snapshot.size() * sizeof(Real)), 0);
}

TEST(AltFeatureTest, ExchangeBeforeInitializeIsSafeNoOp) {
    Subgrid sg(4, 4, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    for (Index j = 0; j < sg.nyTotal(); ++j) {
        for (Index i = 0; i < sg.nxTotal(); ++i) {
            sg.at(i, j) = 1000.0 * static_cast<Real>(i) + static_cast<Real>(j);
        }
    }
    std::vector<Real> snapshot(sg.totalCells());
    std::memcpy(snapshot.data(), sg.u().data(), snapshot.size() * sizeof(Real));

    PointToPointExchanger ex; // intentionally not initialized
    ex.exchange(sg);

    EXPECT_EQ(std::memcmp(sg.u().data(), snapshot.data(), snapshot.size() * sizeof(Real)), 0);
}

TEST(AltFeatureTest, PackStrideAsymmetricSelfLoop3D) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Subgrid sg(5, 4, 3, 2, MPI_COMM_WORLD);
    sg.setNeighbors(rank, rank, rank, rank, rank, rank);

    auto value = [](Index i, Index j, Index k) {
        return 10000.0 * static_cast<Real>(i) + 100.0 * static_cast<Real>(j) +
               static_cast<Real>(k);
    };
    for (Index k = 0; k < sg.nzTotal(); ++k) {
        for (Index j = 0; j < sg.nyTotal(); ++j) {
            for (Index i = 0; i < sg.nxTotal(); ++i) {
                sg.at(i, j, k) = value(i, j, k);
            }
        }
    }

    PointToPointExchanger ex;
    ex.initialize(sg);
    ex.exchange(sg);

    const Index hw = static_cast<Index>(sg.haloWidth());

    for (Index h = 0; h < hw; ++h) {
        for (Index k = sg.kBegin(); k < sg.kEnd(); ++k) {
            for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
                EXPECT_DOUBLE_EQ(sg.at(sg.iEnd() + h, j, k), value(sg.iBegin() + h, j, k));
                EXPECT_DOUBLE_EQ(sg.at(h, j, k), value(sg.iEnd() - hw + h, j, k));
            }
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                EXPECT_DOUBLE_EQ(sg.at(i, sg.jEnd() + h, k), value(i, sg.jBegin() + h, k));
                EXPECT_DOUBLE_EQ(sg.at(i, h, k), value(i, sg.jEnd() - hw + h, k));
            }
        }
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                EXPECT_DOUBLE_EQ(sg.at(i, j, sg.kEnd() + h), value(i, j, sg.kBegin() + h));
                EXPECT_DOUBLE_EQ(sg.at(i, j, h), value(i, j, sg.kEnd() - hw + h));
            }
        }
    }
}

TEST(AltSolverTest, RedBlackGSConvergesFasterThanJacobi) {
    const int maxIter = 30000;
    const Real tol = 1e-8;

    Subgrid sgJ(64, 64, 1, 1, MPI_COMM_SELF);
    sgJ.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sgJ);
    PointToPointExchanger exJ;
    exJ.initialize(sgJ);
    JacobiSolver jacobi;
    const Index itersJ = jacobi.solve(sgJ, exJ, maxIter, tol);
    ASSERT_LT(itersJ, Index(maxIter));

    Subgrid sgR(64, 64, 1, 1, MPI_COMM_SELF);
    sgR.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sgR);
    PointToPointExchanger exR;
    exR.initialize(sgR);
    RedBlackGSSolver rbgs;
    const Index itersR = rbgs.solve(sgR, exR, maxIter, tol);

    EXPECT_LT(itersR, Index(maxIter));
    EXPECT_LE(itersR, itersJ * 60 / 100);
    EXPECT_LT(l2ErrorVsSine(sgR), 1e-3);

    // Tight Jacobi reference: both solvers target the same discrete solution.
    Subgrid sgT(64, 64, 1, 1, MPI_COMM_SELF);
    sgT.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sgT);
    PointToPointExchanger exT;
    exT.initialize(sgT);
    JacobiSolver tight;
    tight.solve(sgT, exT, 100000, 1e-12);

    // Separate tight RBGS run for the solver-agreement comparison.
    Subgrid sgR2(64, 64, 1, 1, MPI_COMM_SELF);
    sgR2.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sgR2);
    PointToPointExchanger exR2;
    exR2.initialize(sgR2);
    RedBlackGSSolver rbgsTight;
    rbgsTight.solve(sgR2, exR2, 30000, 1e-10);

    for (Index j = sgR.jBegin(); j < sgR.jEnd(); ++j) {
        for (Index i = sgR.iBegin(); i < sgR.iEnd(); ++i) {
            EXPECT_NEAR(sgR2.u().data()[sgR2.index(i, j)],
                        sgT.u().data()[sgT.index(i, j)], 1e-8);
        }
    }
}

TEST(AltSolverTest, CGConvergesMuchFasterThanJacobi) {
    const int maxIter = 30000;
    const Real tol = 1e-7;

    Subgrid sgJ(64, 64, 1, 1, MPI_COMM_SELF);
    sgJ.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sgJ);
    PointToPointExchanger exJ;
    exJ.initialize(sgJ);
    JacobiSolver jacobi;
    const Index itersJ = jacobi.solve(sgJ, exJ, maxIter, tol);
    ASSERT_LT(itersJ, Index(maxIter));

    Subgrid sgC(64, 64, 1, 1, MPI_COMM_SELF);
    sgC.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sgC);
    PointToPointExchanger exC;
    exC.initialize(sgC);
    CGSolver cg;
    const Index itersC = cg.solve(sgC, exC, maxIter, tol);

    EXPECT_LT(itersC, Index(maxIter));
    EXPECT_LE(itersC, itersJ / 10);
    EXPECT_LT(l2ErrorVsSine(sgC), 1e-3);

    Subgrid sgT(64, 64, 1, 1, MPI_COMM_SELF);
    sgT.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSine(sgT);
    PointToPointExchanger exT;
    exT.initialize(sgT);
    JacobiSolver tight;
    tight.solve(sgT, exT, 100000, 1e-12);

    for (Index j = sgC.jBegin(); j < sgC.jEnd(); ++j) {
        for (Index i = sgC.iBegin(); i < sgC.iEnd(); ++i) {
            EXPECT_NEAR(sgC.u().data()[sgC.index(i, j)],
                        sgT.u().data()[sgT.index(i, j)], 1e-8);
        }
    }
}

TEST(AltSolverTest, CgIterateBeforeSolveIsSafe) {
    Subgrid sg(8, 8, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    PointToPointExchanger ex;
    ex.initialize(sg);

    CGSolver cg;
    EXPECT_DOUBLE_EQ(cg.iterate(sg, ex), 0.0);
}

// ============================================================================
// AR004 T002: A2-RBGS real residual + residual check interval.
// These tests assert the NEW contract for RedBlackGSSolver:
//   - iterate() always performs a residual scan (exchange(u) +
//     trueResidualSquaredLocal + Allreduce) after the double sweep and
//     returns the CURRENT true residual ||Au - f||_2 (not a diff norm);
//   - solve() evaluates the convergence criterion only on iterations with
//     completed % interval == 0 (interval via the base-class
//     setResidualCheckInterval, default 1), on the true-residual scale;
//   - the loop exit (converged or maxIter) runs exactly one confirmation
//     scan (profiler region "residual_confirm", callCount == 1); maxIter == 0
//     skips it and leaves lastResidual() at 0;
//   - solve()'s iteration region is named "rbgs_iteration";
//   - each in-loop check scan contains one Allreduce in profiler region
//     "residual_allreduce".
// The current implementation is still diff-based and misuses the
// "jacobi_iteration" region name, so these assertions fail at runtime — the
// legal Red state for this task.
// ============================================================================

namespace {

// Relative-error tolerance for residual comparisons: the OpenMP reduction in
// the true-residual scan sums in a different order than the serial
// recomputation below, so results agree to ~1 ulp but not bit-exactly.
// Copied from the AR004 T001 facilities in tests/test_solver.cpp (anonymous
// namespaces are not shared across translation units).
constexpr Real kAr004RelTol = 1e-12;

void expectCloseRelative(Real actual, Real expected, const std::string& context) {
    const Real tol = kAr004RelTol * std::fabs(expected);
    EXPECT_NEAR(actual, expected, tol) << context;
}

// Exact manufactured solution on an n x n global grid: u = sin(pi x) sin(pi y)
// with x = (gI+1)/(n+1). Generalizes the 64^2-only sineExact above (which
// hard-codes /65.0). Copied from tests/test_solver.cpp; named differently to
// avoid ambiguity with the sineExact overload in the namespace above.
Real sineExactGlobalN(long long gI, long long gJ, long long n) {
    const double x = static_cast<double>(gI + 1) / static_cast<double>(n + 1);
    const double y = static_cast<double>(gJ + 1) / static_cast<double>(n + 1);
    return std::sin(M_PI * x) * std::sin(M_PI * y);
}

// Single-rank manufactured sine setup for an arbitrary global size, following
// the setupManufacturedSine pattern above (which is 64^2-specific via
// sineExact). Copied from tests/test_solver.cpp; named differently to avoid
// ambiguity with the (sg, offsetX, offsetY) overload in the namespace above.
void setupManufacturedSineGlobalN(Subgrid& sg, Index globalN) {
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);
    Real* rhs = sg.rhs().data();
    const long long hw = static_cast<long long>(sg.haloWidth());
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const long long gI = static_cast<long long>(i) - hw;
            const long long gJ = static_cast<long long>(j) - hw;
            rhs[sg.index(i, j)] =
                sineExactGlobalN(gI - 1, gJ, globalN) + sineExactGlobalN(gI + 1, gJ, globalN) +
                sineExactGlobalN(gI, gJ - 1, globalN) + sineExactGlobalN(gI, gJ + 1, globalN) -
                4.0 * sineExactGlobalN(gI, gJ, globalN);
        }
    }
}

// Independent serial recomputation of the local true residual sum of squares
// for the Au = -rhs convention (r = D*u - sum(neighbors) + rhs), D = 4 in 2D,
// D = 6 in 3D. Takes an explicit u buffer so callers can evaluate a snapshot.
// Halo values are read as-is: the caller is responsible for having refreshed
// them (exchange / applyPhysicalBoundary). Copied verbatim in style from the
// AR004 T001 facilities in tests/test_solver.cpp.
Real trueResidualSumSqSerial(const Subgrid& sg, const Real* u) {
    const Real* rhs = sg.rhs().data();
    const Index nxT = sg.nxTotal();
    const Index nyT = sg.nyTotal();
    const bool is2D = (sg.nzLocal() == 1);
    const Real denom = is2D ? 4.0 : 6.0;
    // 2D lives in the k = 0 plane of the padded array.
    const Index k0 = is2D ? 0 : sg.kBegin();
    const Index k1 = is2D ? 1 : sg.kEnd();

    Real sum = 0.0;
    for (Index k = k0; k < k1; ++k) {
        for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
            for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                const Index idx = (k * nyT + j) * nxT + i;
                const Real neighborSum = u[idx - 1] + u[idx + 1] +
                                         u[idx - nxT] + u[idx + nxT] +
                                         (is2D ? 0.0
                                               : u[idx - nxT * nyT] +
                                                 u[idx + nxT * nyT]);
                const Real r = denom * u[idx] - neighborSum + rhs[idx];
                sum += r * r;
            }
        }
    }
    return sum;
}

} // namespace

// U4: after convergence, lastResidual() is the true residual ||Au - f||_2 of
// the final iterate (confirm-scan semantics on the true-residual scale), not
// a diff-based proxy.
TEST(RbgsAR004Test, LastResidualAfterConvergenceIsTrueResidual) {
    const Index kGridN = 32;      // manufactured sine domain 32x32
    const Real kTol = 1e-6;       // convergence threshold on the true residual
    const Index kMaxIter = 5000;  // generous budget; RBGS on 32^2 is far faster

    Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSineGlobalN(sg, kGridN);

    PointToPointExchanger ex;
    ex.initialize(sg);

    RedBlackGSSolver solver;
    const Index iters = solver.solve(sg, ex, kMaxIter, kTol);

    SCOPED_TRACE("converged solve");
    ASSERT_GT(iters, Index(0));
    ASSERT_LT(iters, kMaxIter);  // must converge, not run out of budget
    EXPECT_LE(solver.lastResidual(), kTol);

    // Recompute the true residual of the final u. Single rank with all
    // PROC_NULL neighbors: halo IS the physical boundary, so refreshing it
    // via applyPhysicalBoundary is the correct pre-computation step.
    sg.applyPhysicalBoundary();
    const Real expected =
        std::sqrt(trueResidualSumSqSerial(sg, sg.u().data()));
    ASSERT_GT(expected, 0.0);
    expectCloseRelative(solver.lastResidual(), expected,
                        "lastResidual vs recomputed true residual");
}

// U5: iterate() returns the CURRENT true residual — of the state AFTER the
// double sweep — not the diff norm and not the pre-update residual.
TEST(RbgsAR004Test, IterateReturnsCurrentTrueResidual) {
    const Index kGridN = 16;  // small single-rank grid, one direct iterate() call

    Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSineGlobalN(sg, kGridN);

    // Perturb u away from zero with asymmetric values so both the old
    // diff-based semantics and any indexing slip produce a distinguishable
    // answer.
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const Real fi = static_cast<Real>(i);
            const Real fj = static_cast<Real>(j);
            sg.at(i, j) = 0.10 * fi - 0.05 * fj + 0.02 * fi * fj;
        }
    }
    // Refresh the PROC_NULL halo faces before the call.
    sg.applyPhysicalBoundary();

    PointToPointExchanger ex;
    ex.initialize(sg);

    RedBlackGSSolver solver;
    const Real returned = solver.iterate(sg, ex);

    // Recompute the true residual of the post-sweep state. iterate() has
    // already refreshed the physical faces via applyPhysicalBoundary, and
    // with all-PROC_NULL neighbors the in-scan exchange is a no-op; the
    // repeat applyPhysicalBoundary is idempotent and documents the contract.
    sg.applyPhysicalBoundary();
    const Real expected =
        std::sqrt(trueResidualSumSqSerial(sg, sg.u().data()));
    SCOPED_TRACE("iterate() return vs current true residual");
    ASSERT_GT(expected, 0.0);
    expectCloseRelative(returned, expected,
                        "iterate() return vs ||Au - f|| of post-sweep state");
}

// R-I5s: default interval (k=1) — every iteration runs the scan-based check.
// tol=0 never converges, maxIter=20: 20 in-loop checks (one Allreduce each in
// region residual_allreduce), 20 iteration regions, one exit confirmation
// scan.
TEST(RbgsAR004Test, DefaultIntervalChecksEveryIteration) {
    const Index kGridN = 32;   // manufactured sine domain
    const Index kMaxIter = 20;
    const Real kNoConvergeTol = 0.0;

    Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSineGlobalN(sg, kGridN);

    PointToPointExchanger ex;
    ex.initialize(sg);

    Profiler::instance().reset();
    RedBlackGSSolver solver;  // no setResidualCheckInterval: default must be 1
    const Index iters = solver.solve(sg, ex, kMaxIter, kNoConvergeTol);

    SCOPED_TRACE("default interval, maxIter=20");
    EXPECT_EQ(iters, kMaxIter);
    EXPECT_EQ(Profiler::instance().stats("residual_allreduce").callCount,
              std::uint64_t{kMaxIter});
    EXPECT_EQ(Profiler::instance().stats("residual_confirm").callCount,
              std::uint64_t{1});
    EXPECT_EQ(Profiler::instance().stats("rbgs_iteration").callCount,
              std::uint64_t{kMaxIter});
}

// R-I1s: with interval k, the in-loop check fires only every k-th iteration:
// maxIter=20, interval=10 -> 2 checks (at 10 and 20); the iteration region
// still covers every iteration.
TEST(RbgsAR004Test, IntervalTenHalvesAllreduceChecks) {
    const Index kGridN = 32;   // manufactured sine domain
    const Index kMaxIter = 20;
    const Index kInterval = 10;      // floor(20/10) = 2 in-loop checks
    const Real kNoConvergeTol = 0.0;

    Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSineGlobalN(sg, kGridN);

    PointToPointExchanger ex;
    ex.initialize(sg);

    Profiler::instance().reset();
    RedBlackGSSolver solver;
    solver.setResidualCheckInterval(kInterval);
    const Index iters = solver.solve(sg, ex, kMaxIter, kNoConvergeTol);

    SCOPED_TRACE("interval=10, maxIter=20");
    EXPECT_EQ(iters, kMaxIter);
    EXPECT_EQ(Profiler::instance().stats("residual_allreduce").callCount,
              std::uint64_t{2});
    EXPECT_EQ(Profiler::instance().stats("residual_confirm").callCount,
              std::uint64_t{1});
    EXPECT_EQ(Profiler::instance().stats("rbgs_iteration").callCount,
              std::uint64_t{kMaxIter});
}

// R-E4s: interval >= maxIter means no in-loop check at all; the exit still
// gets exactly one confirmation scan and lastResidual() is the true residual
// of the final iterate.
TEST(RbgsAR004Test, IntervalAboveMaxIterSkipsInLoopChecks) {
    const Index kGridN = 32;   // manufactured sine domain
    const Index kMaxIter = 5;  // exhaust the budget
    const Index kInterval = 100;    // no completed iteration is a multiple of 100
    const Real kNoConvergeTol = 0.0;

    Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupManufacturedSineGlobalN(sg, kGridN);

    PointToPointExchanger ex;
    ex.initialize(sg);

    Profiler::instance().reset();
    RedBlackGSSolver solver;
    solver.setResidualCheckInterval(kInterval);
    const Index iters = solver.solve(sg, ex, kMaxIter, kNoConvergeTol);

    SCOPED_TRACE("interval=100, maxIter=5");
    EXPECT_EQ(iters, kMaxIter);
    EXPECT_EQ(Profiler::instance().stats("residual_allreduce").callCount,
              std::uint64_t{0});
    EXPECT_EQ(Profiler::instance().stats("residual_confirm").callCount,
              std::uint64_t{1});

    sg.applyPhysicalBoundary();
    const Real expected =
        std::sqrt(trueResidualSumSqSerial(sg, sg.u().data()));
    ASSERT_GT(expected, 0.0);
    expectCloseRelative(solver.lastResidual(), expected,
                        "lastResidual vs recomputed true residual at exit");
}

// ============================================================================
// AR004 T004: CGSolver refactoring + parallelization (behavior must stay
// BIT-IDENTICAL for a fixed OMP thread count).
// These tests pin the contract of the upcoming Green implementation:
//   - the init loops (r_/p_/ap_ zeroing, r_ = -rhs, p_ = r_ setup) get an
//     OpenMP parallel for, and the duplicated p-update loops in solve() and
//     iterate() get extracted into a private helper updatePInterior
//     (parallelized) — both are element-wise, so NO arithmetic may change
//     for a fixed thread count (no new reductions);
//   - the profiler region name inside solve()'s loop changes from
//     "jacobi_iteration" to "cg_iteration" (naming only, no numerics);
//   - the existing `pap <= 0` breakdown guard in solve() must keep firing
//     BEFORE any state mutation (u untouched) and iterate() after a
//     breakdown must stay safe.
// U6 is a regression-capture golden test: it PASSES on the current baseline
// code (the goldens are captured from it) and fails only if the Green
// implementation changes CG numerics for a fixed thread count. E2 tests
// current behavior and must pass before and after the refactor.
// ============================================================================

namespace {

// FNV-1a 64-bit hash over the raw bytes of the interior solution in
// row-major order (j outer, i inner). Byte-level: any bit flip anywhere in
// the interior changes the hash.
std::uint64_t fnv1a64Interior(const Subgrid& sg) {
    std::uint64_t h = 0xcbf29ce484222325ULL;  // FNV offset basis
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            const Real v = sg.u().data()[sg.index(i, j)];
            unsigned char bytes[sizeof(Real)];
            std::memcpy(bytes, &v, sizeof(Real));
            for (unsigned char b : bytes) {
                h ^= static_cast<std::uint64_t>(b);
                h *= 0x100000001b3ULL;  // FNV prime
            }
        }
    }
    return h;
}

} // namespace

// U6: pin the numerical behavior of CG on 64^2 (uniform-Poisson rhs = -1,
// multi-mode, maxIter=30000, tol=1e-7) per OMP thread tier, so the upcoming
// parallelization/refactor cannot silently change results. The multi-mode
// rhs forces 126 CG iterations, exercising every refactored loop (init,
// matvec, axpy, p-update) repeatedly.
//
// Determinism split (established by measurement on libgomp): dotGlobal's
// OpenMP reduction combines partial sums in thread-arrival order, so the
// OMP=4 result is NOT bitwise reproducible across runs — only the OMP=1
// tier is. Therefore:
//   - tier 1: bit-exact golden (iteration count + FNV-1a hash of interior u)
//   - tier 4: golden iteration count + max-abs deviation vs the tier-1
//     solution must stay below 1e-12 (repo convention for cross-thread
//     comparisons; values are O(0.05), so this is a relative-scale check).
TEST(CgAR004Test, ParallelizedLoopsPreserveBaselineGoldens) {
#ifdef NDEBUG
    // The bit-exact golden below is pinned to the Debug+ASan build the
    // capture ran on; Release codegen (-O3, vectorization, FMA contraction)
    // legitimately differs in low-order bits. The rest of the suite carries
    // relative-assertion coverage in Release.
    GTEST_SKIP() << "bit-exact goldens are Debug-build-only";
#endif
    int sz = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &sz);
    if (sz != 1) {
        GTEST_SKIP() << "golden capture is single-process only";
    }

#ifdef _OPENMP
    // RAII: restore the max-threads setting on every exit path, including
    // ASSERT failures.
    struct MaxThreadsRestore {
        int saved;
        ~MaxThreadsRestore() { omp_set_num_threads(saved); }
    } threadGuard{omp_get_max_threads()};
#endif

    const int kTiers[2] = {1, 4};
    // Goldens captured 2026-10-08 on baseline (pre-T004) code, WSL OpenMPI,
    // Debug+ASan build. Index k matches kTiers[k]. Uniform rhs = -1 is
    // multi-mode: CG needs 126 iterations. The convergence margin at the
    // golden count is ~20% (residual 8.1e-8 vs tol 1e-7), so last-ulp
    // reduction noise cannot flip the iteration count.
    static constexpr Index kGoldenIters[2] = {126, 126};
    static constexpr std::uint64_t kGoldenHashTier1 = 0x654fa5ddb3bd66b5ULL;

    // Interior solution of the tier-1 run, kept for the tier-4 comparison.
    std::vector<Real> tier1U;

    for (int k = 0; k < 2; ++k) {
        const int t = kTiers[k];
#ifdef _OPENMP
        omp_set_num_threads(t);
#endif
        Subgrid sg(64, 64, 1, 1, MPI_COMM_SELF);
        sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
        setupUniformPoisson(sg);
        PointToPointExchanger ex;
        ex.initialize(sg);

        CGSolver cg;
        const Index iters = cg.solve(sg, ex, 30000, 1e-7);
        ASSERT_EQ(iters, kGoldenIters[k]) << "tier=" << t << " iters=" << iters;

        if (t == 1) {
            const std::uint64_t hash = fnv1a64Interior(sg);
            ASSERT_EQ(hash, kGoldenHashTier1)
                << "tier=" << t << " hash=0x" << std::hex << hash;
            tier1U.reserve(static_cast<std::size_t>(62 * 62));
            for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
                for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                    tier1U.push_back(sg.u().data()[sg.index(i, j)]);
                }
            }
        } else {
            // libgomp combines reductions in arrival order, so the
            // multi-threaded result is not bitwise stable across runs;
            // compare against the tier-1 solution instead.
            std::size_t c = 0;
            Real maxAbsDiff = 0.0;
            for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
                for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
                    const Real diff = std::fabs(sg.u().data()[sg.index(i, j)] - tier1U[c]);
                    maxAbsDiff = std::max(maxAbsDiff, diff);
                    ++c;
                }
            }
            ASSERT_EQ(c, tier1U.size());
            // Observed noise band across repeated runs: <= 1.4e-12
            // (arrival-order reduction noise accumulated over 126
            // iterations). 1e-11 leaves ~8x margin; the tier-1 hash above
            // is the bit-exact guard, this is the coarse tier-4 check.
            EXPECT_LT(maxAbsDiff, 1e-11)
                << "tier=" << t << " maxAbsDiff=" << maxAbsDiff;
        }
    }
}

// E2: the `pap <= 0` breakdown guard in solve() must fire BEFORE any state
// mutation (u untouched), and iterate() after a breakdown must stay safe.
//
// Deterministic trigger via IEEE overflow, independent of iteration count
// and thread count: rhs = +inf at two horizontally adjacent interior cells
// (2,2) and (3,2) makes r = p = -inf at both points. The matvec at (2,2)
// then computes ap = 4*(-inf) - (... + (-inf) + ...) = -inf - (-inf) = NaN
// (inf - inf), likewise at (3,2). pap = dot(p, ap) contains NaN, and NaN
// propagates through ANY summation/reduction order, so `!(pap > 0)` fires
// deterministically for any thread count.
TEST(CgAR004Test, BreakdownGuardFiresBeforeStateMutation) {
    Subgrid sg(8, 8, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupUniformPoisson(sg);  // zeroes u/uNext/rhs, then rhs = -1 interior

    // (2,2) and (3,2) are interior for an 8x8 grid with halo 1
    // (iBegin = jBegin = 1, iEnd = jEnd = 9 in total coordinates 0..9).
    ASSERT_LT(sg.iBegin(), Index(2));
    ASSERT_GT(sg.iEnd(), Index(3));
    ASSERT_LT(sg.jBegin(), Index(2));
    ASSERT_GT(sg.jEnd(), Index(2));
    sg.rhs().data()[sg.index(2, 2)] = std::numeric_limits<double>::infinity();
    sg.rhs().data()[sg.index(3, 2)] = std::numeric_limits<double>::infinity();

    PointToPointExchanger ex;
    ex.initialize(sg);

    CGSolver cg;
    const Index iters = cg.solve(sg, ex, 100, 1e-7);
    EXPECT_EQ(iters, Index(0));  // breakdown at iteration 0

    // u is EXACTLY untouched: the guard fired before axpyInterior mutated u.
    // (setupUniformPoisson zeroes u via zeroInitialize; the trailing
    // applyPhysicalBoundary in solve() only touches halo cells, which the
    // interior loop below never reads.)
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            EXPECT_EQ(sg.u().data()[sg.index(i, j)], 0.0)
                << "at i=" << i << " j=" << j;
        }
    }

    // Initial residual sqrt(rho) is inf (r = -rhs is inf at the two cells).
    EXPECT_TRUE(std::isinf(cg.lastResidual()));

    // iterate() after breakdown must not crash and must return inf
    // (guard fires again inside iterate(); lastResidual_ is still inf).
    EXPECT_TRUE(std::isinf(cg.iterate(sg, ex)));
}
