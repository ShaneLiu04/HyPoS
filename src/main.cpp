#include "core/types.hpp"
#include "core/exception.hpp"
#include "grid/grid.hpp"
#include "grid/subgrid.hpp"
#include "grid/partition.hpp"
#include "solver/solver.hpp"
#include "comm/halo_exchanger.hpp"
#include "perf/timer.hpp"
#include "perf/profiler.hpp"
#include "perf/reporter.hpp"
#include "io/io_backend.hpp"
#include "utils/mpi_env.hpp"
#include "utils/logger.hpp"
#include "utils/cmdline_parser.hpp"
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <system_error>
#include <vector>

// For M_PI on Windows/MSVC
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hypo {

void printUsage(const std::string& programName) {
    std::cout << "Usage: " << programName << " [OPTIONS]\n\n"
              << "Grid Options:\n"
              << "  --nx, --ny, --nz <int>       Global grid dimensions (default: 1024 1024 1)\n"
              << "  --halo-width <int>           Halo layer width (default: 1)\n\n"
              << "Solver Options:\n"
              << "  --solver <string>            Solver type: jacobi, red_black_gs, cg (default: jacobi)\n"
              << "  --max-iter <int>             Maximum iterations (default: 10000)\n"
              << "  --tol <double>               Convergence tolerance (default: 1e-6)\n"
              << "  --bc <string>                Physical BC: dirichlet, neumann (default: dirichlet)\n\n"
              << "Parallel Options:\n"
              << "  --omp-threads <int>          OpenMP threads per process (default: all cores)\n"
              << "  --comm-mode <string>         Communication mode: p2p, collective (default: p2p)\n\n"
              << "Performance Options:\n"
              << "  --enable-profiling           Enable detailed performance profiling\n"
              << "  --overlap-comm               Enable communication-computation overlap\n\n"
              << "I/O Options:\n"
              << "  --output-format <string>     Output format: json, csv, vtk, binary (default: json)\n"
              << "  --output-dir <path>          Output directory (default: ./output)\n"
              << "  --save-interval <int>        Save intermediate results every N steps (0=none)\n\n"
              << "Other:\n"
              << "  --help, -h                   Show this help message\n";
}

void setupProblem(Subgrid& subgrid, const Grid& globalGrid) {
    // Set up a test problem: Poisson equation with known analytical solution
    // For testing, use f(x,y) = 2*pi^2*sin(pi*x)*sin(pi*y) with solution u = sin(pi*x)*sin(pi*y)
    // Discretized on [0,1]x[0,1]

    const Real dx = globalGrid.dx;
    const Real dy = globalGrid.dy;

    Real* u = subgrid.u().data();
    Real* rhs = subgrid.rhs().data();

    // For 2D: need global coordinates of each local cell
    // For simplicity, we compute the RHS directly using the analytical f
    // and set boundary conditions (Dirichlet = 0)

    for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
        Real y = j * dy; // global y coordinate (approximate, assumes uniform grid starting at 0)
        for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
            Real x = i * dx;
            Index idx = subgrid.index(i, j);
            // f = -2*pi^2*sin(pi*x)*sin(pi*y) for Poisson with u = sin(pi*x)*sin(pi*y)
            rhs[idx] = -2.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
            u[idx] = 0.0; // initial guess
        }
    }

    // Apply Dirichlet BC on halos (value = 0)
    subgrid.applyDirichletBC(0.0);
}

} // namespace hypo

int main(int argc, char* argv[]) {
    using namespace hypo;

    // MPI initialization with thread support for OpenMP (RAII managed)
    MPIEnv mpiEnv(argc, argv, MPI_THREAD_FUNNELED);
    const int rank = mpiEnv.rank();
    const int size = mpiEnv.size();

    // Setup logger
    Logger::instance().setRank(rank);
    Logger::instance().setLogLevel(LogLevel::INFO);

    // Parse command line
    CommandLineParser parser;
    parser.parse(argc, argv);

    if (parser.has("help") || parser.has("h")) {
        if (rank == 0) {
            printUsage(argc > 0 ? argv[0] : "hypos");
        }
        return 0;
    }

    // Grid configuration
    Grid grid;
    grid.nx = static_cast<Index>(parser.get<int>("nx", 1024));
    grid.ny = static_cast<Index>(parser.get<int>("ny", 1024));
    grid.nz = static_cast<Index>(parser.get<int>("nz", 1));
    grid.haloWidth = parser.get<int>("halo-width", 1);
    grid.dx = 1.0 / static_cast<Real>(grid.nx - 1);
    grid.dy = 1.0 / static_cast<Real>(grid.ny - 1);
    grid.dz = 1.0 / static_cast<Real>(std::max(grid.nz, static_cast<Index>(1)));

    // Solver configuration
    Index maxIter = static_cast<Index>(parser.get<int>("max-iter", 10000));
    Real tolerance = parser.get<double>("tol", 1e-6);
    std::string solverName = parser.get<std::string>("solver", "jacobi");
    std::string commMode = parser.get<std::string>("comm-mode", "p2p");
    std::string outputFormat = parser.get<std::string>("output-format", "json");
    std::string outputDir = parser.get<std::string>("output-dir", "./output");
    std::string bcType = parser.get<std::string>("bc", "dirichlet");
    int saveInterval = parser.get<int>("save-interval", 0);
    bool enableProfiling = parser.has("enable-profiling");
    bool overlapComm = parser.has("overlap-comm");

    int ompThreads = parser.get<int>("omp-threads", 0);
    if (ompThreads > 0) {
#ifdef _OPENMP
        omp_set_num_threads(ompThreads);
#endif
    }

    if (rank == 0) {
        std::string threadInfo;
#ifdef _OPENMP
        threadInfo = ", OMP threads: " + std::to_string(omp_get_max_threads());
#else
        threadInfo = ", OMP disabled";
#endif
        HYPOS_INFO("HyPoS — Hybrid Poisson Solver");
        HYPOS_INFO("Grid: " << grid.nx << "x" << grid.ny << "x" << grid.nz
                   << ", MPI procs: " << size
                   << threadInfo);
    }

    // Partition grid
    UniformPartition partition;
    SubgridInfo info = partition.partition(grid, MPI_COMM_WORLD, rank);

    // Create Cartesian communicator and get neighbors
    int dims[3], periods[3] = {0, 0, 0};
    partition.getTopologyDims(dims[0], dims[1], dims[2]);
    MPI_Comm cartComm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &cartComm);

    int left, right, down, up, back, front;
    MPI_Cart_shift(cartComm, 0, 1, &left, &right);
    MPI_Cart_shift(cartComm, 1, 1, &down, &up);
    MPI_Cart_shift(cartComm, 2, 1, &back, &front);

    // Create subdomain
    Subgrid subgrid(info.nxLocal, info.nyLocal, info.nzLocal, grid.haloWidth, cartComm);
    subgrid.setNeighbors(left, right, down, up, back, front);
    subgrid.setOffsets(info.offsetX, info.offsetY, info.offsetZ);
    subgrid.setBoundaryCondition(bcType == "neumann" ? BoundaryCondition::Neumann
                                                      : BoundaryCondition::Dirichlet);

    // Setup test problem
    setupProblem(subgrid, grid);

    // Create solver and exchanger
    std::unique_ptr<PoissonSolver> solver;
    if (solverName == "jacobi") {
        solver = std::make_unique<JacobiSolver>(overlapComm);
    } else if (solverName == "red_black_gs") {
        solver = std::make_unique<RedBlackGSSolver>();
    } else if (solverName == "cg") {
        solver = std::make_unique<CGSolver>();
    } else {
        HYPOS_ERROR("Unknown solver: " << solverName);
        return 1;
    }

    if (overlapComm && solverName != "jacobi" && rank == 0) {
        HYPOS_WARN("--overlap-comm is only supported by the jacobi solver; ignored");
    }

    if (bcType != "dirichlet" && bcType != "neumann") {
        HYPOS_ERROR("Unknown boundary condition: " << bcType);
        return 1;
    }

    std::unique_ptr<HaloExchanger> exchanger;
    if (commMode == "p2p") {
        exchanger = std::make_unique<PointToPointExchanger>();
    } else if (commMode == "collective") {
        exchanger = std::make_unique<CollectiveExchanger>();
    } else {
        HYPOS_ERROR("Unknown comm mode: " << commMode);
        return 1;
    }
    exchanger->initialize(subgrid);

    std::error_code dirError;
    std::filesystem::create_directories(outputDir, dirError);
    if (dirError) {
        HYPOS_WARN("Cannot create output directory: " << outputDir
                   << " (" << dirError.message() << ")");
    }

    // Setup I/O
    std::unique_ptr<IOBackend> io;
    if (outputFormat == "vtk") {
        io = std::make_unique<VTKIOBackend>(grid.dx, grid.dy, grid.dz);
    } else if (outputFormat == "binary") {
        io = std::make_unique<BinaryIOBackend>();
    }

    // Snapshot writer shared by the final output and --save-interval saves:
    // every rank writes its piece; rank 0 writes the parallel index.
    const auto writeSolution = [&](Index step) {
        if (!io) {
            return;
        }
        const std::string base = outputDir + "/solution_" + std::to_string(step);
        io->write(subgrid, base + "_r" + std::to_string(rank), static_cast<int>(step));

        const PieceExtent localPiece{info.offsetX, info.offsetY, info.offsetZ,
                                     subgrid.nxLocal(), subgrid.nyLocal(), subgrid.nzLocal()};
        std::vector<Index> flatPiece;
        if (rank == 0) {
            flatPiece.resize(static_cast<std::size_t>(size) * 6);
        }
        MPI_Gather(&localPiece, 6, MPI_UNSIGNED_LONG,
                   flatPiece.data(), 6, MPI_UNSIGNED_LONG, 0, cartComm);

        if (rank == 0) {
            std::vector<PieceExtent> pieces(static_cast<std::size_t>(size));
            for (int r = 0; r < size; ++r) {
                const Index* p = &flatPiece[static_cast<std::size_t>(r) * 6];
                pieces[static_cast<std::size_t>(r)] =
                    PieceExtent{p[0], p[1], p[2], p[3], p[4], p[5]};
            }
            io->writeParallelIndex(grid, base, static_cast<int>(step), pieces);
        }
    };

    if (saveInterval > 0 && io) {
        solver->setProgressCallback([&writeSolution, saveInterval](Index iteration) {
            if (iteration % static_cast<Index>(saveInterval) == 0) {
                writeSolution(iteration);
            }
        });
    }

    // Performance tracking
    Timer totalTimer;
    Timer commTimer;
    Timer computeTimer;
    totalTimer.start();

    // Main solve loop
    Index actualIter = 0;
    {
        HYPOS_PROFILE("solver_total");
        actualIter = solver->solve(subgrid, *exchanger, maxIter, tolerance);
    }
    totalTimer.stop();

    double totalTime = totalTimer.elapsedSeconds();

    // Compute performance metrics (estimate)
    double iterTime = totalTime / std::max(actualIter, Index(1));
    // Fused stencil cost per cell: 2D five-point update (5 flops) + residual
    // accumulation (3) = 8; 3D seven-point update (8) + residual (3) = 11.
    const double flopsPerCell = grid.is2D() ? 8.0 : 11.0;
    const double flopsPerIter = flopsPerCell * static_cast<double>(grid.nx * grid.ny * grid.nz);
    double flopsPerSec = flopsPerIter * actualIter / totalTime;

    // Gather timing stats from all ranks
    double localMaxTime = totalTime;
    double globalMaxTime = 0.0;
    MPI_Reduce(&localMaxTime, &globalMaxTime, 1, MPI_DOUBLE, MPI_MAX, 0, cartComm);

    // Report
    RunConfig config;
    config.nx = grid.nx;
    config.ny = grid.ny;
    config.nz = grid.nz;
    config.mpiProcs = size;
#ifdef _OPENMP
    config.ompThreads = omp_get_max_threads();
#else
    config.ompThreads = 1;
#endif
    config.solver = solverName;
    config.partition = "uniform";
    config.commMode = commMode;
    config.maxIter = maxIter;
    config.tolerance = tolerance;
    config.overlapComm = overlapComm;

    PerformanceMetrics metrics;
    metrics.totalTimeSec = totalTime;
    metrics.iterTimeMs = iterTime * 1000.0;
    metrics.iterations = actualIter;
    metrics.finalResidual = solver->lastResidual();
    metrics.flopsPerSec = flopsPerSec;

    const Profiler& profiler = Profiler::instance();
    const double haloPostSec = profiler.stats("halo_exchange").totalSeconds;
    const double haloWaitSec = profiler.stats("halo_wait").totalSeconds;
    const double commSec = haloPostSec + haloWaitSec;
    const double iterCount = static_cast<double>(std::max(actualIter, Index(1)));
    metrics.commTimeMs = commSec / iterCount * 1000.0;
    if (overlapComm && commSec > 0.0) {
        metrics.overlapRatio = std::clamp(1.0 - haloWaitSec / commSec, 0.0, 1.0);
    } else {
        metrics.overlapRatio = 0.0;
    }

    const double computeSec = profiler.stats("stencil_interior").totalSeconds +
                              profiler.stats("stencil_boundary").totalSeconds;
    metrics.computeTimeMs = computeSec / iterCount * 1000.0;
    metrics.commOverheadRatio = (metrics.iterTimeMs > 0.0)
                                    ? metrics.commTimeMs / metrics.iterTimeMs
                                    : 0.0;

    Reporter reporter(config);
    reporter.setMetrics(metrics);

    if (rank == 0) {
        HYPOS_INFO("Total time: " << totalTime << " s, "
                   << "Iterations: " << actualIter
                   << ", Est. FLOP/s: " << flopsPerSec);

        // Output performance report
        if (outputFormat == "json" || outputFormat == "csv") {
            reporter.writeToFile(outputDir + "/performance_report." + outputFormat, outputFormat);
        }

        // Output profile report
        if (enableProfiling) {
            std::cout << Profiler::instance().report() << std::endl;
        }
    }

    // Output final field (all ranks write their piece; rank 0 writes the index)
    writeSolution(actualIter);

    // Shut down resources owned by the run before MPI is finalized:
    // persistent exchange requests must be freed while MPI is still active.
    solver.reset();
    exchanger.reset();
    io.reset();
    MPI_Comm_free(&cartComm);
    return 0;
}
