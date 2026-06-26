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
#include <cmath>
#include <memory>
#include <string>
#include <iostream>

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
              << "  --solver <string>            Solver type: jacobi (default: jacobi)\n"
              << "  --max-iter <int>             Maximum iterations (default: 10000)\n"
              << "  --tol <double>               Convergence tolerance (default: 1e-6)\n\n"
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

    // MPI initialization with thread support for OpenMP
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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
        MPI_Finalize();
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
        HYPOS_INFO("HyPoS — Hybrid Poisson Solver");
        HYPOS_INFO("Grid: " << grid.nx << "x" << grid.ny << "x" << grid.nz
                   << ", MPI procs: " << size
#ifdef _OPENMP
                   << ", OMP threads: " << omp_get_max_threads()
#else
                   << ", OMP disabled"
#endif
        );
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

    // Setup test problem
    setupProblem(subgrid, grid);

    // Create solver and exchanger
    std::unique_ptr<PoissonSolver> solver;
    if (solverName == "jacobi") {
        solver = std::make_unique<JacobiSolver>();
    } else {
        HYPOS_ERROR("Unknown solver: " << solverName);
        MPI_Finalize();
        return 1;
    }

    std::unique_ptr<HaloExchanger> exchanger;
    if (commMode == "p2p") {
        exchanger = std::make_unique<PointToPointExchanger>();
    } else if (commMode == "collective") {
        exchanger = std::make_unique<CollectiveExchanger>();
    } else {
        HYPOS_ERROR("Unknown comm mode: " << commMode);
        MPI_Finalize();
        return 1;
    }
    exchanger->initialize(subgrid);

    // Setup I/O
    std::unique_ptr<IOBackend> io;
    if (outputFormat == "vtk") {
        io = std::make_unique<VTKIOBackend>(grid.dx, grid.dy, grid.dz);
    } else if (outputFormat == "binary") {
        io = std::make_unique<BinaryIOBackend>();
    }

    // Performance tracking
    Timer totalTimer;
    Timer commTimer;
    Timer computeTimer;
    totalTimer.start();

    // Main solve loop
    Index actualIter = 0;
    Real finalResidual = 0.0;
    {
        HYPOS_PROFILE("solver_total");
        actualIter = solver->solve(subgrid, *exchanger, maxIter, tolerance);
    }
    totalTimer.stop();

    double totalTime = totalTimer.elapsedSeconds();

    // Compute performance metrics (estimate)
    double iterTime = totalTime / std::max(actualIter, Index(1));
    double flopsPerIter = 2.0 * grid.nx * grid.ny * grid.nz; // 5-point stencil + residual
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
    metrics.finalResidual = finalResidual;
    metrics.flopsPerSec = flopsPerSec;

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

    // Output final field
    if (io && rank == 0) {
        io->write(subgrid, outputDir + "/solution", 0);
    }

    MPI_Comm_free(&cartComm);
    MPI_Finalize();
    return 0;
}
