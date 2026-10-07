#include <gtest/gtest.h>
#include <mpi.h>
#include "perf/timer.hpp"
#include "perf/profiler.hpp"
#include "solver/solver.hpp"
#include "grid/subgrid.hpp"
#include "comm/halo_exchanger.hpp"
#include <thread>
#include <chrono>

using namespace hypo;

TEST(PerfTest, TimerAccuracy) {
    Timer t;
    t.start();
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    t.stop();

    double elapsed = t.elapsedMilliseconds();
    EXPECT_GT(elapsed, 90.0);  // Should be at least 90ms
    EXPECT_LT(elapsed, 200.0); // Should be less than 200ms
}

TEST(PerfTest, ProfilerRegions) {
    Profiler::instance().reset();

    {
        ProfileScope scope("test_region");
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    RegionStats stats = Profiler::instance().stats("test_region");
    EXPECT_EQ(stats.callCount, 1);
    EXPECT_GT(stats.totalSeconds, 0.005);
    EXPECT_LT(stats.totalSeconds, 0.5);
}

namespace {

void setupUniformRhs(Subgrid& sg) {
    sg.zeroInitialize();
    sg.applyDirichletBC(0.0);

    Real* rhs = sg.rhs().data();
    for (Index j = sg.jBegin(); j < sg.jEnd(); ++j) {
        for (Index i = sg.iBegin(); i < sg.iEnd(); ++i) {
            rhs[sg.index(i, j)] = -1.0;
        }
    }
}

} // namespace

TEST(PerfTest, SolveRecordsProfilerRegionsOffMode) {
    Profiler::instance().reset();

    Subgrid sg(16, 16, 1, 1, MPI_COMM_SELF);
    sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL);
    setupUniformRhs(sg);

    PointToPointExchanger ex;
    ex.initialize(sg);

    JacobiSolver solver;
    solver.solve(sg, ex, 10, 0.0);

    RegionStats halo = Profiler::instance().stats("halo_exchange");
    EXPECT_GE(halo.callCount, 10u);
    EXPECT_GT(halo.totalSeconds, 0.0);

    RegionStats interior = Profiler::instance().stats("stencil_interior");
    EXPECT_GE(interior.callCount, 10u);
    EXPECT_GT(interior.totalSeconds, 0.0);

    RegionStats boundary = Profiler::instance().stats("stencil_boundary");
    EXPECT_GE(boundary.callCount, 10u);
    EXPECT_GT(boundary.totalSeconds, 0.0);

    RegionStats residual = Profiler::instance().stats("residual_allreduce");
    EXPECT_GE(residual.callCount, 10u);
    EXPECT_GT(residual.totalSeconds, 0.0);

    RegionStats iteration = Profiler::instance().stats("jacobi_iteration");
    EXPECT_GE(iteration.callCount, 10u);
    EXPECT_GT(iteration.totalSeconds, 0.0);

    RegionStats wait = Profiler::instance().stats("halo_wait");
    EXPECT_EQ(wait.callCount, 0u);
}

TEST(PerfTest, SolveRecordsProfilerRegionsOverlapMode) {
    Profiler::instance().reset();

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Subgrid sg(16, 16, 1, 1, MPI_COMM_WORLD);
    sg.setNeighbors(rank, rank, rank, rank);
    setupUniformRhs(sg);

    PointToPointExchanger ex;
    ex.initialize(sg);

    JacobiSolver solver(true);
    solver.solve(sg, ex, 10, 0.0);

    RegionStats halo = Profiler::instance().stats("halo_exchange");
    EXPECT_GE(halo.callCount, 10u);
    EXPECT_GT(halo.totalSeconds, 0.0);

    RegionStats wait = Profiler::instance().stats("halo_wait");
    EXPECT_GE(wait.callCount, 10u);
    EXPECT_GE(wait.totalSeconds, 0.0);

    RegionStats interior = Profiler::instance().stats("stencil_interior");
    EXPECT_GE(interior.callCount, 10u);

    RegionStats boundary = Profiler::instance().stats("stencil_boundary");
    EXPECT_GE(boundary.callCount, 10u);
}
