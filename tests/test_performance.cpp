#include <gtest/gtest.h>
#include <mpi.h>
#include "perf/timer.hpp"
#include "perf/profiler.hpp"
#include "solver/solver.hpp"
#include "grid/subgrid.hpp"
#include "comm/halo_exchanger.hpp"
#include <thread>
#include <chrono>
#include <atomic>
#include <cstdint>
#include <vector>

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

// ============================================================================
// AR004 T007: A4 profiler rework (design §4.2.2-7, D7; U8 and U11).
// These tests pin the NEW contract for the reworked Profiler:
//   - per-thread ThreadData registry: beginRegion/endRegion operate on
//     thread-local data with ZERO locking on the hot path, so interleaved
//     regions from different threads are never misattributed;
//   - stats()/report() aggregate across ALL registered threads (sum of
//     total/callCount, min/max across threads);
//   - reset() clears per-thread stats but keeps the registry;
//   - ThreadData is heap-allocated and never freed after registration, so a
//     joined thread's stats remain readable (the post-join assertions in
//     MultithreadedAggregationIsExact exercise this "ThreadData outlives its
//     thread" lifecycle contract);
//   - single-thread scenarios stay isomorphic to the old semantics;
//   - region names: RBGS solves profile "rbgs_iteration" (landed in T002)
//     and CG solves profile "cg_iteration" (landed in T004); the old
//     "jacobi_iteration" name must NOT appear for them (Jacobi keeps it).
// The current implementation keeps ONE shared active_ stack (endRegion even
// pops it without holding the mutex), so interleaved multithreaded
// begin/end misattributes regions and garbles call counts —
// MultithreadedAggregationIsExact therefore fails at runtime, the legal Red
// state for this task. IterationRegionNamesMatchSolver only captures the
// already-landed region-name contract (it passes today) and must keep
// passing after the profiler rework.
// ============================================================================

// U8: profiler multithreaded aggregation must be EXACT.
//  - Round A (distinct per-thread names): with the CURRENT shared-stack
//    implementation, interleaved endRegion calls pop whichever entry
//    another thread pushed, misattributing counts ACROSS NAMES — each
//    thread's own region ends up with a wrong count. Per-thread data
//    fixes this; each name must see exactly N calls.
//  - Round B (shared name): the sum across threads must aggregate to
//    4N (cross-thread sum semantics of the reworked stats()).
//  - Counts accumulate across rounds without reset; reset() zeroes them.
TEST(ProfilerAR004Test, MultithreadedAggregationIsExact) {
    constexpr int kThreads = 4;
    constexpr int kIterations = 200;

    const auto runRound = [](const std::string& nameForThread) {
        std::atomic<int> startGate{0};
        std::atomic<int> lockstepPhase{0};
        std::vector<std::thread> threads;
        for (int t = 0; t < kThreads; ++t) {
            threads.emplace_back([&, t]() {
                // Barrier-ish spin so all workers enter roughly together and
                // their regions interleave on the shared-stack bug.
                startGate.fetch_add(1, std::memory_order_acq_rel);
                while (startGate.load(std::memory_order_acquire) < kThreads) {
                }

                const std::string name =
                    nameForThread.empty() ? "u8_region" : nameForThread + std::to_string(t);
                volatile double sink = 0.0;
                for (int i = 0; i < kIterations; ++i) {
                    Profiler::instance().beginRegion(name);
                    sink += 1.0;  // tiny work; no sleep
                    // Hold all 4 threads inside the region EVERY iteration:
                    // unconditional lockstep makes concurrent endRegion calls
                    // certain, so the shared-stack misattribution bites
                    // deterministically (a start-only lockstep lets threads
                    // drift apart and occasionally pair up correctly).
                    lockstepPhase.fetch_add(1, std::memory_order_acq_rel);
                    while (lockstepPhase.load(std::memory_order_acquire) <
                           kThreads * (i + 1)) {
                    }
                    Profiler::instance().endRegion();
                }
                (void)sink;
            });
        }
        for (auto& th : threads) {
            th.join();
        }
    };

    Profiler::instance().reset();

    // Round A: distinct names per thread — must each see exactly N calls.
    runRound("u8_mt_");
    for (int t = 0; t < kThreads; ++t) {
        // Workers are joined before these reads: post-join readability of a
        // departed thread's stats (ThreadData outlives its thread).
        EXPECT_EQ(Profiler::instance().stats("u8_mt_" + std::to_string(t)).callCount,
                  std::uint64_t{kIterations})
            << "thread " << t << " region misattributed";
    }

    // Round B: shared name — cross-thread sum must be exact.
    runRound("");
    EXPECT_EQ(Profiler::instance().stats("u8_region").callCount,
              std::uint64_t{kThreads * kIterations});
    EXPECT_GE(Profiler::instance().stats("u8_region").totalSeconds, 0.0);

    // Single-thread isomorphism with the old semantics.
    {
        volatile double sink = 0.0;
        for (int i = 0; i < 50; ++i) {
            Profiler::instance().beginRegion("u8_single");
            sink += 1.0;
            Profiler::instance().endRegion();
        }
        (void)sink;
    }
    EXPECT_EQ(Profiler::instance().stats("u8_single").callCount,
              std::uint64_t{50});

    // Accumulation across rounds without reset(): round A again doubles the
    // per-thread-name counts; round B's shared total stays put.
    runRound("u8_mt_");
    for (int t = 0; t < kThreads; ++t) {
        EXPECT_EQ(Profiler::instance().stats("u8_mt_" + std::to_string(t)).callCount,
                  std::uint64_t{2 * kIterations})
            << "thread " << t << " region lost counts across rounds";
    }
    EXPECT_EQ(Profiler::instance().stats("u8_region").callCount,
              std::uint64_t{kThreads * kIterations});

    // reset() clears per-thread stats but keeps the registry.
    Profiler::instance().reset();
    for (int t = 0; t < kThreads; ++t) {
        EXPECT_EQ(Profiler::instance().stats("u8_mt_" + std::to_string(t)).callCount,
                  std::uint64_t{0});
    }
    EXPECT_EQ(Profiler::instance().stats("u8_region").callCount,
              std::uint64_t{0});
    EXPECT_EQ(Profiler::instance().stats("u8_single").callCount,
              std::uint64_t{0});
}

// U11: region-name contract — RBGS profiles "rbgs_iteration", CG profiles
// "cg_iteration" (renames already landed in T002/T004; capture semantics,
// passes today), and the old "jacobi_iteration" name must not appear for
// them. Jacobi keeps "jacobi_iteration". This test must keep passing after
// the profiler rework.
TEST(ProfilerAR004Test, IterationRegionNamesMatchSolver) {
    const Index kGridN = 32;
    const Index kMaxIter = 20;
    const Real kNoConvergeTol = 0.0;

    // Red-Black Gauss-Seidel: "rbgs_iteration", never "jacobi_iteration".
    {
        Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
        sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
                        MPI_PROC_NULL);
        setupUniformRhs(sg);

        PointToPointExchanger ex;
        ex.initialize(sg);

        Profiler::instance().reset();
        RedBlackGSSolver solver;
        const Index iters = solver.solve(sg, ex, kMaxIter, kNoConvergeTol);

        SCOPED_TRACE("RedBlackGSSolver region names");
        EXPECT_GT(iters, Index(0));
        EXPECT_EQ(Profiler::instance().stats("rbgs_iteration").callCount,
                  std::uint64_t(iters));
        EXPECT_EQ(Profiler::instance().stats("jacobi_iteration").callCount,
                  std::uint64_t(0));
    }

    // Conjugate Gradient: "cg_iteration", never "jacobi_iteration".
    {
        Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
        sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
                        MPI_PROC_NULL);
        setupUniformRhs(sg);

        PointToPointExchanger ex;
        ex.initialize(sg);

        Profiler::instance().reset();
        CGSolver solver;
        const Index iters = solver.solve(sg, ex, kMaxIter, kNoConvergeTol);

        SCOPED_TRACE("CGSolver region names");
        EXPECT_GT(iters, Index(0));
        EXPECT_EQ(Profiler::instance().stats("cg_iteration").callCount,
                  std::uint64_t(iters));
        EXPECT_EQ(Profiler::instance().stats("jacobi_iteration").callCount,
                  std::uint64_t(0));
    }

    // Jacobi keeps "jacobi_iteration" (preserved behavior).
    {
        Subgrid sg(kGridN, kGridN, 1, 1, MPI_COMM_SELF);
        sg.setNeighbors(MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
                        MPI_PROC_NULL);
        setupUniformRhs(sg);

        PointToPointExchanger ex;
        ex.initialize(sg);

        Profiler::instance().reset();
        JacobiSolver solver;
        const Index iters = solver.solve(sg, ex, kMaxIter, kNoConvergeTol);

        SCOPED_TRACE("JacobiSolver region names");
        EXPECT_GT(iters, Index(0));
        EXPECT_EQ(Profiler::instance().stats("jacobi_iteration").callCount,
                  std::uint64_t(iters));
    }
}
