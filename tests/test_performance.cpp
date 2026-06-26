#include <gtest/gtest.h>
#include "perf/timer.hpp"
#include "perf/profiler.hpp"
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
