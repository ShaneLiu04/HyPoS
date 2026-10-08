#pragma once

#include "perf/timer.hpp"
#include <string>
#include <unordered_map>
#include <vector>
#include <mutex>

namespace hypo {

/**
 * @brief Accumulated statistics for one named region.
 */
struct RegionStats {
    double totalSeconds = 0.0;
    std::uint64_t callCount = 0;
    double minSeconds = 0.0;
    double maxSeconds = 0.0;
};

/**
 * @brief Lightweight hierarchical profiler (AR004 A4 rework).
 *
 * Concurrency model:
 *   - beginRegion()/endRegion() operate on THREAD-LOCAL data only: the hot
 *     path takes NO lock. Each thread's ThreadData is heap-allocated on
 *     first use and registered in the singleton's registry; registration
 *     happens exactly once per thread under a mutex.
 *   - ThreadData is intentionally NEVER freed: a departed thread's stats
 *     remain readable after join (the registry would otherwise dangle).
 *     The singleton reclaims everything at process teardown — a documented
 *     trade-off, not a leak.
 *   - stats()/report() take the registry mutex and AGGREGATE across all
 *     registered threads (total/callCount summed, min/max folded).
 *   - reset() clears every thread's stats and active stack but keeps the
 *     registry. The caller must ensure no thread is inside a region while
 *     reset() runs (single-threaded phase boundary is the intended use).
 */
class Profiler {
public:
    static Profiler& instance();

    /**
     * @brief Start a named region on the calling thread.
     */
    void beginRegion(const std::string& name);

    /**
     * @brief End the most recently started region on the calling thread.
     */
    void endRegion();

    /**
     * @brief Get stats for a region, aggregated over all threads.
     */
    RegionStats stats(const std::string& name) const;

    /**
     * @brief Reset all per-thread stats (registry is kept).
     */
    void reset() noexcept;

    /**
     * @brief Dump all stats to string (JSON-like), aggregated over threads.
     */
    std::string report() const;

private:
    Profiler() = default;

    struct ThreadData {
        std::vector<std::pair<std::string, Timer>> active;
        std::unordered_map<std::string, RegionStats> stats;
    };

    ThreadData& threadData() const;

    static void foldInto(RegionStats& agg, const RegionStats& s);

    mutable std::mutex registryMutex_;
    mutable std::vector<ThreadData*> registry_;
};

/**
 * @brief RAII scope guard for profiling.
 */
class ProfileScope {
public:
    explicit ProfileScope(const std::string& name);
    ~ProfileScope();

    ProfileScope(const ProfileScope&) = delete;
    ProfileScope& operator=(const ProfileScope&) = delete;

private:
    std::string name_;
};

#define HYPOS_PROFILE(name) hypo::ProfileScope _hypo_profile_scope(name)

} // namespace hypo
