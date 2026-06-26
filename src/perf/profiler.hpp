#pragma once

#include "perf/timer.hpp"
#include <string>
#include <unordered_map>
#include <vector>
#include <mutex>

namespace hypo {

/**
 * @brief Lightweight hierarchical profiler.
 * Records elapsed time per named region. Thread-safe via atomic counters.
 */
struct RegionStats {
    double totalSeconds = 0.0;
    std::uint64_t callCount = 0;
    double minSeconds = 0.0;
    double maxSeconds = 0.0;
};

class Profiler {
public:
    static Profiler& instance();

    /**
     * @brief Start a named region.
     */
    void beginRegion(const std::string& name);

    /**
     * @brief End the most recently started region.
     */
    void endRegion();

    /**
     * @brief Get stats for a region.
     */
    RegionStats stats(const std::string& name) const;

    /**
     * @brief Reset all stats.
     */
    void reset() noexcept;

    /**
     * @brief Dump all stats to string (JSON-like).
     */
    std::string report() const;

private:
    Profiler() = default;

    mutable std::mutex mutex_;
    std::unordered_map<std::string, RegionStats> stats_;
    std::vector<std::pair<std::string, Timer>> active_;
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
