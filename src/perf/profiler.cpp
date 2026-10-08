#include "perf/profiler.hpp"
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace hypo {

Profiler& Profiler::instance() {
    static Profiler p;
    return p;
}

Profiler::ThreadData& Profiler::threadData() const {
    // One ThreadData per thread, heap-allocated and registered exactly once.
    // Never freed: stats of a joined thread must stay readable, and the
    // registry would dangle otherwise. Reclaimed at process teardown.
    thread_local ThreadData& td = [this]() -> ThreadData& {
        ThreadData* p = new ThreadData();
        std::lock_guard<std::mutex> lock(registryMutex_);
        registry_.push_back(p);
        return *p;
    }();
    return td;
}

void Profiler::beginRegion(const std::string& name) {
    // Hot path: thread-local only, no lock.
    ThreadData& td = threadData();
    td.active.push_back({name, Timer()});
    td.active.back().second.start();
}

void Profiler::endRegion() {
    // Hot path: thread-local only, no lock.
    ThreadData& td = threadData();
    if (td.active.empty()) return;

    auto pair = td.active.back();
    td.active.pop_back();
    pair.second.stop();
    double elapsed = pair.second.elapsedSeconds();

    auto& s = td.stats[pair.first];
    s.totalSeconds += elapsed;
    s.callCount += 1;
    if (s.callCount == 1) {
        s.minSeconds = s.maxSeconds = elapsed;
    } else {
        s.minSeconds = std::min(s.minSeconds, elapsed);
        s.maxSeconds = std::max(s.maxSeconds, elapsed);
    }
}

void Profiler::foldInto(RegionStats& agg, const RegionStats& s) {
    agg.totalSeconds += s.totalSeconds;
    const bool firstContributor = (agg.callCount == 0);
    agg.callCount += s.callCount;
    if (s.callCount > 0) {
        if (firstContributor) {
            agg.minSeconds = s.minSeconds;
            agg.maxSeconds = s.maxSeconds;
        } else {
            agg.minSeconds = std::min(agg.minSeconds, s.minSeconds);
            agg.maxSeconds = std::max(agg.maxSeconds, s.maxSeconds);
        }
    }
}

RegionStats Profiler::stats(const std::string& name) const {
    std::lock_guard<std::mutex> lock(registryMutex_);
    RegionStats agg;
    for (const ThreadData* td : registry_) {
        auto it = td->stats.find(name);
        if (it == td->stats.end()) continue;
        foldInto(agg, it->second);
    }
    return agg;
}

void Profiler::reset() noexcept {
    // Caller contract: no thread is inside a region while reset() runs.
    std::lock_guard<std::mutex> lock(registryMutex_);
    for (ThreadData* td : registry_) {
        td->stats.clear();
        td->active.clear();
    }
}

std::string Profiler::report() const {
    std::lock_guard<std::mutex> lock(registryMutex_);
    std::unordered_map<std::string, RegionStats> combined;
    for (const ThreadData* td : registry_) {
        for (const auto& kv : td->stats) {
            foldInto(combined[kv.first], kv.second);
        }
    }

    std::ostringstream oss;
    oss << "{\n";
    bool first = true;
    for (const auto& kv : combined) {
        if (!first) oss << ",\n";
        first = false;
        const auto& s = kv.second;
        oss << "  \"" << kv.first << "\": {\n"
            << "    \"total_sec\": " << std::fixed << std::setprecision(6) << s.totalSeconds << ",\n"
            << "    \"calls\": " << s.callCount << ",\n"
            << "    \"avg_ms\": " << (s.callCount > 0 ? (s.totalSeconds / s.callCount * 1000.0) : 0.0) << ",\n"
            << "    \"min_ms\": " << s.minSeconds * 1000.0 << ",\n"
            << "    \"max_ms\": " << s.maxSeconds * 1000.0 << "\n"
            << "  }";
    }
    oss << "\n}\n";
    return oss.str();
}

ProfileScope::ProfileScope(const std::string& name) : name_(name) {
    Profiler::instance().beginRegion(name_);
}

ProfileScope::~ProfileScope() {
    Profiler::instance().endRegion();
}

} // namespace hypo
