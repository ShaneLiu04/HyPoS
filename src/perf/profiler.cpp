#include "perf/profiler.hpp"
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace hypo {

Profiler& Profiler::instance() {
    static Profiler p;
    return p;
}

void Profiler::beginRegion(const std::string& name) {
    std::lock_guard<std::mutex> lock(mutex_);
    active_.push_back({name, Timer()});
    active_.back().second.start();
}

void Profiler::endRegion() {
    if (active_.empty()) return;

    auto pair = active_.back();
    active_.pop_back();
    pair.second.stop();
    double elapsed = pair.second.elapsedSeconds();

    std::lock_guard<std::mutex> lock(mutex_);
    auto& s = stats_[pair.first];
    s.totalSeconds += elapsed;
    s.callCount += 1;
    if (s.callCount == 1) {
        s.minSeconds = s.maxSeconds = elapsed;
    } else {
        s.minSeconds = std::min(s.minSeconds, elapsed);
        s.maxSeconds = std::max(s.maxSeconds, elapsed);
    }
}

RegionStats Profiler::stats(const std::string& name) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = stats_.find(name);
    if (it == stats_.end()) return RegionStats{};
    return it->second;
}

void Profiler::reset() noexcept {
    std::lock_guard<std::mutex> lock(mutex_);
    stats_.clear();
    active_.clear();
}

std::string Profiler::report() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::ostringstream oss;
    oss << "{\n";
    bool first = true;
    for (const auto& kv : stats_) {
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
