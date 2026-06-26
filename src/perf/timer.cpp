#include "perf/timer.hpp"

namespace hypo {

void Timer::start() noexcept {
    start_ = std::chrono::steady_clock::now();
    running_ = true;
}

void Timer::stop() noexcept {
    end_ = std::chrono::steady_clock::now();
    running_ = false;
}

double Timer::elapsedSeconds() const noexcept {
    auto end = running_ ? std::chrono::steady_clock::now() : end_;
    return std::chrono::duration<double>(end - start_).count();
}

double Timer::elapsedMilliseconds() const noexcept {
    return elapsedSeconds() * 1000.0;
}

void Timer::restart() noexcept {
    start();
}

} // namespace hypo
