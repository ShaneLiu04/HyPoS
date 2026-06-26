#pragma once

#include <chrono>
#include <string>

namespace hypo {

/**
 * @brief High-precision timer using std::chrono::steady_clock.
 * Also supports RDTSC-style reading when available.
 */
class Timer {
public:
    Timer() = default;

    void start() noexcept;
    void stop() noexcept;

    /**
     * @brief Elapsed time in seconds since start().
     */
    double elapsedSeconds() const noexcept;

    /**
     * @brief Elapsed time in milliseconds.
     */
    double elapsedMilliseconds() const noexcept;

    /**
     * @brief Reset and restart.
     */
    void restart() noexcept;

private:
    std::chrono::steady_clock::time_point start_;
    std::chrono::steady_clock::time_point end_;
    bool running_ = false;
};

} // namespace hypo
