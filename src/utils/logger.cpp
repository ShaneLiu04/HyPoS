#include "utils/logger.hpp"
#include <iostream>
#include <iomanip>
#include <ctime>

namespace hypo {

Logger& Logger::instance() {
    static Logger instance;
    return instance;
}

void Logger::log(LogLevel level, const std::string& message) {
    if (level < minLevel_) return;

    std::ostringstream oss;
    oss << "[" << timestamp() << "] "
        << "[RANK " << rank_ << "] "
        << "[" << levelString(level) << "] "
        << message;

    std::cout << oss.str() << std::endl;
}

std::string Logger::levelString(LogLevel level) const {
    switch (level) {
        case LogLevel::TRACE: return "TRACE";
        case LogLevel::DEBUG: return "DEBUG";
        case LogLevel::INFO:  return "INFO";
        case LogLevel::WARN:  return "WARN";
        case LogLevel::ERROR: return "ERROR";
        case LogLevel::FATAL: return "FATAL";
        default: return "UNKNOWN";
    }
}

std::string Logger::timestamp() const {
    auto now = std::chrono::system_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;
    auto time = std::chrono::system_clock::to_time_t(now);
    std::tm tm;
#ifdef _WIN32
    localtime_s(&tm, &time);
#else
    localtime_r(&time, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S")
        << "." << std::setfill('0') << std::setw(3) << ms.count();
    return oss.str();
}

} // namespace hypo
