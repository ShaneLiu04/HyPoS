#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <chrono>

namespace hypo {

enum class LogLevel { TRACE, DEBUG, INFO, WARN, ERROR, FATAL };

/**
 * @brief Simple logger (single-threaded use; not thread-safe).
 * Root process logs to stdout; others can be silenced or redirected.
 */
class Logger {
public:
    static Logger& instance();

    void setLogLevel(LogLevel level) { minLevel_ = level; }
    void setRank(int rank) { rank_ = rank; }
    void enablePerRankOutput(bool enable) { perRankOutput_ = enable; }

    void log(LogLevel level, const std::string& message);

    bool isEnabled(LogLevel level) const { return level >= minLevel_; }

private:
    Logger() = default;

    LogLevel minLevel_ = LogLevel::INFO;
    int rank_ = 0;
    bool perRankOutput_ = false;

    std::string levelString(LogLevel level) const;
    std::string timestamp() const;
};

#define HYPOS_LOG(level, msg)                                                  \
    do {                                                                       \
        if (hypo::Logger::instance().isEnabled(level)) {                       \
            std::ostringstream _oss;                                           \
            _oss << msg;                                                       \
            hypo::Logger::instance().log(level, _oss.str());                   \
        }                                                                      \
    } while (0)

#define HYPOS_TRACE(msg) HYPOS_LOG(hypo::LogLevel::TRACE, msg)
#define HYPOS_DEBUG(msg) HYPOS_LOG(hypo::LogLevel::DEBUG, msg)
#define HYPOS_INFO(msg)  HYPOS_LOG(hypo::LogLevel::INFO, msg)
#define HYPOS_WARN(msg)  HYPOS_LOG(hypo::LogLevel::WARN, msg)
#define HYPOS_ERROR(msg) HYPOS_LOG(hypo::LogLevel::ERROR, msg)
#define HYPOS_FATAL(msg) HYPOS_LOG(hypo::LogLevel::FATAL, msg)

} // namespace hypo
