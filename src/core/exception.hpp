#pragma once

#include <stdexcept>
#include <string>

namespace hypo {

class HyPoSException : public std::runtime_error {
public:
    explicit HyPoSException(const std::string& msg)
        : std::runtime_error(msg) {}
};

class MPIException : public std::runtime_error {
public:
    explicit MPIException(const std::string& msg)
        : std::runtime_error(msg) {}
};

class ConfigException : public HyPoSException {
public:
    explicit ConfigException(const std::string& msg)
        : HyPoSException(msg) {}
};

} // namespace hypo

#ifndef NDEBUG
#  define HYPOS_ASSERT(cond)                                                   \
     do {                                                                      \
         if (!(cond)) {                                                        \
             throw hypo::HyPoSException("Assertion failed: " #cond           \
                                        " at " __FILE__ ":"                   \
                                        + std::to_string(__LINE__));          \
         }                                                                     \
     } while (0)
#else
#  define HYPOS_ASSERT(cond) ((void)0)
#endif
