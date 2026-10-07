#pragma once

#include "core/types.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <stdexcept>

namespace hypo {

/**
 * @brief Simple command-line argument parser.
 * Supports --key value, --key=value, --flag, and --help.
 */
class CommandLineParser {
public:
    void parse(int argc, char* argv[]);

    bool has(const std::string& key) const;

    template <typename T>
    T get(const std::string& key, const T& defaultValue) const;

    void setHelpText(const std::string& text) { helpText_ = text; }
    void printHelp() const;

    const std::vector<std::string>& positionalArgs() const { return positional_; }

private:
    std::unordered_map<std::string, std::string> args_;
    std::vector<std::string> positional_;
    std::string helpText_;
    std::string programName_;

    template <typename T>
    T convert(const std::string& str) const;
};

// Explicit specializations for common types (defined in cmdline_parser.cpp)
template <> int    CommandLineParser::get<int>   (const std::string& key, const int& defaultValue)    const;
template <> double CommandLineParser::get<double>(const std::string& key, const double& defaultValue) const;
template <> std::string CommandLineParser::get<std::string>(const std::string& key, const std::string& defaultValue) const;
template <> bool   CommandLineParser::get<bool>  (const std::string& key, const bool& defaultValue)   const;

} // namespace hypo
