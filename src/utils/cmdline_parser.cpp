#include "utils/cmdline_parser.hpp"
#include <algorithm>
#include <cstdlib>

namespace hypo {

void CommandLineParser::parse(int argc, char* argv[]) {
    if (argc > 0) {
        programName_ = argv[0];
    }
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.empty() || arg[0] != '-') {
            positional_.push_back(arg);
            continue;
        }
        // Remove leading dashes
        std::size_t dashCount = 0;
        while (dashCount < arg.size() && arg[dashCount] == '-') ++dashCount;
        std::string key = arg.substr(dashCount);

        // Check for --help
        if (key == "help" || key == "h") {
            printHelp();
            std::exit(0);
        }

        // Check for =value
        std::size_t eqPos = key.find('=');
        if (eqPos != std::string::npos) {
            std::string val = key.substr(eqPos + 1);
            key = key.substr(0, eqPos);
            args_[key] = val;
        } else if (i + 1 < argc && argv[i + 1][0] != '-') {
            // Next arg is value
            args_[key] = argv[i + 1];
            ++i;
        } else {
            // Flag (no value)
            args_[key] = "true";
        }
    }
}

bool CommandLineParser::has(const std::string& key) const {
    return args_.find(key) != args_.end();
}

void CommandLineParser::printHelp() const {
    std::cout << helpText_ << std::endl;
}

// Explicit instantiations
// int
template <>
int CommandLineParser::convert<int>(const std::string& str) const {
    return std::stoi(str);
}

// double
template <>
double CommandLineParser::convert<double>(const std::string& str) const {
    return std::stod(str);
}

// string
template <>
std::string CommandLineParser::convert<std::string>(const std::string& str) const {
    return str;
}

// bool
template <>
bool CommandLineParser::convert<bool>(const std::string& str) const {
    std::string lower;
    lower.reserve(str.size());
    for (char c : str) {
        lower += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    if (lower == "true" || lower == "1" || lower == "yes" || lower == "on") return true;
    if (lower == "false" || lower == "0" || lower == "no" || lower == "off") return false;
    throw std::invalid_argument("Cannot convert to bool: " + str);
}

template <>
int CommandLineParser::get<int>(const std::string& key, const int& defaultValue) const {
    auto it = args_.find(key);
    if (it == args_.end()) return defaultValue;
    return convert<int>(it->second);
}

template <>
double CommandLineParser::get<double>(const std::string& key, const double& defaultValue) const {
    auto it = args_.find(key);
    if (it == args_.end()) return defaultValue;
    return convert<double>(it->second);
}

template <>
std::string CommandLineParser::get<std::string>(const std::string& key, const std::string& defaultValue) const {
    auto it = args_.find(key);
    if (it == args_.end()) return defaultValue;
    return convert<std::string>(it->second);
}

template <>
bool CommandLineParser::get<bool>(const std::string& key, const bool& defaultValue) const {
    auto it = args_.find(key);
    if (it == args_.end()) return defaultValue;
    return convert<bool>(it->second);
}

} // namespace hypo
