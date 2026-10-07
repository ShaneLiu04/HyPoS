#include "perf/reporter.hpp"
#include <sstream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>

namespace hypo {

Reporter::Reporter(const RunConfig& config)
    : config_(config)
    , runId_(generateRunId())
{
}

std::string Reporter::generateRunId() const {
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::tm tm;
#ifdef _WIN32
    localtime_s(&tm, &time);
#else
    localtime_r(&time, &tm);
#endif
    std::ostringstream oss;
    oss << "hypo_" << std::put_time(&tm, "%Y%m%d_%H%M%S");
    return oss.str();
}

std::string Reporter::toJson() const {
    std::ostringstream oss;
    oss << "{\n";
    oss << "  \"run_id\": \"" << runId_ << "\",\n";
    oss << "  \"config\": {\n";
    oss << "    \"nx\": " << config_.nx << ",\n";
    oss << "    \"ny\": " << config_.ny << ",\n";
    oss << "    \"nz\": " << config_.nz << ",\n";
    oss << "    \"mpi_procs\": " << config_.mpiProcs << ",\n";
    oss << "    \"omp_threads\": " << config_.ompThreads << ",\n";
    oss << "    \"solver\": \"" << config_.solver << "\",\n";
    oss << "    \"partition\": \"" << config_.partition << "\",\n";
    oss << "    \"comm_mode\": \"" << config_.commMode << "\",\n";
    oss << "    \"max_iter\": " << config_.maxIter << ",\n";
    oss << "    \"tolerance\": " << config_.tolerance << ",\n";
    oss << "    \"overlap_comm\": " << (config_.overlapComm ? "true" : "false") << "\n";
    oss << "  },\n";
    oss << "  \"performance\": {\n";
    oss << "    \"total_time_sec\": " << std::fixed << std::setprecision(4) << metrics_.totalTimeSec << ",\n";
    oss << std::defaultfloat;
    oss << "    \"iter_time_ms\": " << metrics_.iterTimeMs << ",\n";
    oss << "    \"compute_time_ms\": " << metrics_.computeTimeMs << ",\n";
    oss << "    \"comm_time_ms\": " << metrics_.commTimeMs << ",\n";
    oss << "    \"comm_overhead_ratio\": " << metrics_.commOverheadRatio << ",\n";
    oss << "    \"overlap_ratio\": " << metrics_.overlapRatio << ",\n";
    oss << "    \"flops_per_sec\": " << std::scientific << metrics_.flopsPerSec << ",\n";
    oss << "    \"iterations\": " << metrics_.iterations << ",\n";
    oss << "    \"final_residual\": " << metrics_.finalResidual << "\n";
    oss << "  }\n";
    oss << "}\n";
    return oss.str();
}

std::string Reporter::toCsv() const {
    std::ostringstream oss;
    oss << runId_ << ","
        << config_.nx << "," << config_.ny << "," << config_.nz << ","
        << config_.mpiProcs << "," << config_.ompThreads << ","
        << config_.solver << "," << config_.commMode << ","
        << config_.maxIter << "," << config_.tolerance << ","
        << (config_.overlapComm ? "1" : "0") << ","
        << metrics_.totalTimeSec << ","
        << metrics_.iterTimeMs << ","
        << metrics_.commOverheadRatio << ","
        << metrics_.iterations << ","
        << metrics_.finalResidual << "\n";
    return oss.str();
}

void Reporter::writeToFile(const std::string& path, const std::string& format) const {
    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("Cannot open file for writing: " + path);
    }
    if (format == "json") {
        ofs << toJson();
    } else if (format == "csv") {
        ofs << toCsv();
    } else {
        ofs << toJson();
    }
}

} // namespace hypo
