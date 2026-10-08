#pragma once

#include "core/types.hpp"
#include <string>
#include <vector>
#include <unordered_map>

namespace hypo {

/**
 * @brief Performance report generator.
 * Aggregates timing, iteration counts, and scaling metrics into JSON/CSV.
 */
struct RunConfig {
    Index nx = 1024, ny = 1024, nz = 1;
    int mpiProcs = 1;
    int ompThreads = 1;
    std::string solver = "jacobi";
    std::string partition = "uniform";
    std::string commMode = "p2p";
    Index maxIter = kDefaultMaxIter;
    Real tolerance = kDefaultTolerance;
    bool overlapComm = false;
    int residualCheckInterval = 1;
};

struct PerformanceMetrics {
    double totalTimeSec = 0.0;
    double iterTimeMs = 0.0;
    double computeTimeMs = 0.0;
    double commTimeMs = 0.0;
    double commOverheadRatio = 0.0;
    double overlapRatio = 0.0;
    double flopsPerSec = 0.0;
    Index iterations = 0;
    Real finalResidual = 0.0;
};

class Reporter {
public:
    explicit Reporter(const RunConfig& config);

    void setMetrics(const PerformanceMetrics& metrics) { metrics_ = metrics; }

    /**
     * @brief Generate JSON report string.
     */
    std::string toJson() const;

    /**
     * @brief Generate CSV line.
     */
    std::string toCsv() const;

    /**
     * @brief Write report to file.
     */
    void writeToFile(const std::string& path, const std::string& format) const;

private:
    RunConfig config_;
    PerformanceMetrics metrics_;
    std::string runId_;

    std::string generateRunId() const;
};

} // namespace hypo
