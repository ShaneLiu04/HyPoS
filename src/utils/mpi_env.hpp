#pragma once

#include <mpi.h>
#include <string>
#include "core/exception.hpp"

namespace hypo {

/**
 * @brief RAII wrapper for MPI initialization and finalization.
 * Must be created before any other MPI-dependent objects.
 */
class MPIEnv {
public:
    explicit MPIEnv(int& argc, char**& argv, int requiredThreadLevel = MPI_THREAD_FUNNELED);
    ~MPIEnv();

    MPIEnv(const MPIEnv&) = delete;
    MPIEnv& operator=(const MPIEnv&) = delete;

    int rank() const noexcept { return rank_; }
    int size() const noexcept { return size_; }
    int providedThreadLevel() const noexcept { return providedThreadLevel_; }
    bool isRoot() const noexcept { return rank_ == 0; }

    MPI_Comm comm() const noexcept { return comm_; }

    /**
     * @brief Synchronize all processes. Convenience wrapper.
     */
    void barrier() const;

    /**
     * @brief Log an error string from MPI error code.
     */
    static std::string errorString(int errorCode);

    /**
     * @brief Check MPI return code and throw MPIException on error.
     */
    static void check(int mpiError, const std::string& context);

private:
    MPI_Comm comm_ = MPI_COMM_WORLD;
    int rank_ = 0;
    int size_ = 1;
    int providedThreadLevel_ = MPI_THREAD_SINGLE;
    bool finalized_ = false;
};

} // namespace hypo
