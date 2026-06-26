#include "utils/mpi_env.hpp"
#include <cstring>

namespace hypo {

MPIEnv::MPIEnv(int& argc, char**& argv, int requiredThreadLevel) {
    int provided = 0;
    int err = MPI_Init_thread(&argc, &argv, requiredThreadLevel, &provided);
    if (err != MPI_SUCCESS) {
        throw MPIException("MPI_Init_thread failed: " + errorString(err));
    }
    providedThreadLevel_ = provided;
    comm_ = MPI_COMM_WORLD;
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
}

MPIEnv::~MPIEnv() {
    if (!finalized_) {
        MPI_Finalize();
        finalized_ = true;
    }
}

void MPIEnv::barrier() const {
    MPI_Barrier(comm_);
}

std::string MPIEnv::errorString(int errorCode) {
    char errStr[MPI_MAX_ERROR_STRING];
    int len = 0;
    MPI_Error_string(errorCode, errStr, &len);
    return std::string(errStr, len);
}

void MPIEnv::check(int mpiError, const std::string& context) {
    if (mpiError != MPI_SUCCESS) {
        throw MPIException(context + ": " + errorString(mpiError));
    }
}

} // namespace hypo
