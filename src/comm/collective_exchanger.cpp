#include "comm/halo_exchanger.hpp"

namespace hypo {

// This file is intentionally left with minimal content.
// The CollectiveExchanger implementation is included in p2p_exchanger.cpp
// as a reference fallback to avoid linker issues with incomplete implementations.
// In a production build, this file would contain a proper MPI_Neighbor_allgatherv
// implementation using a dedicated cartesian communicator.

} // namespace hypo
