#include "comm/halo_exchanger.hpp"
#include "utils/logger.hpp"

namespace hypo {

void CollectiveExchanger::initialize(Subgrid& subgrid) {
    delegate_.initialize(subgrid);

    int rank = 0;
    MPI_Comm_rank(subgrid.comm(), &rank);
    if (rank == 0) {
        HYPOS_WARN("CollectiveExchanger currently falls back to PointToPointExchanger "
                   "(MPI_Neighbor_allgatherv is a roadmap item)");
    }
}

void CollectiveExchanger::exchange(Subgrid& subgrid, Real* data) {
    delegate_.exchange(subgrid, data);
}

void CollectiveExchanger::beginExchange(Subgrid& subgrid, Real* data) {
    delegate_.beginExchange(subgrid, data);
}

void CollectiveExchanger::endExchange(Subgrid& subgrid, Real* data) {
    delegate_.endExchange(subgrid, data);
}

} // namespace hypo
