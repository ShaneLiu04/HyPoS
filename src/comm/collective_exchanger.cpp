#include "comm/halo_exchanger.hpp"
#include "utils/logger.hpp"

namespace hypo {

void CollectiveExchanger::initialize(Subgrid& subgrid) {
    delegate_.initialize(subgrid);
    HYPOS_WARN("CollectiveExchanger currently falls back to PointToPointExchanger "
               "(MPI_Neighbor_allgatherv is a roadmap item)");
}

void CollectiveExchanger::exchange(Subgrid& subgrid) {
    delegate_.exchange(subgrid);
}

void CollectiveExchanger::beginExchange(Subgrid& subgrid) {
    delegate_.beginExchange(subgrid);
}

void CollectiveExchanger::endExchange(Subgrid& subgrid) {
    delegate_.endExchange(subgrid);
}

} // namespace hypo
