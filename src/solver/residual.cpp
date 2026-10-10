#include "solver/residual.hpp"
#include "comm/halo_exchanger.hpp"
#include <cmath>
#include <mpi.h>

namespace hypo {

Real trueResidualSquaredLocal(const Subgrid& subgrid) noexcept {
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    const Real* __restrict__ u = subgrid.u().data();
    const Real* __restrict__ rhs = subgrid.rhs().data();

    Real sum = 0.0;

    if (subgrid.nzLocal() == 1) {
        // 2D: field lives in the k=0 plane, five-point operator.
        const Real denom = 4.0;
        #pragma omp parallel for reduction(+:sum) schedule(static)
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            #pragma omp simd reduction(+:sum)
            for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                const Index idx = j * nxT + i;
                const Real neighborSum = u[idx - 1] + u[idx + 1] +
                                         u[idx - nxT] + u[idx + nxT];
                const Real r = denom * u[idx] - neighborSum + rhs[idx];
                sum += r * r;
            }
        }
    } else {
        // 3D: seven-point operator.
        const Real denom = 6.0;
        #pragma omp parallel for reduction(+:sum) schedule(static)
        for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
            for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                #pragma omp simd reduction(+:sum)
                for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                    const Index idx = (k * nyT + j) * nxT + i;
                    const Real neighborSum = u[idx - 1] + u[idx + 1] +
                                             u[idx - nxT] + u[idx + nxT] +
                                             u[idx - nxT * nyT] +
                                             u[idx + nxT * nyT];
                    const Real r = denom * u[idx] - neighborSum + rhs[idx];
                    sum += r * r;
                }
            }
        }
    }

    return sum;
}

Real globalTrueResidual(Subgrid& subgrid, HaloExchanger& exchanger) noexcept {
    exchanger.exchange(subgrid);
    const Real localSquared = trueResidualSquaredLocal(subgrid);
    Real globalSquared = 0.0;
    MPI_Allreduce(&localSquared, &globalSquared, 1, MPI_DOUBLE, MPI_SUM,
                  subgrid.comm());
    return std::sqrt(globalSquared);
}

void residualFieldLocal(const Subgrid& subgrid, Real* residual) noexcept {
    const Index nxT = subgrid.nxTotal();
    const Index nyT = subgrid.nyTotal();
    const Real* __restrict__ u = subgrid.u().data();
    const Real* __restrict__ rhs = subgrid.rhs().data();

    if (subgrid.nzLocal() == 1) {
        #pragma omp parallel for schedule(static)
        for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
            #pragma omp simd
            for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                const Index idx = j * nxT + i;
                const Real neighborSum = u[idx - 1] + u[idx + 1] +
                                         u[idx - nxT] + u[idx + nxT];
                residual[idx] = neighborSum - 4.0 * u[idx] - rhs[idx];
            }
        }
    } else {
        #pragma omp parallel for schedule(static)
        for (Index k = subgrid.kBegin(); k < subgrid.kEnd(); ++k) {
            for (Index j = subgrid.jBegin(); j < subgrid.jEnd(); ++j) {
                #pragma omp simd
                for (Index i = subgrid.iBegin(); i < subgrid.iEnd(); ++i) {
                    const Index idx = (k * nyT + j) * nxT + i;
                    const Real neighborSum = u[idx - 1] + u[idx + 1] +
                                             u[idx - nxT] + u[idx + nxT] +
                                             u[idx - nxT * nyT] +
                                             u[idx + nxT * nyT];
                    residual[idx] = neighborSum - 6.0 * u[idx] - rhs[idx];
                }
            }
        }
    }
}

} // namespace hypo
