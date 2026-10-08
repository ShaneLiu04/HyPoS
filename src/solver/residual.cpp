#include "solver/residual.hpp"

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

} // namespace hypo
