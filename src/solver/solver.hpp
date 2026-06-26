#pragma once

#include "core/types.hpp"
#include "grid/subgrid.hpp"
#include <mpi.h>

namespace hypo {

class HaloExchanger;

/**
 * @brief Abstract interface for Poisson equation solvers.
 */
class PoissonSolver {
public:
    virtual ~PoissonSolver() = default;

    /**
     * @brief Solve the Poisson equation on the given subdomain.
     * @param subgrid Local subdomain with current u, uNext, and rhs.
     * @param exchanger Halo exchange handler for boundary synchronization.
     * @param maxIter Maximum number of iterations.
     * @param tolerance Convergence criterion on L2 residual.
     * @return Number of iterations performed.
     */
    virtual Index solve(Subgrid& subgrid,
                        HaloExchanger& exchanger,
                        Index maxIter,
                        Real tolerance) = 0;

    /**
     * @brief Perform one iteration (for external convergence monitoring).
     * @return L2 residual after this iteration.
     */
    virtual Real iterate(Subgrid& subgrid, HaloExchanger& exchanger) = 0;

    virtual std::string name() const = 0;
};

/**
 * @brief Classic Jacobi solver with OpenMP acceleration.
 */
class JacobiSolver : public PoissonSolver {
public:
    Index solve(Subgrid& subgrid,
                HaloExchanger& exchanger,
                Index maxIter,
                Real tolerance) override;

    Real iterate(Subgrid& subgrid, HaloExchanger& exchanger) override;

    std::string name() const override { return "jacobi"; }

private:
    /**
     * @brief Compute one Jacobi update on the interior, returning residual.
     */
    Real jacobiUpdate(Subgrid& subgrid) const;
};

} // namespace hypo
