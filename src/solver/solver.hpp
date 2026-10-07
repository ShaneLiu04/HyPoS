#pragma once

#include "core/types.hpp"
#include "grid/subgrid.hpp"
#include <mpi.h>
#include <string>

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

    /**
     * @brief Global L2 residual after the most recent iteration.
     * Returns 0 before any iteration has been performed.
     */
    virtual Real lastResidual() const noexcept { return 0.0; }
};

/**
 * @brief Classic Jacobi solver with OpenMP acceleration.
 * Supports optional communication-computation overlap: the interior box is
 * updated while halo messages are in flight, followed by the boundary band.
 */
class JacobiSolver : public PoissonSolver {
public:
    explicit JacobiSolver(bool overlapComm = false);

    Index solve(Subgrid& subgrid,
                HaloExchanger& exchanger,
                Index maxIter,
                Real tolerance) override;

    Real iterate(Subgrid& subgrid, HaloExchanger& exchanger) override;

    std::string name() const override { return "jacobi"; }

    bool overlapEnabled() const noexcept { return overlapComm_; }

    Real lastResidual() const noexcept override { return lastResidual_; }

private:
    /**
     * @brief Update the given interior region (writes uNext, reads u).
     * @param k0,k1 Active layer range: [0,1) for 2D, [kBegin,kEnd) for 3D.
     */
    void updateRegion(Subgrid& subgrid,
                      Index i0, Index i1,
                      Index j0, Index j1,
                      Index k0, Index k1) const;

    /**
     * @brief Independent full-interior residual scan (identical for both
     * overlap modes so iteration counts match bit-exactly).
     */
    Real residualPass(Subgrid& subgrid) const;

    bool overlapComm_ = false;
    Real lastResidual_ = 0.0;
};

} // namespace hypo
