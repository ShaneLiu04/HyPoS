#pragma once

// FP4 (AR008 design D5): the four CG kernels extracted as free functions
// so MGPreconditionedCGSolver can reuse them without copying ~120 lines.
// The bodies are moved verbatim out of CGSolver's member functions —
// CGSolver's call sites stay bit-identical (U6 golden guard), and the
// preconditioned consumer passes its own buffers.
//
// FP2/FP3 (AR009 design D2/D3): cgDotGlobal is split into cgDotLocal
// (pure kernel) + cgBlockingAllreduce (single Allreduce behind the
// "cg_blocking_allreduce" profiler region) so the P3 probe can count
// blocking reductions; the composition stays bit-identical for cg. The
// pipelined solver instead packs [rho, m] into ONE MPI_Iallreduce per
// iteration (cgIallreduce2Start/Wait, "pcg_iallreduce" region).

#include "core/types.hpp"

#include <mpi.h>

namespace hypo {

class HaloExchanger;
class Subgrid;

/** @brief Global dot product over interior points (MPI_Allreduce SUM).
 *  FP2 (AR009): composition of cgDotLocal + cgBlockingAllreduce —
 *  bit-identical to the pre-split body (U6 golden guard). */
Real cgDotGlobal(const Subgrid& subgrid, const Real* a, const Real* b);

/** @brief Local sum of products over interior points (no reduction); the
 *  loop body is moved verbatim out of the pre-split cgDotGlobal. */
Real cgDotLocal(const Subgrid& subgrid, const Real* a, const Real* b);

/** @brief Single MPI_Allreduce(SUM) of one scalar behind the
 *  "cg_blocking_allreduce" profiler region. */
Real cgBlockingAllreduce(const Subgrid& subgrid, Real local);

/** @brief Issue ONE packed MPI_Iallreduce of local[2] = [rho, m] behind
 *  the "pcg_iallreduce" profiler region (count = 2 doubles). `global` is
 *  the caller-provided recv buffer (MPI forbids send/recv aliasing and a
 *  request cannot carry the pointer, so the kernel takes it explicitly —
 *  implementation-level refinement of design §4.3, logged in tasks.md). */
void cgIallreduce2Start(const Subgrid& subgrid,
                        const Real local[2],
                        Real global[2],
                        MPI_Request* req);

/** @brief Complete a packed Iallreduce started by cgIallreduce2Start; on
 *  return the caller's global[2] buffer holds the two sums. */
void cgIallreduce2Wait(MPI_Request* req);

/** @brief y = M x interior segment; halo exchange + physical boundary of x. */
void cgMatvec(Subgrid& subgrid,
              HaloExchanger& exchanger,
              Real* x,
              Real* y);

/** @brief y += alpha * x on interior points. */
void cgAxpyInterior(const Subgrid& subgrid, Real alpha, const Real* x, Real* y);

/** @brief p = v + beta * p on interior points (v/p parameterized: CGSolver
 *  passes (r_, p_), the preconditioned CG passes (z, p)). */
void cgUpdatePInterior(const Subgrid& subgrid, const Real* v, Real* p, Real beta);

} // namespace hypo
