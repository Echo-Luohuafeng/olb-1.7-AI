/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/


#ifndef ADJOINT_LBSOLVER_H
#define ADJOINT_LBSOLVER_H

#include "functors/lattice/latticeFpop3D.h"
#include "solver/lbSolver.h"


namespace olb {

namespace opti {

/** Tags different simulation modes: compute either reference simulation or
 * perform primal or dual (adjoint) simulation.
 */
enum class SolverMode : int {Reference, Primal, Dual};

/** Base class for solvers that solve both primal and dual problems.
 * Implementation is close to old Solver3D implementation.
 * So far, LATTICES is expected to hold exactly one lattice
 * (only 3D and "Navier-Stokes" has been tested).
 */
template<
  typename T,
  typename PARAMETERS,
  typename LATTICES,
  SolverMode MODE
>
class AdjointLbSolverBase : public LbSolver<T,PARAMETERS,LATTICES>
{
private:
  mutable OstreamManager            clout {std::cout, "AdjointLbSolverBase"};

public:
  using DESCRIPTOR = typename LATTICES::values_t::template get<0>;

  AdjointLbSolverBase(utilities::TypeIndexedSharedPtrTuple<PARAMETERS> params)
   : AdjointLbSolverBase::LbSolver(params)
  { }

protected:
  /// Helper for dual solver: init external fields from primal solution
  void loadPrimalPopulations()
  {
    const auto& params = this->parameters(names::Opti());
    auto lattice = std::get<0>(this->_sLattices);

    lattice->template defineField<descriptors::F>(
      this->geometry(), 1, *params.fpop);
    lattice->template defineField<descriptors::F>(
      this->geometry(), params.controlMaterial, *params.fpop);

    lattice->template defineField<descriptors::DJDF>(
      this->geometry(), 1, *(params.dObjectiveDf));
    lattice->template defineField<descriptors::DJDF>(
      this->geometry(), params.controlMaterial, *(params.dObjectiveDf));
    lattice->template defineField<descriptors::DJDALPHA>(
      this->geometry(), params.controlMaterial, *(params.dObjectiveDcontrol));

    // update fields if gpu used
    lattice->template setProcessingContext<Array<descriptors::F>>(ProcessingContext::Simulation);
    lattice->template setProcessingContext<Array<descriptors::DJDF>>(ProcessingContext::Simulation);
    lattice->template setProcessingContext<Array<descriptors::DJDALPHA>>(ProcessingContext::Simulation);
  }
};

}

}  // namespace olb

#endif
