/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Louis Kronberg, Pavel Eichler, Stephan Simonis
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

#ifndef DYNAMICS_COLLISION_CUM_H
#define DYNAMICS_COLLISION_CUM_H

#include "cum.h"

namespace olb {

namespace collision {
  struct CUM {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "CUM";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {

    static_assert(DESCRIPTOR::d==3 && DESCRIPTOR::q == 27, "Cumulant Dynamics only implemented in D3Q27");
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      const V omega = parameters.template get<descriptors::OMEGA>();
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      V uSqr = cum<DESCRIPTOR>::cumCollision(cell, omega, rho, u);
      return {rho, uSqr};
    };
  };
};
}
}
#endif
