/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz, 2023 Shota Ito
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

//This file contains the AdvectionDiffusionTemperature Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_ADVECTION_DIFFUSION_TEMPERATURE_BOUNDARY_3D_HH
#define SET_ADVECTION_DIFFUSION_TEMPERATURE_BOUNDARY_3D_HH

#include "setAdvectionDiffusionTemperatureBoundary3D.h"

namespace olb {

///Initialising the setAdvectionDiffusionTemperatureBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setAdvectionDiffusionTemperatureBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry, int material)
{
  setAdvectionDiffusionTemperatureBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, superGeometry.getMaterialIndicator(material));
}

///Initialising the setAdvectionDiffusionTemperatureBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setAdvectionDiffusionTemperatureBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setAdvectionDiffusionTemperatureBoundary");
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); iCloc++) {
    setAdvectionDiffusionTemperatureBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getBlock(iCloc),
        indicator->getBlockIndicatorF(iCloc));
  }
}


/// Set AdvectionDiffusionTemperatureBoundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setAdvectionDiffusionTemperatureBoundary(BlockLattice<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator)
{
  using namespace boundaryhelper;
  const auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = 1;
  std::vector<int> discreteNormal(4,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);

      if(discreteNormal[0] == 2 || discreteNormal[0] == 4) { // internal corner or internal edge -> probably curved boundary
        throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
      }

      if (discreteNormal[0] == 0) { // flat //set momenta, dynamics and postProcessors for indicated cells
        _block.defineDynamics({iX, iY, iZ}, DirectionOrientationMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
            AdvectionDiffusionBoundariesDynamics,MixinDynamics,momenta::EquilibriumBoundaryTuple>::construct(Vector<int,3>(discreteNormal.data() + 1)));
      }
      else if (discreteNormal[0] == 1) {    // corner //set momenta, dynamics and postProcessors on indicated boundery corner cells
        _block.defineDynamics({iX, iY, iZ}, NormalMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
          AdvectionDiffusionCornerDynamics3D,MixinDynamics,momenta::EquilibriumBoundaryTuple
        >::construct(Vector<int,3>(discreteNormal.data() + 1)));
      }
      else if (discreteNormal[0] == 3) {    // edge //set momenta, dynamics and postProcessors on indicated boundary edge cells
        _block.defineDynamics({iX, iY, iZ}, NormalSpecialMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
          AdvectionDiffusionEdgesDynamics,MixinDynamics,momenta::EquilibriumBoundaryTuple
        >::construct(Vector<int,3>(discreteNormal.data() + 1)));
      }
    }
  });
}


}//namespace olb
#endif

