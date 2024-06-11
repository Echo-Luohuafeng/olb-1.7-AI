/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Fedor Bukreev, Adrian Kummerländer,
 *                2024 Marc Heinzelmann
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

//This file contains the RobinBoundary

#ifndef SET_ROBIN_BOUNDARY_3D_HH
#define SET_ROBIN_BOUNDARY_3D_HH

#include "setRobinBoundary3D.h"

#include "boundary/postprocessor/robinBoundaryLatticePostProcessor3D.h"


namespace olb {

///Initialising the setRobinBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setRobinBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega, SuperGeometry<T,3>& superGeometry, int material)
{
  setRobinBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material));
}

///Initialising the setRobinBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setRobinBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setRobinBoundary");
  bool includeOuterCells = false;
  bool useOtherStrategy = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); iCloc++) {
    setRobinBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getBlock(iCloc), indicator->getBlockIndicatorF(iCloc), omega, includeOuterCells, useOtherStrategy);
  }
}


/// Set RobinBoundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setRobinBoundary(BlockLattice<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator, T omega, bool includeOuterCells, bool useOtherStrategy)
{
  using namespace boundaryhelper;
  const auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(4,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
      Dynamics<T,DESCRIPTOR>* dynamics = nullptr;

      if(discreteNormal[0] == 0){ //flat
        if(!useOtherStrategy){
          _block.addPostProcessor(
            typeid(stage::PostStream), {iX,iY,iZ},
            promisePostProcessorForNormal<T, DESCRIPTOR, robinBoundaryLatticePostProcessor3D>(
              Vector <int,3> (discreteNormal.data()+1)
            )
          );
        } else{
          _block.addPostProcessor(
            typeid(stage::PostStream), {iX,iY,iZ},
            promisePostProcessorForNormal<T, DESCRIPTOR, robinBoundaryLatticePostProcessor3Dother>(
              Vector <int,3> (discreteNormal.data()+1)
            )
          );
        }
        dynamics = _block.template getDynamics<AdvectionDiffusionBGKdynamics<T,DESCRIPTOR>>();
      }
      setBoundary(_block, iX, iY, iZ, dynamics);
    }
  });
}


}//namespace olb
#endif
