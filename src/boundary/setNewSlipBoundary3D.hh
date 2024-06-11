/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Fedor Bukreev, Adrian Kummerländer
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

#ifndef SET_NEW_SLIP_BOUNDARY_3D_HH
#define SET_NEW_SLIP_BOUNDARY_3D_HH

#include "dynamics/dynamics.h"

namespace olb {

template <typename T, typename DESCRIPTOR, int discreteNormalX, int discreteNormalY, int discreteNormalZ>
struct SlipBoundaryPostProcessor3D
{
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform{
    int mirrorDirection0;
    int mirrorDirection1;
    int mirrorDirection2;
    int mult = 2 / (discreteNormalX*discreteNormalX + discreteNormalY*discreteNormalY + discreteNormalZ*discreteNormalZ);
    int reflectionPop[DESCRIPTOR::q] = {0};
    for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
      reflectionPop[iPop] = 0;
      // iPop are the directions which pointing into the fluid, discreteNormal is pointing outwarts
      int scalarProduct = descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ;
      if (scalarProduct < 0) {
        // bounce back for the case discreteNormalX = discreteNormalY = discreteNormalZ = 1, that is mult=0
        if (mult == 0) {
          mirrorDirection0 = -descriptors::c<DESCRIPTOR>(iPop,0);
          mirrorDirection1 = -descriptors::c<DESCRIPTOR>(iPop,1);
          mirrorDirection2 = -descriptors::c<DESCRIPTOR>(iPop,2);
        }
        else {
          mirrorDirection0 = descriptors::c<DESCRIPTOR>(iPop,0) - mult*scalarProduct*discreteNormalX;
          mirrorDirection1 = descriptors::c<DESCRIPTOR>(iPop,1) - mult*scalarProduct*discreteNormalY;
          mirrorDirection2 = descriptors::c<DESCRIPTOR>(iPop,2) - mult*scalarProduct*discreteNormalZ;
        }
        // run through all lattice directions and look for match of direction
        for (int i = 1; i < DESCRIPTOR::q; i++) {
          if (descriptors::c<DESCRIPTOR>(i,0)==mirrorDirection0
              && descriptors::c<DESCRIPTOR>(i,1)==mirrorDirection1
              && descriptors::c<DESCRIPTOR>(i,2)==mirrorDirection2) {
            reflectionPop[iPop] = i;
            break;
          }
        }
      }
    }
    for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
      if (reflectionPop[iPop]!=0) {
        cell[iPop] = cell[reflectionPop[iPop]];
      }
    }
  }
};

template<typename T,typename DESCRIPTOR>
void setNewSlipBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry, int material)
{
  setNewSlipBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
void setNewSlipBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setInterpolatedVelocityBoundary");
  int _overlap = 1;
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setNewSlipBoundary<T,DESCRIPTOR>(sLattice.getBlock(iC), indicator->getBlockIndicatorF(iC),includeOuterCells);
  }
  addPoints2CommBC(sLattice,std::forward<decltype(indicator)>(indicator), _overlap);
}


template <typename T, typename DESCRIPTOR>
void setNewSlipBoundary(BlockLattice<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  using namespace boundaryhelper;
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(4,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
      Dynamics<T,DESCRIPTOR>* dynamics = nullptr;

      if(util::fabs(discreteNormal[1])+util::fabs(discreteNormal[2])+util::fabs(discreteNormal[3]) == T{1.}){
        _block.addPostProcessor(
            typeid(stage::PostCollide), {iX,iY,iZ},
            promisePostProcessorForNormal<T, DESCRIPTOR, SlipBoundaryPostProcessor3D>(
              Vector<int,3>(discreteNormal.data()+1)
            )
          );
        dynamics = _block.template getDynamics<NoCollideDynamics<T,DESCRIPTOR>>();
      } else {
        dynamics = _block.template getDynamics<BounceBack<T,DESCRIPTOR>>();
      }
      setBoundary(_block, iX, iY, iZ, dynamics);
    }
  });
}

}

#endif
