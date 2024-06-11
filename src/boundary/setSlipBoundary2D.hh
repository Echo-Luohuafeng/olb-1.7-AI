/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz
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

//This file contains the Slip Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_SLIP_BOUNDARY_2D_HH
#define SET_SLIP_BOUNDARY_2D_HH

#include "setSlipBoundary2D.h"



namespace olb {
template <typename,typename,int NX, int NY>
struct FullSlipBoundaryPostProcessor2D{
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return -1;
  }


  template <typename CELL, typename V = typename CELL::value_t>
  void apply(CELL& x_b) any_platform {
    using DESCRIPTOR = typename CELL::descriptor_t;
    int reflectionPop[DESCRIPTOR::q];
    int mirrorDirection0;
    int mirrorDirection1;
    int mult = 2 / (NX*NX + NY*NY);
    reflectionPop[0] =0;
    for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
      reflectionPop[iPop] = 0;
      // iPop are the directions which pointing into the fluid, discreteNormal is pointing outwarts
      int scalarProduct = descriptors::c<DESCRIPTOR>(iPop,0)*NX + descriptors::c<DESCRIPTOR>(iPop,1)*NY;
      if ( scalarProduct < 0) {
        // bounce back for the case discreteNormalX = discreteNormalY = 1, that is mult=1
        if (mult == 1) {
          mirrorDirection0 = -descriptors::c<DESCRIPTOR>(iPop,0);
          mirrorDirection1 = -descriptors::c<DESCRIPTOR>(iPop,1);
        }
        else {
          mirrorDirection0 = descriptors::c<DESCRIPTOR>(iPop,0) - mult*scalarProduct*NX;
          mirrorDirection1 = descriptors::c<DESCRIPTOR>(iPop,1) - mult*scalarProduct*NY;
        }

        // run through all lattice directions and look for match of direction
        for (int i = 1; i < DESCRIPTOR::q; i++) {
          if (descriptors::c<DESCRIPTOR>(i,0)==mirrorDirection0
              && descriptors::c<DESCRIPTOR>(i,1)==mirrorDirection1) {
            reflectionPop[iPop] = i;
            break;
          }
        }
      }
    }
    for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
      if (reflectionPop[iPop]!=0) {
        //do reflection
        x_b[iPop] = x_b[reflectionPop[iPop]];
      }
    }
  }

};

///Initialising the SlipBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setSlipBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry, int material)
{
  setSlipBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material));
}

///Initialising the SlipBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setSlipBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setSlipBoundary");
  int _overlap = 1;
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setSlipBoundary<T, DESCRIPTOR>(sLattice.getBlock(iCloc),
                                   indicator->getBlockIndicatorF(iCloc), includeOuterCells);
  }
  //defined inside setLocalVelocityBoundary2D.h/hh
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC<T, DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}

///Set slip boundary on blockLattice domain
template<typename T, typename DESCRIPTOR>
void setSlipBoundary(BlockLattice<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& indicator, bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setSlipBoundary");
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(3, 0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY}) >= margin
        && indicator(iX, iY)) {
      discreteNormal = indicator.getBlockGeometry().getStatistics().getType(iX, iY);
      if (discreteNormal[1]!=0 || discreteNormal[2]!=0) {
        //set slip boundary on indicated cells
        bool _output = false;
        if (_output) {
          clout << "setSlipBoundary<" << discreteNormal[1] << ","<< discreteNormal[2] << ">("  << iX << ", "<< iX << ", " << iY << ", " << iY << " )" << std::endl;
        }
        block.addPostProcessor(
            typeid(stage::PostStream), {iX,iY},
            boundaryhelper::promisePostProcessorForNormal<T, DESCRIPTOR, FullSlipBoundaryPostProcessor2D>(
              Vector <int,2> (discreteNormal.data()+1)
            )
          );
      }
      else {
        clout << "Warning: Could not addSlipBoundary (" << iX << ", " << iY << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<"), set to bounceBack" << std::endl;
        block.template defineDynamics<BounceBack>({iX, iY});
      }
    }
  });
}


}//namespace olb

#endif
