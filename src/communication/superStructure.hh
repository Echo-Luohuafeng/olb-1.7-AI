/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Peter Weisbrod, Albert Mink, Mathias J. Krause
 *                2021 Adrian Kummerlaender
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

#ifndef SUPER_STRUCTURE_HH
#define SUPER_STRUCTURE_HH

#include "communication/superStructure.h"

namespace olb {


template<typename T, unsigned D>
SuperStructure<T,D>::SuperStructure(CuboidGeometry<T,D>& cuboidGeometry,
                                    LoadBalancer<T>& loadBalancer,
                                    int overlap)
  : _cuboidGeometry(cuboidGeometry),
    _loadBalancer(loadBalancer),
    _overlap(overlap),
    clout(std::cout, "SuperGeometry" + std::to_string(D) + "D")
{
}

template<typename T, unsigned D>
SuperStructure<T,D>::SuperStructure(int overlap)
  : SuperStructure(*(new CuboidGeometry<T,D> ()),
                   *(new LoadBalancer<T> ()),
                   overlap)
{ }

template<typename T, unsigned D>
CuboidGeometry<T,D>& SuperStructure<T,D>::getCuboidGeometry()
{
  return _cuboidGeometry;
}

template<typename T, unsigned D>
CuboidGeometry<T,D> const& SuperStructure<T,D>::getCuboidGeometry() const
{
  return _cuboidGeometry;
}

template<typename T, unsigned D>
int SuperStructure<T,D>::getOverlap()
{
  return _overlap;
}

template<typename T, unsigned D>
int SuperStructure<T,D>::getOverlap() const
{
  return _overlap;
}

template<typename T, unsigned D>
LoadBalancer<T>& SuperStructure<T,D>::getLoadBalancer()
{
  return _loadBalancer;
}

template<typename T, unsigned D>
LoadBalancer<T> const& SuperStructure<T,D>::getLoadBalancer() const
{
  return _loadBalancer;
}

template<typename T, unsigned D>
template <typename F>
void SuperStructure<T,D>::forCorePhysLocations(F f) const
{
  using loc = typename PhysR<T,D>::value_t;
  Vector<loc,D> minPhysR = _cuboidGeometry.getMinPhysR();
  Vector<loc,D> maxPhysR = _cuboidGeometry.getMaxPhysR();
  const loc L = _cuboidGeometry.getMotherCuboid().getDeltaR();
  for (loc iX=minPhysR[0]; iX < maxPhysR[0]; iX+=L) {
    for (loc iY=minPhysR[1]; iY < maxPhysR[1]; iY+=L) {
      if constexpr (D == 3) {
        for (loc iZ=minPhysR[2]; iZ < maxPhysR[2]; iZ+=L) {
          if constexpr (std::is_invocable_v<F, PhysR<T,D>>) {
            f({iX,iY,iZ});
          } else {
            f(iX,iY,iZ);
          }
        }
      } else {
        if constexpr (std::is_invocable_v<F, PhysR<T,D>>) {
          f({iX,iY});
        } else {
          f(iX,iY);
        }
      }
    }
  }
};

template<typename T, unsigned D>
template <typename F>
void SuperStructure<T,D>::forCorePhysLocations(PhysR<T,D> min, PhysR<T,D> max, F f) const
{
  using loc = typename PhysR<T,D>::value_t;
  Vector<loc,D> minPhysR = _cuboidGeometry.getMinPhysR();
  Vector<loc,D> maxPhysR = _cuboidGeometry.getMaxPhysR();
  const loc L = _cuboidGeometry.getMotherCuboid().getDeltaR();
  for (loc iX=std::max(minPhysR[0],min[0]); iX < std::min(maxPhysR[0],max[0]+L); iX+=L) {
    for (loc iY=std::max(minPhysR[1],min[1]); iY < std::min(maxPhysR[1],max[1]+L); iY+=L) {
      if constexpr (D == 3) {
        for (loc iZ=std::max(minPhysR[2],min[2]); iZ < std::min(maxPhysR[2],max[2]+L); iZ+=L) {
          if constexpr (std::is_invocable_v<F, PhysR<T,D>>) {
            f({iX,iY,iZ});
          } else {
            f(iX,iY,iZ);
          }
        }
      } else {
        if constexpr (std::is_invocable_v<F, PhysR<T,D>>) {
          f({iX,iY});
        } else {
          f(iX,iY);
        }
      }
    }
  }
};

template<typename T, unsigned D>
template <typename F>
void SuperStructure<T,D>::forCoreSpatialLocations(F f) const
{
  forCorePhysLocations([&](PhysR<T,D> physLoc){
    LatticeR<D+1> latticeR;
    _cuboidGeometry.getLatticeR( latticeR.data(), physLoc.data() );
    if constexpr (std::is_invocable_v<F, LatticeR<D+1>>) {
      f(latticeR);
    } else {
      if constexpr (D == 3) {
        f(latticeR[0],latticeR[1],latticeR[2],latticeR[3]);
      } else {
        f(latticeR[0],latticeR[1],latticeR[2]);
      }
    }
  });
};

template<typename T, unsigned D>
template <typename F>
void SuperStructure<T,D>::forCoreSpatialLocations(PhysR<T,D> min, PhysR<T,D> max, F f) const
{
  forCorePhysLocations(min, max, [&](PhysR<T,D> physLoc){
    LatticeR<D+1> latticeR;
    _cuboidGeometry.getLatticeR( latticeR.data(), physLoc.data() );
    if constexpr (std::is_invocable_v<F, LatticeR<D+1>>) {
      f(latticeR);
    } else {
      if constexpr (D == 3) {
        f(latticeR[0],latticeR[1],latticeR[2],latticeR[3]);
      } else {
        f(latticeR[0],latticeR[1],latticeR[2]);
      }
    }
  });
};

}

#endif
