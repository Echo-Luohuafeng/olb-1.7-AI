/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Thomas Henn, Mathias J. Krause
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

#ifndef MATERIALSTLBOUNDARY3D_HH
#define MATERIALSTLBOUNDARY3D_HH

#include <set>
#include "materialSTLBoundary3D.h"

namespace olb {



template<typename T, template<typename U> class PARTICLETYPE>
MaterialSTLBoundary3D<T, PARTICLETYPE>::MaterialSTLBoundary3D(
  SuperGeometry<T,3>& sg,
  std::set<int> materials,
  STLreader<T> &stlReader)
  : Boundary3D<T, PARTICLETYPE>(),
    _sg(sg),
    _materials(materials.begin(),materials.end()),
    _stlReader(stlReader)
{

}
/*
template<typename T, template<typename U> class PARTICLETYPE>
MaterialBoundary3D<T, PARTICLETYPE>::MaterialBoundary3D(
  SuperGeometry<T,3>& sg,
  std::set<int> materials)
  : Boundary3D<T, PARTICLETYPE>(),
    _sg(sg),
    _materials(materials.begin(),materials.end()
    )
{
}



template<typename T, template<typename U> class PARTICLETYPE>
void MaterialBoundary3D<T, PARTICLETYPE>::applyBoundary(
  typename std::deque<PARTICLETYPE<T> >::iterator& p,
  ParticleSystem3D<T, PARTICLETYPE>& psSys)
{
  int latticeR[3] = { 0 };
  _sg.getCuboidGeometry().get(p->getCuboid()).getFloorLatticeR(latticeR, &p->getPos()[0]);
  // Read only access to the material numbers of nodes around particle position
  const BlockGeometry<T,3>& bg = _sg.getBlockGeometry(
      _sg.getLoadBalancer().loc(p->getCuboid()));
  // + overlap is because of lower boundaries, latticeR has to be shifted up
  int iX = latticeR[0];
  int iY = latticeR[1];
  int iZ = latticeR[2];
  for (_matIter = _materials.begin(); _matIter != _materials.end(); _matIter++) {
    if (bg.get({iX, iY, iZ}) == *_matIter ||
        bg.get({iX, iY+1, iZ}) == *_matIter ||
        bg.get({iX, iY, iZ+1}) == *_matIter ||
        bg.get({iX, iY+1, iZ+1}) == *_matIter ||
        bg.get({iX+1, iY, iZ}) == *_matIter ||
        bg.get({iX+1, iY+1, iZ}) == *_matIter ||
        bg.get({iX+1, iY, iZ+1}) == *_matIter ||
        bg.get({iX+1, iY+1, iZ+1}) == *_matIter
       ) {
      p->setActive(false);
      return;
    }
  }
}

}*/

template<typename T, template<typename U> class PARTICLETYPE>
void MaterialSTLBoundary3D<T, PARTICLETYPE>::applyBoundary(
  typename std::deque<PARTICLETYPE<T> >::iterator& p,
  ParticleSystem3D<T, PARTICLETYPE>& psSys)
{

    /*Particle3D<T>& particle= &p;
    auto &position = particle->getPos();
    std::cout << position << std::endl;*/

  int latticeR[3] = { 0 };
  _sg.getCuboidGeometry().get(p->getCuboid()).getFloorLatticeR(latticeR, &p->getPos()[0]);
  //std::cout <<"particle position " <<  p->getPos()[0] << '\t' << p->getPos()[1] << '\t' << p->getPos()[2] << std::endl;
 // std::cout << "output stlReader " << _stlReader.signedDistance(Vector<T,3>(p->getPos())) << std::endl;

  // Read only access to the material numbers of nodes around particle position

  const BlockGeometry<T,3>& bg = _sg.getBlockGeometry(
      _sg.getLoadBalancer().loc(p->getCuboid()));
  // + overlap is because of lower boundaries, latticeR has to be shifted up
  int iX = latticeR[0];
  int iY = latticeR[1];
  int iZ = latticeR[2];
  for (_matIter = _materials.begin(); _matIter != _materials.end(); _matIter++) {
    if (bg.get({iX, iY, iZ}) == *_matIter ||
        bg.get({iX, iY+1, iZ}) == *_matIter ||
        bg.get({iX, iY, iZ+1}) == *_matIter ||
        bg.get({iX, iY+1, iZ+1}) == *_matIter ||
        bg.get({iX+1, iY, iZ}) == *_matIter ||
        bg.get({iX+1, iY+1, iZ}) == *_matIter ||
        bg.get({iX+1, iY, iZ+1}) == *_matIter ||
        bg.get({iX+1, iY+1, iZ+1}) == *_matIter
       ) {
    if(p->getActive()){
            //std::vector<std::vector<T>> direction//std::cout << "signed function called" << std::endl;
                        //for ( auto i=0; i <
//if ( std::abs(_stlReader.signedDistance(Vector<T,3>(p->getPos()))) < p->getRad())
     T distance = 0.;
 //std::cout <<"new distance function" << std::endl;
 _stlReader.distance(distance, Vector<T,3>(p->getPos()), p->getVel());
     if (  std::abs(distance) < p->getRad())
  {
      //std::cout << "signed function called ended" << std::endl;
      p->setActive(false);
      /*std::cout << "Particle " <<p->getID() << " deposited at material number " << *_matIter << std::endl;
      std::cout <<"Particle position " <<  p->getPos()[0] << '\t' << p->getPos()[1] << '\t' << p->getPos()[2] << std::endl;
      std::cout << "Output stlReader " << std::abs(distance) << std::endl;
      std::cout << "Particle Radius " << p->getRad() << std::endl;*/

  }

            }
      //p->setActive(false);
      return;
    }
  }
}

}

#endif /* MATERIALSTLBOUNDARY3D_HH */
