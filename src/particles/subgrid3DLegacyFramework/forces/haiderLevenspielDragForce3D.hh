/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Marie-Luise Maier, Mathias J. Krause, Sascha Janz
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

/// Schiller and Naumann drag force
#ifndef HAIDERLEVENSPIELDRAGFORCE_3D_HH
#define HAIDERLEVENSPIELDRAGFORCE_3D_HH

#include <cmath>
#include "haiderLevenspielDragForce3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
HaiderLevenspielDragForce3D<T, PARTICLETYPE, DESCRIPTOR>::HaiderLevenspielDragForce3D(SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& getVel, T dynVisc, T physDensity, UnitConverter<T, DESCRIPTOR> converter):
  Force3D<T, PARTICLETYPE>(), _getVel(getVel), _dynVisc(dynVisc),  _physDensity(physDensity), _converter(converter)
{
//  this->_name = "HaiderLevenspielDragForce3D";
}

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void HaiderLevenspielDragForce3D<T, PARTICLETYPE, DESCRIPTOR>::applyForce(
  typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
  ParticleSystem3D<T, PARTICLETYPE>& pSys)
{

  Vector<T, 3> force = {T(0), T(0), T(0)} ;
  std::cout << _converter.getConversionFactorForce() << std::endl;
  T fluidVel[3] = {0., 0., 0.};
  std::cout << _physDensity << "phys density" <<std::endl;//je fyzikalni
  std::cout << _dynVisc << "phys density" <<std::endl;//je fyzikalni
  _getVel(fluidVel, &p->getPos()[0], p->getCuboid());
  std::cout << "Fluidvel: " << fluidVel[0] << " " << fluidVel[1] << " " << fluidVel[2] << std::endl;
  T rel_velo = util::sqrt(util::pow(p->getVel()[0]-fluidVel[0], 2.) + util::pow(p->getVel()[1]-fluidVel[1], 2.) + util::pow(p->getVel()[2]-fluidVel[2], 2.));
  std::cout << "rel_velo " << rel_velo << std::endl;
  //T fluidVelAbs = util::sqrt(std::pow(fluidVel[0], 2.) + util::pow(fluidVel[1], 2.) + util::pow(fluidVel[2], 2.));
  //T partVelAbs = util::sqrt(std::pow(p->getVel()[0], 2.) + util::pow(p->getVel()[1], 2.) + util::pow(p->getVel()[2], 2.));
  T particleRe = 2. * p->getRad() * rel_velo * _physDensity / _dynVisc;
std::cout << "radius " << p->getRad() << std::endl;
  std::cout << "particleRe" << particleRe << std::endl;
  T coeffCD;
   if(particleRe < 0.000001){
    coeffCD= 0;
  } else {
    coeffCD = 24./particleRe*(1. + p->getA()*util::pow(particleRe, p->getB()))+ p->getC()/(1.+p->getD()/particleRe);
  }
  std::cout << "coeffCD " << coeffCD << std::endl;


  for (int i = 0; i < 3; i++) {

    force[i] = -1. * 0.5* coeffCD * M_PI  * p->getRad()*p->getRad()*_physDensity*rel_velo*(p->getVel()[i]-fluidVel[i]);
 std::cout <<"force " << i << " " << force[i] << std::endl;
    p->getForce()[i] += force[i] ;

  }
std::cout << "apply force finished " << std::endl;
}
}

#endif // HaiderLevenspielDragForce3D
