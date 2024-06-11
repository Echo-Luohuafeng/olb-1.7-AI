/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Davide Dapelo
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

#ifndef HAIDER_LEVENSPIEL_PARTICLE_3D_HH
#define HAIDER_LEVENSPIEL_PARTICLE_3D_HH

#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <deque>

#include "HaiderLevenspielParticle3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {


template<typename T>
HaiderLevenspielParticle3D<T>::HaiderLevenspielParticle3D()
  : Particle3D<T>::Particle3D()
{ }
/*
template<typename T>
HaiderLevenspielParticle3D<T>::HaiderLevenspielParticle3D(std::vector<T> pos, T mas, T rad)
  : Particle3D<T>::Particle3D(pos, mas, rad)
{ }*/

template<typename T>
HaiderLevenspielParticle3D<T>::HaiderLevenspielParticle3D(const HaiderLevenspielParticle3D<T>& p)
  : Particle3D<T>::Particle3D(p)
{ }

template<typename T>
HaiderLevenspielParticle3D<T>::HaiderLevenspielParticle3D(std::vector<T> pos, std::vector<T> vel, T mas,
    T rad)
  : Particle3D<T>::Particle3D(pos, vel, mas, rad)
{ }

template<typename T>
inline T& HaiderLevenspielParticle3D<T>::getsphericity()
{
  return sphericity;
}

template<typename T>
inline const T& HaiderLevenspielParticle3D<T>::getsphericity() const
{
  return sphericity;
}
/*
template<typename T>
inline std::vector<T>& HaiderLevenspielParticle3D<T>::getTorque()
{
  return _torque;
}

template<typename T>
inline const std::vector<T>& HaiderLevenspielParticle3D<T>::getTorque() const
{
  return _torque;
}*/

template<typename T>
HaiderLevenspielParticle3D<T>::HaiderLevenspielParticle3D(std::vector<T> pos, T mas, T rad, T _volume , T _surface)
  : Particle3D<T>::Particle3D( pos,  mas,  std::pow(3*_volume/(4*M_PI),0.3333333333))/*{pos, mas, rad}*///, sphericity {_volume}, surface {_surface}
{
sphericity = (4*M_PI*std::pow(3*_volume/(4*M_PI),0.66666666666))/_surface; //surface of volume equivalent sphere divided by particle surface
//std::cout << "sphericity " <<sphericity << std::endl;
surface = _surface;
volume = _volume;
//std::cout <<"getmass " << this->getMass() << std::endl;
//std::cout << "surface " << surface << std::endl;
 A = util::exp(2.3288-6.4581*sphericity + 2.4486*sphericity*sphericity);
 B = 0.0964 + 0.5565*sphericity;
 C = util::exp(4.905-13.8944*sphericity + 18.4222*sphericity*sphericity-10.2599*sphericity*sphericity*sphericity);
 D = util::exp(1.4681+12.2584*sphericity-20.7322*sphericity*sphericity+15.8855*sphericity*sphericity*sphericity);
//set the radius to volume equivalent raidus
//std::cout << "radius " <<this->getRad() << std::endl;
/*std::cout << "A " <<this->getA() << std::endl;
std::cout << "B " <<this->getB() << std::endl;

std::cout << "C " <<this->getC() << std::endl;
std::cout << "D " <<this->getD() << std::endl;*/



}
/*
template<typename T>
void HaiderLevenspielParticle3D<T>::serialize(T serial[])
{
  for (int i = 0; i < 3; i++) {
    serial[i] = this->_pos[i];
    serial[i + 3] = this->_vel[i];
    serial[i + 6] = this->_force[i];
  }
  serial[9] = this->_mas;
  serial[10] = this->_rad;
  serial[11] = this->_cuboid;
  serial[12] = (double) this->_active;
  serial[13] = (double) _aVel[0];
  serial[14] = (double) _aVel[1];
  serial[15] = (double) _aVel[2];
  serial[16] = (double) _torque[0];
  serial[17] = (double) _torque[1];
  serial[18] = (double) _torque[2];
}

template<typename T>
void HaiderLevenspielParticle3D<T>::unserialize(T* data)
{
  for (int i = 0; i < 3; i++) {
    this->_pos[i] = data[i];
    this->_vel[i] = data[i + 3];
    this->_force[i] = data[i + 6];
  }
  this->_mas = data[9];
  this->_rad = data[10];
  this->_cuboid = int(data[11]);
  this->_active = (bool) data[12];
  _aVel[0] = (bool) data[13];
  _aVel[1] = (bool) data[14];
  _aVel[2] = (bool) data[15];
  _torque[0] = (bool) data[16];
  _torque[1] = (bool) data[17];
  _torque[2] = (bool) data[18];
}
*/

/////////////////////////////////////// SimulateParticles<T, HaiderLevenspielParticle3D> ///////////////////////////////////////

template<typename T>
SimulateParticles<T,HaiderLevenspielParticle3D>::SimulateParticles(ParticleSystem3D<T, HaiderLevenspielParticle3D>* ps)
  : _pSys(ps)
{ }

template<typename T>
inline void SimulateParticles<T,HaiderLevenspielParticle3D>::simulate(T dT, bool scale)
{
 //std::cout << "SimulateParticles<T,HaiderLevenspielParticle3D>::simulate" << std::endl;
  //_pSys->resetMag();
  _pSys->computeForce();
  //std::cout << "SimulateParticles<T,HaiderLevenspielParticle3D>::simulate after compute force" << std::endl;
  _pSys->explicitEuler(dT, scale);
 // _pSys->analytical(dT, scale);
  //_pSys->rungeKutta4(dT);
  //std::cout << "SimulateParticles<T,HaiderLevenspielParticle3D>::simulate after explicit euler" << std::endl;
/*
#ifdef CollisionModels
  _pSys->partialElasticImpact(0.67);
#endif*/
//std::cout << "SimulateParticles<T,HaiderLevenspielParticle3D>::simulate after compute force" << std::endl;
}

template<typename T>
inline void SimulateParticles<T,HaiderLevenspielParticle3D>::simulateWithTwoWayCoupling_Mathias ( T dT,
    ForwardCouplingModel<T,HaiderLevenspielParticle3D>& forwardCoupling,
    BackCouplingModel<T,HaiderLevenspielParticle3D>& backCoupling,
    int material, int subSteps, bool scale )
{
  for (int iSubStep=1; iSubStep<=subSteps; iSubStep++) {
    if (! _pSys->executeForwardCoupling(forwardCoupling) ) {
      std::cout << " on substep " << iSubStep << std::endl;
      singleton::exit(1);
    }
    _pSys->computeForce();
    _pSys->explicitEuler(dT/(T)(subSteps), scale);
    //_pSys->rungeKutta4(dT/(T)(subSteps));
  }
  _pSys->executeBackwardCoupling(backCoupling, material);
}

template<typename T>
inline void SimulateParticles<T,HaiderLevenspielParticle3D>::simulateWithTwoWayCoupling_Davide ( T dT,
    ForwardCouplingModel<T,HaiderLevenspielParticle3D>& forwardCoupling,
    BackCouplingModel<T,HaiderLevenspielParticle3D>& backCoupling,
    int material, int subSteps, bool scale )
{
  for (int iSubStep=1; iSubStep<=subSteps; iSubStep++) {
    if (! _pSys->executeForwardCoupling(forwardCoupling) ) {
      std::cout << " on substep " << iSubStep << std::endl;
      singleton::exit(1);
    }
    _pSys->computeForce();
    _pSys->explicitEuler(dT/(T)(subSteps), scale);
    //_pSys->rungeKutta4(dT/(T)(subSteps));
    _pSys->executeBackwardCoupling(backCoupling, material, subSteps);
  }
}

template<typename T>
inline void SimulateParticles<T,HaiderLevenspielParticle3D>::simulate(T dT, std::set<int> sActivityOfParticle, bool scale)
{
  _pSys->resetMag(sActivityOfParticle);
  _pSys->computeForce(sActivityOfParticle);
  _pSys->explicitEuler(dT, sActivityOfParticle, scale);
  _pSys->integrateTorqueMag(dT, sActivityOfParticle);

#ifdef CollisionModels
  _pSys->partialElasticImpact(0.67);
#endif

#ifdef CollisionModelsCombindedWithMechContactForce
  _pSys->partialElasticImpactForCombinationWithMechContactForce(0.67);
#endif
}

}

#endif
