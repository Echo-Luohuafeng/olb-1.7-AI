/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Mathias J. Krause, Davide Dapelo
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

#ifndef HAIDER_LEVENSPIEL_PARTICLE_3D_H
#define HAIDER_LEVENSPIEL_PARTICLE_3D_H

#include <set>
#include <vector>
#include <list>
#include <deque>
#include <string>
#include <iostream>
#include "particles/subgrid3DLegacyFramework/particle3D.h"




namespace olb {


/////////////////////////////////////////// RotatingParticle3D ///////////////////////////////////////////

/*
 * Rotating Particles
 */
template<typename T>
class HaiderLevenspielParticle3D: public Particle3D<T> {
public:
  HaiderLevenspielParticle3D();
 //HaiderLevenspielParticle3D(std::vector<T> pos, T mas = 1., T rad = 1.);
  HaiderLevenspielParticle3D(std::vector<T> pos, std::vector<T> vel, T mas = 1., T rad = 1.);
  HaiderLevenspielParticle3D(const HaiderLevenspielParticle3D<T>& p);
  HaiderLevenspielParticle3D(std::vector<T> pos, T mas = 1., T rad = 1., T _volume=1., T _surface = 1.);
 inline T& getsphericity();
  inline const T& getsphericity() const;
 /* inline std::vector<T>& getTorque();
  inline const std::vector<T>& getTorque() const;
  void serialize(T serial[]);
  void unserialize(T*);*/
  inline T& getA()
  {
   return A;
  };
    inline T& getB()
  {
   return B;
  };
    inline T& getC()
  {
   return C;
  };
    inline T& getD()
  {
   return D;
  };
  void set_fluidVel(T fluidVel[3])
  {
   _fluidVel[0]=fluidVel[0];
   _fluidVel[1]=fluidVel[1];
   _fluidVel[2]=fluidVel[2];
  }
  T _fluidVel[3] = {};

  static const int serialPartSize = 19;

private:
//set the radius to volume equivalent radius
T sphericity;
 T surface;
 T volume;
 //coeficients
 T A = std::exp(2.3288-6.4581*sphericity + 2.4486*sphericity*sphericity);
 T B = 0.0964 + 0.5565*sphericity;
 T C = std::exp(4.905-13.8944*sphericity + 18.4222*sphericity*sphericity-10.2599*sphericity*sphericity*sphericity);
 T D = std::exp(1.4681+12.2584*sphericity-20.7322*sphericity*sphericity+15.8855*sphericity*sphericity*sphericity);
// T _fluidVel[3] = {};
};

/////////////////////////////////////// SimulateParticles<T, HaiderLevenspielParticle3D> ///////////////////////////////////////

template<typename T>
class SimulateParticles<T, HaiderLevenspielParticle3D> {
public:
  SimulateParticles(ParticleSystem3D<T, HaiderLevenspielParticle3D>* ps);
  inline void simulate(T dT, bool scale = false);
  inline void simulateWithTwoWayCoupling_Mathias ( T dT,
      ForwardCouplingModel<T,HaiderLevenspielParticle3D>& forwardCoupling,
      BackCouplingModel<T,HaiderLevenspielParticle3D>& backCoupling,
      int material, int subSteps = 1, bool scale = false );
  inline void simulateWithTwoWayCoupling_Davide ( T dT,
      ForwardCouplingModel<T,HaiderLevenspielParticle3D>& forwardCoupling,
      BackCouplingModel<T,HaiderLevenspielParticle3D>& backCoupling,
      int material, int subSteps = 1, bool scale = false );
//   multiple collision models
  inline void simulate(T dT, std::set<int> sActivityOfParticle, bool scale = false);

private:
  ParticleSystem3D<T, HaiderLevenspielParticle3D>* _pSys;

};

}
#endif

