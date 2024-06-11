/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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


/* This file contains functions used for the calculation of particle related dynamics.
 *
*/

#ifndef PARTICLE_MOTION_FUNCTIONS_H
#define PARTICLE_MOTION_FUNCTIONS_H


namespace olb {

namespace particles {

namespace dynamics {

/////////// Basic Newtons Movement Functions /////////

// Velocity verlet integration according to
// Verlet1967 (https://doi.org/10.1103/PhysRev.159.98)
// Swope1982 (https://doi.org/10.1063/1.442716)
template<typename T, typename PARTICLETYPE>
void velocityVerletTranslation( Particle<T,PARTICLETYPE>& particle, T delTime, T delTime2,
                                Vector<T,PARTICLETYPE::d> acceleration )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  //Check field existence
  static_assert(PARTICLETYPE::template providesNested<GENERAL,POSITION>(), "Field POSITION has to be provided");
  static_assert(PARTICLETYPE::template providesNested<MOBILITY,VELOCITY>(), "Field VELOCITY has to be provided");
  static_assert(PARTICLETYPE::template providesNested<MOBILITY,ACCELERATION_STRD>(), "Field ACCELERATION_STRD has to be provided");
  //Calculate quantities
  Vector<T,D> position( particle.template getField<GENERAL,POSITION>()
                        + particle.template getField<MOBILITY,VELOCITY>() * delTime
                        + (0.5 * particle.template getField<MOBILITY,ACCELERATION_STRD>() * delTime2) );
  Vector<T,D> avgAcceleration( .5*(particle.template getField<MOBILITY,ACCELERATION_STRD>() + acceleration) );
  Vector<T,D> velocity( particle.template getField<MOBILITY,VELOCITY>() + avgAcceleration * delTime );
  //Update values
  particle.template setField<GENERAL,POSITION>( position );
  particle.template setField<MOBILITY,VELOCITY>( velocity );
  particle.template setField<MOBILITY,ACCELERATION_STRD>( acceleration );
}

template<typename T, typename PARTICLETYPE>
void velocityVerletRotation( Particle<T,PARTICLETYPE>& particle, T delTime, T delTime2,
                                 Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angularAcceleration )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  constexpr unsigned Drot = utilities::dimensions::convert<D>::rotation;
  //Check field existence
  static_assert(PARTICLETYPE::template providesNested<SURFACE,ANGLE>(), "Field SURFACE:ANGLE has to be provided");
  static_assert(PARTICLETYPE::template providesNested<MOBILITY,ANG_VELOCITY>(), "Field MOBILITY:ANG_VELOCITY has to be provided");
  static_assert(PARTICLETYPE::template providesNested<MOBILITY,ANG_ACC_STRD>(), "Field MOBILITY:ANG_ACC_STRD has to be provided");
  //Calculate quantities
  Vector<T,Drot> angle( particle.template getField<SURFACE,ANGLE>()
                        + particle.template getField<MOBILITY,ANG_VELOCITY>() * delTime
                        + (0.5 * particle.template getField<MOBILITY,ANG_ACC_STRD>() * delTime2) );
  Vector<T,Drot> avgAngularAcceleration( .5*(particle.template getField<MOBILITY,ANG_ACC_STRD>() + angularAcceleration) );
  Vector<T,Drot> angularVelocity( particle.template getField<MOBILITY,ANG_VELOCITY>() + avgAngularAcceleration * delTime );
  //Treat full rotation
  angle = util::fmod( angle, 2.*M_PI );
  //Update values
  particle.template setField<SURFACE,ANGLE>(
    utilities::dimensions::convert<D>::serialize_rotation(angle) );
  particle.template setField<MOBILITY,ANG_VELOCITY>(
    utilities::dimensions::convert<D>::serialize_rotation(angularVelocity) );
  particle.template setField<MOBILITY,ANG_ACC_STRD>(
    utilities::dimensions::convert<D>::serialize_rotation(angularAcceleration) );

}

template<typename T, typename PARTICLETYPE>
void velocityVerletRotor( Particle<T,PARTICLETYPE>& particle, T delTime, Vector<T,PARTICLETYPE::d> angVel )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  constexpr unsigned Drot = utilities::dimensions::convert<D>::rotation;
  //Check field existence
  static_assert(PARTICLETYPE::template providesNested<SURFACE,ANGLE>(), "Field SURFACE:ANGLE has to be provided");
  static_assert(PARTICLETYPE::template providesNested<MOBILITY,ANG_VELOCITY>(), "Field MOBILITY:ANG_VELOCITY has to be provided");
  static_assert(PARTICLETYPE::template providesNested<MOBILITY,ANG_ACC_STRD>(), "Field MOBILITY:ANG_ACC_STRD has to be provided");
  //Calculate quantities
  Vector<T,Drot> angle( particle.template getField<SURFACE,ANGLE>() + angVel * delTime);
  Vector<T,Drot> avgAngularAcceleration({0.});
  Vector<T,Drot> angularVelocity( angVel );
  //Treat full rotation
  angle = util::fmod( angle, 2.*M_PI );
  //Update values
  particle.template setField<SURFACE,ANGLE>(
    utilities::dimensions::convert<D>::serialize_rotation(angle) );
  particle.template setField<MOBILITY,ANG_VELOCITY>(
    utilities::dimensions::convert<D>::serialize_rotation(angularVelocity) );
  particle.template setField<MOBILITY,ANG_ACC_STRD>(
    utilities::dimensions::convert<D>::serialize_rotation(avgAngularAcceleration) );

}

template<typename T, typename PARTICLETYPE>
void velocityVerletIntegration( Particle<T,PARTICLETYPE>& particle, T delTime,
                                Vector<T,PARTICLETYPE::d> acceleration,
                                Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angularAcceleration )
{
  //Calculate squared time step
  T delTime2 = delTime*delTime;
  //Calculate cartesian directions
  velocityVerletTranslation( particle, delTime, delTime2, acceleration );
  //Calculate rotational directions
  if constexpr ( access::providesSurface<PARTICLETYPE>() ) {
    velocityVerletRotation( particle, delTime, delTime2, angularAcceleration );
  }
}





/// Euler integration
template<typename T, typename PARTICLETYPE>
void eulerIntegrationTranslation(
  Particle<T,PARTICLETYPE>& particle, T delTime,
  Vector<T,PARTICLETYPE::d> acceleration )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  //Check field existence
  static_assert(PARTICLETYPE::template providesNested<GENERAL,POSITION>(), "Field POSITION has to be provided");
  static_assert(PARTICLETYPE::template providesNested<MOBILITY,VELOCITY>(), "Field VELOCITY has to be provided");
  //Calculate quantities
  Vector<T,D> velocity( particle.template getField<MOBILITY,VELOCITY>() + acceleration * delTime );
  Vector<T,D> position( particle.template getField<GENERAL,POSITION>() + velocity * delTime );
  //Update values
  particle.template setField<GENERAL,POSITION>( position );
  particle.template setField<MOBILITY,VELOCITY>( velocity );
}

template<typename T, typename PARTICLETYPE>
void eulerIntegrationRotation(
  Particle<T,PARTICLETYPE>& particle, T delTime,
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angularAcceleration )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  const unsigned Drot = utilities::dimensions::convert<D>::rotation;
  //Check field existence
  static_assert(PARTICLETYPE::template providesNested<SURFACE,ANGLE>(), "Field SURFACE:ANGLE has to be provided");
  static_assert(PARTICLETYPE::template providesNested<MOBILITY,ANG_VELOCITY>(), "Field MOBILITY:ANG_VELOCITY has to be provided");
  //Calculate quantities
  Vector<T,Drot> angularVelocity( particle.getMobility().getAngularVelocity()
                                  + angularAcceleration * delTime );
  Vector<T,Drot> angle( particle.getSurface().getAngle() + angularVelocity * delTime );
  //Update values
  particle.template setField<SURFACE,ANGLE>( angle );
  particle.template setField<MOBILITY,ANG_VELOCITY>( angularVelocity );
}

template<typename T, typename PARTICLETYPE>
void eulerIntegration( Particle<T,PARTICLETYPE>& particle, T delTime,
                       Vector<T,PARTICLETYPE::d> acceleration,
                       Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angularAcceleration )
{
  //Calculate cartesian directions
  eulerIntegrationTranslation( particle, delTime, acceleration );
  //Calculate rotational directions
  if constexpr ( access::providesSurface<PARTICLETYPE>() ) {
    eulerIntegrationRotation( particle, delTime, angularAcceleration );
  }
}

//from FLUENT theory guide
//TODO: Check and include into unit tests
template<typename T, typename PARTICLETYPE>
void analyticalTranslation(
  Particle<T,PARTICLETYPE>& particle, Vector <T,PARTICLETYPE::d> acceleration, T delTime,
  Vector <T,PARTICLETYPE::d> fluid_vel )
{

  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  //Check field existence
  static_assert(PARTICLETYPE::template providesNested<GENERAL,POSITION>(), "Field POSITION has to be provided");
  static_assert(PARTICLETYPE::template providesNested<MOBILITY,VELOCITY>(), "Field VELOCITY has to be provided");
  //Calculate quant accelearation[0]
  if (acceleration[0]!=0){
    T tau_p0 = 1./(acceleration[0])*(fluid_vel[0]-particle.template getField<MOBILITY,VELOCITY>()[0]);
    Vector<T,D> velocity( fluid_vel+ std::exp((-delTime/tau_p0))*(particle.template getField<MOBILITY,VELOCITY>()- fluid_vel));
    Vector<T,D> position( particle.template getField<GENERAL,POSITION>() + delTime*(fluid_vel) + tau_p0*(1-  std::exp((-delTime/tau_p0)))*(particle.template getField<MOBILITY,VELOCITY>()-fluid_vel));
    //Update values
    particle.template setField<GENERAL,POSITION>( position );
    particle.template setField<MOBILITY,VELOCITY>( velocity );
  } else {
    Vector<T,D> position( particle.template getField<GENERAL,POSITION>() + delTime*particle.template getField<MOBILITY,VELOCITY>());
    //Update values
    particle.template setField<GENERAL,POSITION>( position );
    //particle.template setField<MOBILITY,VELOCITY>( velocity );
  }
}



} //namespace dynamics

} //namespace particles

} //namespace olb


#endif
