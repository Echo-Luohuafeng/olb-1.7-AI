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

#ifndef PARTICLE_DYNAMICS_BASE_HH
#define PARTICLE_DYNAMICS_BASE_HH


namespace olb {

namespace particles {

namespace dynamics {

template <typename T, typename PARTICLETYPE>
std::string& ParticleDynamics<T,PARTICLETYPE>::getName()
{
  return _name;
}

template <typename T, typename PARTICLETYPE>
std::string const& ParticleDynamics<T,PARTICLETYPE>::getName() const
{
  return _name;
}


template<typename T, typename PARTICLETYPE>
NoParticleDynamics<T,PARTICLETYPE>::NoParticleDynamics( T rhoDummy )
{
  this->getName() = "NoParticleDynamics";
}

template<typename T, typename PARTICLETYPE>
void NoParticleDynamics<T,PARTICLETYPE>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize)
{ }



template<typename T, typename PARTICLETYPE, typename PCONDITION>
VerletParticleDynamics<T,PARTICLETYPE,PCONDITION>::VerletParticleDynamics ( )
{
  this->getName() = "VerletParticleDynamics";
}

template<typename T, typename PARTICLETYPE, typename PCONDITION>
void VerletParticleDynamics<T,PARTICLETYPE,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  using namespace particles::access;
  //Check for particle condition
  doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,[&](){
    //Calculate acceleration
    auto acceleration = getAcceleration( particle );
    //Check for angular components
    if constexpr ( providesAngle<PARTICLETYPE>() ) {
      //Calculate angular acceleration
      auto angularAcceleration = getAngAcceleration( particle );
      //Verlet algorithm
      particles::dynamics::velocityVerletIntegration(
        particle, timeStepSize, acceleration, angularAcceleration );
      //Check if rotation matrix provided
      if constexpr ( providesRotationMatrix<PARTICLETYPE>() ) {
        //Update rotation matrix
        updateRotationMatrix( particle );
      }
    }
    else {
      //Verlet algorithm without rotation
      particles::dynamics::velocityVerletTranslation( particle,
        timeStepSize, timeStepSize*timeStepSize, acceleration );
    }
  }); //doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
VerletParticleDynamicsTranslationOnly<T,PARTICLETYPE,PCONDITION>::VerletParticleDynamicsTranslationOnly ( )
{
  this->getName() = "VerletParticleDynamicsTranslationOnly";
}

template<typename T, typename PARTICLETYPE, typename PCONDITION>
void VerletParticleDynamicsTranslationOnly<T,PARTICLETYPE,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  using namespace particles::access;
  //Check for particle condition
  doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,[&](){
    //Calculate acceleration
    auto acceleration = getAcceleration( particle );
    //Verlet algorithm without rotation
    particles::dynamics::velocityVerletTranslation(
      particle, timeStepSize, timeStepSize*timeStepSize, acceleration );
  }); //doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
VerletParticleDynamicsRotationOnly<T,PARTICLETYPE,PCONDITION>::VerletParticleDynamicsRotationOnly ( )
{
  this->getName() = "VerletParticleDynamicsRotationOnly";
}

template<typename T, typename PARTICLETYPE, typename PCONDITION>
void VerletParticleDynamicsRotationOnly<T,PARTICLETYPE,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  using namespace particles::access;
  //Check for particle condition
  doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,[&](){
    //Calculate angular acceleration
    auto angularAcceleration = getAngAcceleration( particle );
    //Verlet algorithm
    particles::dynamics::velocityVerletRotation(
      particle, timeStepSize, timeStepSize*timeStepSize, angularAcceleration );
    //Check if rotation matrix provided
    if constexpr ( providesRotationMatrix<PARTICLETYPE>() ) {
      //Update rotation matrix
      updateRotationMatrix( particle );
    }
  }); //doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>
}

template<typename T, typename PARTICLETYPE, typename PCONDITION>
VerletParticleDynamicsRotor<T,PARTICLETYPE,PCONDITION>::VerletParticleDynamicsRotor ( Vector<T,PARTICLETYPE::d> angVel )
{
  _angVel = angVel;
  this->getName() = "VerletParticleDynamicsRotor";
}

template<typename T, typename PARTICLETYPE, typename PCONDITION>
void VerletParticleDynamicsRotor<T,PARTICLETYPE,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  using namespace particles::access;
  //Check for particle condition
  doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,[&](){
    //Verlet algorithm
    particles::dynamics::velocityVerletRotor(
      particle, timeStepSize, _angVel );
    //Check if rotation matrix provided
    if constexpr ( providesRotationMatrix<PARTICLETYPE>() ) {
      //Update rotation matrix
      updateRotationMatrix( particle );
    }
  }); //doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>
}


template<typename T, typename PARTICLETYPE, bool useCubicBounds, typename PCONDITION>
VerletParticleDynamicsVelocityWallReflection<T,PARTICLETYPE,useCubicBounds,PCONDITION>::
  VerletParticleDynamicsVelocityWallReflection(
    SolidBoundary<T,PARTICLETYPE::d>& solidBoundary )
  : _solidBoundary(solidBoundary)
{
  this->getName() = "VerletParticleDynamicsVelocityWallReflection";
}


template<typename T, typename PARTICLETYPE, bool useCubicBounds, typename PCONDITION>
void VerletParticleDynamicsVelocityWallReflection<T,PARTICLETYPE,useCubicBounds,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  //Execute process of VerletParticleDynamcis
  VerletParticleDynamics<T,PARTICLETYPE,PCONDITION>::process(particle,timeStepSize);
  //Apply wall reflection
  boundaries::velocityWallReflection<useCubicBounds>(particle, _solidBoundary);
}


template<typename T, typename PARTICLETYPE, bool useCubicBounds, typename PCONDITION>
VerletParticleDynamicsWallCapture<T,PARTICLETYPE,useCubicBounds,PCONDITION>::
  VerletParticleDynamicsWallCapture(
    SolidBoundary<T,PARTICLETYPE::d>& solidBoundary )
  : _solidBoundary(solidBoundary)
{
  this->getName() = "VerletParticleDynamicsWallCapture";
}


template<typename T, typename PARTICLETYPE, bool useCubicBounds, typename PCONDITION>
void VerletParticleDynamicsWallCapture<T,PARTICLETYPE,useCubicBounds,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  //Execute process of VerletParticleDynamcis
  VerletParticleDynamics<T,PARTICLETYPE,PCONDITION>::process(particle,timeStepSize);
  //Apply wall capture
  boundaries::wallCapture<useCubicBounds>(particle, _solidBoundary);
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
VerletParticleDynamicsMaterialCapture<T,PARTICLETYPE,PCONDITION>::
  VerletParticleDynamicsMaterialCapture(
    std::shared_ptr<SuperIndicatorMaterial<T,PARTICLETYPE::d>> materialIndicator )
  : _materialIndicator(materialIndicator)
{
  this->getName() = "VerletParticleDynamicsMaterialCapture";
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
void VerletParticleDynamicsMaterialCapture<T,PARTICLETYPE,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  //Execute process of VerletParticleDynamcis
  VerletParticleDynamics<T,PARTICLETYPE,PCONDITION>::process(particle,timeStepSize);
  //Apply material capture
  boundaries::materialCapture(particle, *_materialIndicator);
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
VerletParticleDynamicsMaterialAwareWallCapture<T,PARTICLETYPE,PCONDITION>::
  VerletParticleDynamicsMaterialAwareWallCapture(
    SolidBoundary<T,PARTICLETYPE::d>& solidBoundary,
    std::shared_ptr<SuperIndicatorMaterial<T,PARTICLETYPE::d>> materialIndicator )
  : _solidBoundary(solidBoundary), _materialIndicator(materialIndicator)
{
  this->getName() = "VerletParticleDynamicsMaterialAwareWallCapture";
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
void VerletParticleDynamicsMaterialAwareWallCapture<T,PARTICLETYPE,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  //Execute process of VerletParticleDynamcis
  VerletParticleDynamics<T,PARTICLETYPE,PCONDITION>::process(particle,timeStepSize);
  //Apply wall capture
  boundaries::wallCaptureMaterialAware(particle, _solidBoundary, *_materialIndicator);
}


template<typename T, typename PARTICLETYPE, bool useCubicBounds, typename PCONDITION>
VerletParticleDynamicsEscape<T,PARTICLETYPE,useCubicBounds,PCONDITION>::
  VerletParticleDynamicsEscape(
    SolidBoundary<T,PARTICLETYPE::d>& solidBoundary )
  : _solidBoundary(solidBoundary)
{
  this->getName() = "VerletParticleDynamicsEscape";
}


template<typename T, typename PARTICLETYPE, bool useCubicBounds, typename PCONDITION>
void VerletParticleDynamicsEscape<T,PARTICLETYPE,useCubicBounds,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  //Execute process of VerletParticleDynamcis
  VerletParticleDynamics<T,PARTICLETYPE,PCONDITION>::process(particle,timeStepSize);
  //Apply escape boundary
  boundaries::escape<useCubicBounds>(particle, _solidBoundary);
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
VerletParticleDynamicsMaterialEscape<T,PARTICLETYPE,PCONDITION>::
  VerletParticleDynamicsMaterialEscape(
    std::shared_ptr<SuperIndicatorMaterial<T,PARTICLETYPE::d>> materialIndicator )
  : _materialIndicator(materialIndicator)
{
  this->getName() = "VerletParticleDynamicsMaterialEscape";
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
void VerletParticleDynamicsMaterialEscape<T,PARTICLETYPE,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  //Execute process of VerletParticleDynamcis
  VerletParticleDynamics<T,PARTICLETYPE,PCONDITION>::process(particle,timeStepSize);
  //Apply material escape
  boundaries::materialEscape(particle, *_materialIndicator);
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
VerletParticleDynamicsMaterialAwareEscape<T,PARTICLETYPE,PCONDITION>::
  VerletParticleDynamicsMaterialAwareEscape(
    SolidBoundary<T,PARTICLETYPE::d>& solidBoundary,
    std::shared_ptr<SuperIndicatorMaterial<T,PARTICLETYPE::d>> materialIndicator )
  : _solidBoundary(solidBoundary), _materialIndicator(materialIndicator)
{
  this->getName() = "VerletParticleDynamicsMaterialAwareEscape";
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
void VerletParticleDynamicsMaterialAwareEscape<T,PARTICLETYPE,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  //Execute process of VerletParticleDynamcis
  VerletParticleDynamics<T,PARTICLETYPE,PCONDITION>::process(particle,timeStepSize);
  //Apply escape boundary
  boundaries::escapeMaterialAware(particle, _solidBoundary, *_materialIndicator);
}



template<typename T, typename PARTICLETYPE, typename PCONDITION>
VerletParticleDynamicsMaterialCaptureAndEscape<T,PARTICLETYPE,PCONDITION>::
  VerletParticleDynamicsMaterialCaptureAndEscape(
    std::shared_ptr<SuperIndicatorMaterial<T,PARTICLETYPE::d>> captureMaterialIndicator,
    std::shared_ptr<SuperIndicatorMaterial<T,PARTICLETYPE::d>> escapeMaterialIndicator )
  : _captureMaterialIndicator(captureMaterialIndicator),
    _escapeMaterialIndicator(escapeMaterialIndicator)
{
  this->getName() = "VerletParticleDynamicsMaterialCaptureAndEscape";
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
void VerletParticleDynamicsMaterialCaptureAndEscape<T,PARTICLETYPE,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  //Execute process of VerletParticleDynamcis
  VerletParticleDynamics<T,PARTICLETYPE,PCONDITION>::process(particle,timeStepSize);
  //Apply material escape
  boundaries::materialCaptureAndEscape(particle, *_captureMaterialIndicator,
      *_escapeMaterialIndicator);
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
VerletParticleDynamicsMaterialAwareWallCaptureAndEscape<T,PARTICLETYPE,PCONDITION>::
  VerletParticleDynamicsMaterialAwareWallCaptureAndEscape(
    SolidBoundary<T,PARTICLETYPE::d>& solidBoundary,
    std::shared_ptr<SuperIndicatorMaterial<T,PARTICLETYPE::d>> captureMaterialIndicator,
    std::shared_ptr<SuperIndicatorMaterial<T,PARTICLETYPE::d>> escapeMaterialIndicator )
  : _solidBoundary(solidBoundary),
    _captureMaterialIndicator(captureMaterialIndicator),
    _escapeMaterialIndicator(escapeMaterialIndicator)
{
  this->getName() = "VerletParticleDynamicsMaterialAwareWallCaptureAndEscape";
}


template<typename T, typename PARTICLETYPE, typename PCONDITION>
void VerletParticleDynamicsMaterialAwareWallCaptureAndEscape<T,PARTICLETYPE,PCONDITION>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  //Execute process of VerletParticleDynamcis
  VerletParticleDynamics<T,PARTICLETYPE,PCONDITION>::process(particle,timeStepSize);
  //Apply escape boundary
  boundaries::wallCaptureAndEscapeMaterialAware(particle, _solidBoundary,
      *_captureMaterialIndicator, *_escapeMaterialIndicator);
}



template<typename T, typename PARTICLETYPE>
ParticleDetachmentDynamics<T,PARTICLETYPE>::ParticleDetachmentDynamics(
  SolidBoundary<T,PARTICLETYPE::d>& solidBoundary, Vector<T,PARTICLETYPE::d>& mainFlowDirection,
  T tiltThreshold )
  : _solidBoundary(solidBoundary), _mainFlowDirection(mainFlowDirection), _tiltThreshold(tiltThreshold)
{
  this->getName() = "ParticleDetachmentDynamics";
  static_assert(PARTICLETYPE::template providesNested<descriptors::SURFACE,descriptors::ANGLE>(),
                "Field SURFACE:ANGLE has to be provided");
  static_assert(PARTICLETYPE::template providesNested<descriptors::DYNBEHAVIOUR,descriptors::DETACHING>(),
                "Field DYNBEHAVIOUR:DETACHING has to be provided");
}

template<typename T, typename PARTICLETYPE>
void ParticleDetachmentDynamics<T,PARTICLETYPE>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  using namespace particles::access;
  using namespace descriptors;

  //Check for valid particles condition
  doWhenMeetingCondition<T,PARTICLETYPE,conditions::valid_particles>( particle,[&](){
    //Check current dynamic state
    bool detaching = isDetaching( particle );
    //State: Normal motion
    if (!detaching){
      //Execute process of VerletParticleDynamcis (valid, not detaching, active)
      VerletParticleDynamics<T,PARTICLETYPE,conditions::active_particles>
        ::process(particle,timeStepSize);
    }
    //State: Detaching
    else{
      //Check adhesion threshold (set active, if passed once)
      bool isAdhering = interaction::checkAdhesion( _solidBoundary, _mainFlowDirection, particle );
      //Perform detachment dynamics (valid, detaching, adhering)
      if (!isAdhering){
        //Calculate angular acceleration
        auto angularAcceleration = getAngAcceleration( particle );
        //Verlet algorithm
        particles::dynamics::velocityVerletRotation(
          particle, timeStepSize, timeStepSize*timeStepSize, angularAcceleration );
        //Check if rotation matrix provided and update
        if constexpr ( providesRotationMatrix<PARTICLETYPE>() ) {
          updateRotationMatrix( particle );
        }
        //Handle detachment from surface
        interaction::handleDetachment( _solidBoundary, _mainFlowDirection, particle );
        //Reevaluate state: Check, whether state has to be changed (e.g. detachment finished)
        interaction::evaluateDetachmentState( _solidBoundary, particle, _tiltThreshold );
      } //if (!isAdhering)
    }
  }); //doWhenMeetingCondition<T,PARTICLETYPE,conditions::valid_particles>
}





//TODO: REWORK FOLLOWING DYNAMICS ACCORDING TO OBOVE ONES



//TODO: remove for release
//TODO: rework according to CubicBoundsCheck
template<typename T, typename PARTICLETYPE>
VerletParticleDynamicsCubicBoundsAdhesion<T,PARTICLETYPE>::VerletParticleDynamicsCubicBoundsAdhesion (
  PhysR<T,PARTICLETYPE::d>& domainMin,
  PhysR<T,PARTICLETYPE::d>& domainMax
)
  : _domainMin(domainMin), _domainMax(domainMax)
{
  this->getName() = "VerletParticleDynamicsCubicBoundsAdhesion";
}


/**
 * - In order to move, both normal and at least one of the tangential
 *   adhesion force components have to be surpassed, thus
 *    adhesionSuperior = adhesionSuperior || bothTangentalSuperior.
 * - Carefull, when debugging: Dont't use << bool, but rather explicit
 *   if (bool){ << true } outputs.
 */

//TODO: remove for release
template<typename T, typename PARTICLETYPE>
void VerletParticleDynamicsCubicBoundsAdhesion<T,PARTICLETYPE>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  using namespace particles::access;
  //Calculate acceleration
  auto acceleration = getAcceleration( particle );
  //Note position before calculating movement (and store in dynStatePre)
  DynState<T,PARTICLETYPE> dynStatePre( getPosition(particle) );

  //Check for angular components
  if constexpr ( providesAngle<PARTICLETYPE>() ) {
    //Calculate angular acceleration
    auto angularAcceleration = getAngAcceleration( particle );
    //Note angle before calculating rotation (and store in dynStatePre)
    dynStatePre.angle = getAngle( particle );
    //Verlet algorithm
    particles::dynamics::velocityVerletIntegration(
      particle, timeStepSize, acceleration, angularAcceleration );
  } else {
    //Verlet algorithm
    particles::dynamics::velocityVerletIntegration(
      particle, timeStepSize, acceleration);
  }

  //Retrieve force and adhesion threshold
  auto force = getForce( particle );
  auto adhesionThreshold = getAdhesion( particle );

  //Check domain contact and potentially reset to positionPre and velocity zero
  bool adhesionSuperior = false;
  doAtCubicBoundPenetration( particle, _domainMin, _domainMax,
  [&](unsigned iDim, Vector<T,PARTICLETYPE::d>& normal, T distToBound ) {
    //Reset (only!) penetration direction
    resetDirection( particle, dynStatePre.position, iDim );
    //Check adhesion in normal direction
    //-only consider contributions pointing from the surface (normalForce > 0)
    // to ensure, that pressing on particles does not release them
    T normalForce = force[iDim]*normal[iDim];
    if (normalForce > 0 && normalForce < adhesionThreshold[0] ){
      adhesionSuperior = true;
    }
    //Check adhesion in all tangential directions
    // - info: checks in both directions in both positive and negative direction
    bool bothTangentalSuperior = true;
    for (int iDimT=1; iDimT<PARTICLETYPE::d; ++iDimT) {
      int jDim = (iDim+iDimT)%PARTICLETYPE::d;
      if (std::abs(force[jDim]) >= adhesionThreshold[1]){  //Tangential direction (absolute)
        bothTangentalSuperior = false;
      }
    }
    //Combine normal and tangential evaluation
    adhesionSuperior = adhesionSuperior || bothTangentalSuperior;
  });

  if constexpr ( providesAngle<PARTICLETYPE>() ) {
    //Check whether adhesion is superior
    if (adhesionSuperior){
      resetMovement( particle, dynStatePre.position, dynStatePre.angle );
    }
    //Update rotation matrix
    updateRotationMatrix( particle );
  } else {
    //Check whether adhesion is superior
    if (adhesionSuperior){ resetMovement( particle, dynStatePre.position); }
  }

}


//TODO: remove for release
template<typename T, typename PARTICLETYPE, typename DEPOSITION_MODEL>
VerletParticleDynamicsCubicBoundsDeposition<T,PARTICLETYPE,DEPOSITION_MODEL>
  ::VerletParticleDynamicsCubicBoundsDeposition (
    PhysR<T,PARTICLETYPE::d>& domainMin,
    PhysR<T,PARTICLETYPE::d>& domainMax,
    DEPOSITION_MODEL& depositionModel
)
  : _domainMin(domainMin), _domainMax(domainMax),
    _depositionModel( depositionModel )
{
  this->getName() = "VerletParticleDynamicsCubicBoundsDeposition";
}

//TODO: remove for release
template<typename T, typename PARTICLETYPE, typename DEPOSITION_MODEL>
void VerletParticleDynamicsCubicBoundsDeposition<T,PARTICLETYPE,DEPOSITION_MODEL>
  ::process (Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  using namespace particles::access;
  static_assert( providesActive<PARTICLETYPE>(), "Field ACTIVE has to be provided");

  //Check if active
  if (particle.template getField<descriptors::DYNBEHAVIOUR,descriptors::ACTIVE>()) {

    auto force = getForce( particle );
    auto velocity = getVelocity( particle );

    //Check domain contact and check deposition
    bool deposition = false;
    doAtCubicBoundPenetration( particle, _domainMin, _domainMax,
    [&](unsigned iDim, Vector<T,PARTICLETYPE::d>& normal, T distToBound ) {
      //Check deposition
      deposition = deposition || _depositionModel.checkDeposition(velocity,normal);
      //Wall contact treatment
      if (!deposition) {
        auto radius = getRadius( particle );
        T penetrationDepth = -distToBound;
        T velNormal = -normal[iDim]*velocity[iDim]; //Necessary for damping
        T forceNormal = _depositionModel.contactForceSphereHalfSpaceDampened(
                          radius, penetrationDepth, velNormal );
        force[iDim] = forceNormal*normal[iDim];
      }
    });
    particle.template setField<descriptors::FORCING,descriptors::FORCE>( force );

    //Run verlet or deactivate depending on deposition
    if (!deposition) {
      //Calculate acceleration
      auto acceleration = getAcceleration( particle );
      //Note position before calculating movement
      auto positionPre = getPosition( particle );
      //Check for angular components
      if constexpr ( providesAngle<PARTICLETYPE>() ) {
        //Calculate angular acceleration
        auto angularAcceleration = getAngAcceleration( particle );
        //Verlet algorithm
        particles::dynamics::velocityVerletIntegration>(
          particle, timeStepSize, acceleration, angularAcceleration );
        //Update rotation matrix
        updateRotationMatrix( particle );
      }
      else {
        //Verlet algorithm without angle
        particles::dynamics::velocityVerletIntegration>(
          particle, timeStepSize, acceleration );
      }
    }
    else {
      particle.template setField<descriptors::DYNBEHAVIOUR,descriptors::ACTIVE>( false );
      particle.template setField<descriptors::MOBILITY,descriptors::VELOCITY>( 0. );
      particle.template setField<descriptors::MOBILITY,descriptors::ACCELERATION_STRD>( 0. );
      if constexpr ( providesAngle<PARTICLETYPE>() ) {
        particle.template setField<descriptors::MOBILITY,descriptors::ANG_VELOCITY>( 0. );
        particle.template setField<descriptors::MOBILITY,descriptors::ANG_ACC_STRD>( 0. );
      }
    }
  }
}

} //namespace dynamics

} //namespace particles

} //namespace olb

#endif
