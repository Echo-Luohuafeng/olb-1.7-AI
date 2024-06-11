/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2022 Florian Raichle
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

#ifndef ADSORPTION_COUPLING_POST_PROCESSOR_3D_HH
#define ADSORPTION_COUPLING_POST_PROCESSOR_3D_HH

#include "adsorptionCoupling3D.h"
#include "../dynamics/advectionDiffusionForces.hh"
#include "../dynamics/advectionDiffusionForces.h"
namespace olb {
///All in one adsorption coupling
///
template<typename T, typename NSDESCRIPTOR, typename ADEDESCRIPTOR>
AdsorptionFullCouplingPostProcessor3D<T, NSDESCRIPTOR, ADEDESCRIPTOR>::AdsorptionFullCouplingPostProcessor3D(
      int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int iC_,
      std::vector<BlockStructureD<3> *> partners_,
      AdsorptionReaction<T, ADEDESCRIPTOR> *adsorptionReaction_,
      std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T,NSDESCRIPTOR,ADEDESCRIPTOR>>> forces_)
      : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), iC(iC_),
        particleLattice(static_cast<BlockLattice<T, ADEDESCRIPTOR> *>(partners_[0])),
        soluteLattice(static_cast<BlockLattice<T, ADEDESCRIPTOR> *>(partners_[1])),
        loadingLattice(static_cast<BlockLattice<T, ADEDESCRIPTOR> *>(partners_[2])),
        adsorptionReaction(adsorptionReaction_),
        forces(forces_) {
    this->getName() = "AdsorptionFullCouplingPostProcessor3D";
  }

template<typename T, typename NSDESCRIPTOR, typename ADEDESCRIPTOR>
void AdsorptionFullCouplingPostProcessor3D<T, NSDESCRIPTOR, ADEDESCRIPTOR>::process(BlockLattice<T, NSDESCRIPTOR> &blockLattice) {
    processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
  }

template<typename T, typename NSDESCRIPTOR, typename ADEDESCRIPTOR>
void AdsorptionFullCouplingPostProcessor3D<T, NSDESCRIPTOR, ADEDESCRIPTOR>::processSubDomain(BlockLattice<T, NSDESCRIPTOR> &blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) {

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if (util::intersect(
      x0, x1, y0, y1, z0, z1,
      x0_, x1_, y0_, y1_, z0_, z1_,
      newX0, newX1, newY0, newY1, newZ0, newZ1)) {
    for (int iX = newX0; iX <= newX1; ++iX) {
      for (int iY = newY0; iY <= newY1; ++iY) {
        for (int iZ = newZ0; iZ <= newZ1; ++iZ) {
          //computation of particle velocity for particle and load lattices
          auto vel = blockLattice.get(iX, iY, iZ).template getField<descriptors::VELOCITY>();
          T forceValue[3] = {0.,0.,0.};
          if (forces.begin() != forces.end()) {
            auto velXp = blockLattice.get(iX+1, iY, iZ).template getField<descriptors::VELOCITY>();
            auto velXn = blockLattice.get(iX-1, iY, iZ).template getField<descriptors::VELOCITY>();
            auto velYp = blockLattice.get(iX, iY+1, iZ).template getField<descriptors::VELOCITY>();
            auto velYn = blockLattice.get(iX, iY-1, iZ).template getField<descriptors::VELOCITY>();
            auto velZp = blockLattice.get(iX, iY, iZ+1).template getField<descriptors::VELOCITY>();
            auto velZn = blockLattice.get(iX, iY, iZ-1).template getField<descriptors::VELOCITY>();
            T velGrad[3] = {0., 0., 0.};
            velGrad[0] = 0.5*(vel[0]*(velXp[0] - velXn[0]) + vel[1]*(velYp[0] - velYn[0]) + vel[2]*(velZp[0] - velZn[0]));
            velGrad[1] = 0.5*(vel[0]*(velXp[1] - velXn[1]) + vel[1]*(velYp[1] - velYn[1]) + vel[2]*(velZp[1] - velZn[1]));
            velGrad[2] = 0.5*(vel[0]*(velXp[2] - velXn[2]) + vel[1]*(velYp[2] - velYn[2]) + vel[2]*(velZp[2] - velZn[2]));

            int latticeR[4] = {iC, iX, iY, iZ};
            auto nsCell = blockLattice.get(iX, iY, iZ);
            auto adCell = particleLattice->get(iX, iY, iZ);
            for (AdvectionDiffusionForce3D<T, NSDESCRIPTOR, ADEDESCRIPTOR>& f : forces) {
              f.applyForce(forceValue, &nsCell, &adCell, vel.data(), latticeR);
              }

            // compute new particle velocity under action of forces
            Vector<T,ADEDESCRIPTOR::d> newVel;
            for (int i=0; i < ADEDESCRIPTOR::d; i++) {
              newVel[i] = vel[i] + forceValue[i] - velGrad[i];
              }
            particleLattice->get(iX, iY, iZ).template setField<descriptors::VELOCITY>(newVel);
            loadingLattice->get(iX, iY, iZ).template setField<descriptors::VELOCITY>(newVel);
          }
          else {   // set particle velocity to the carrier fluid velocity
            particleLattice->get(iX, iY, iZ).template setField<descriptors::VELOCITY>(vel);
            loadingLattice->get(iX, iY, iZ).template setField<descriptors::VELOCITY>(vel);
            }

          //setting of the solute velocity to the carrier fluid velocity
          soluteLattice->get(iX, iY, iZ).template setField<descriptors::VELOCITY>(vel);

          //computation of the adsorption source terms for solute and load lattices
          T soluteConcentration = soluteLattice->get(iX, iY, iZ).computeRho();
          T particleConcentration = particleLattice->get(iX, iY, iZ).computeRho();
          T particleLoading = loadingLattice->get(iX, iY, iZ).computeRho();
          Vector<T, 2> reactionRates = adsorptionReaction->getReactionRate(soluteConcentration, particleLoading, particleConcentration);
          loadingLattice->get(iX, iY, iZ).template setField<descriptors::SOURCE>(reactionRates[1]);
          soluteLattice->get(iX, iY, iZ).template setField<descriptors::SOURCE>(-reactionRates[0]);
        }
      }
    }
  }
}

template<typename T, typename NSDESCRIPTOR, typename ADEDESCRIPTOR>
AdsorptionFullCouplingPostProcessorGenerator3D<T, NSDESCRIPTOR, ADEDESCRIPTOR>::
  AdsorptionFullCouplingPostProcessorGenerator3D(AdsorptionReaction<T, ADEDESCRIPTOR> *reaction,
                                             int x0, int x1, int y0, int y1, int z0, int z1)
  : LatticeCouplingGenerator3D<T, NSDESCRIPTOR>(x0, x1, y0, y1, z0, z1), reaction(reaction)
  {}
template<typename T, typename NSDESCRIPTOR, typename ADEDESCRIPTOR>
PostProcessor3D<T, NSDESCRIPTOR>* AdsorptionFullCouplingPostProcessorGenerator3D<T, NSDESCRIPTOR, ADEDESCRIPTOR>::generate(
    std::vector<BlockStructureD<3> *> partners) const {
  return new AdsorptionFullCouplingPostProcessor3D<T, NSDESCRIPTOR, ADEDESCRIPTOR>(
      this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, this->iC, partners, reaction, ADforces);
}

template<typename T, typename NSDESCRIPTOR, typename ADEDESCRIPTOR>
LatticeCouplingGenerator3D<T, NSDESCRIPTOR>* AdsorptionFullCouplingPostProcessorGenerator3D<T, NSDESCRIPTOR, ADEDESCRIPTOR>::clone() const {
  return new AdsorptionFullCouplingPostProcessorGenerator3D<T, NSDESCRIPTOR, ADEDESCRIPTOR>(*this);
}

template<typename T, typename NSDESCRIPTOR, typename ADEDESCRIPTOR>
void AdsorptionFullCouplingPostProcessorGenerator3D<T, NSDESCRIPTOR, ADEDESCRIPTOR>::addForce(
  AdvectionDiffusionForce3D<T,NSDESCRIPTOR,ADEDESCRIPTOR> &force)
{
  ADforces.push_back(force);
}
}
#endif
