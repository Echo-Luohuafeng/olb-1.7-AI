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

#ifndef ADSORPTION_COUPLING_POST_PROCESSOR_3D_H
#define ADSORPTION_COUPLING_POST_PROCESSOR_3D_H

#include "core/blockStructure.h"
#include "core/postProcessing.h"
#include "core/util.h"
#include "dynamics/latticeDescriptors.h"
#include "utilities/omath.h"
#include "../dynamics/advectionDiffusionForces.hh"
#include "../dynamics/advectionDiffusionForces.h"

#include "adsorptionReaction.h"

namespace olb {

/** Coupling post processor for adsorption on moving particles.
 *
 * @tparam T
 * @tparam NSDESCRIPTOR solute lattice
 * @tparam ADEDESCRIPTOR particle lattice
 */
  template<typename T, typename NSDESCRIPTOR, typename ADEDESCRIPTOR>
  class AdsorptionCouplingPostProcessor3D : public LocalPostProcessor3D<T, NSDESCRIPTOR> {
  public:
      /** Coupling post processor for adsorption on moving particles.
       *
       * Coupling between particle lattice, solute lattice and loading lattice.
       * Add coupling to fluid lattice with NSLattice.addLatticeCoupling() with the partner lattices in that order.
       *
       * @param partners_
       * @param k
       * @param isotherm_ lambda function that provides a relationship between particle loading and concentration
       */
      AdsorptionCouplingPostProcessor3D(
              int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
              std::vector<BlockStructureD<3> *> partners_,
              AdsorptionReaction<T, ADEDESCRIPTOR> *adsorptionReaction_);

      int extent() const override {
        return 0;
      }

      int extent(int whichDirection) const override {
        return 0;
      }

      void process(BlockLattice<T, NSDESCRIPTOR> &blockLattice) override;

      void processSubDomain(BlockLattice<T, NSDESCRIPTOR> &blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);

  private:
      int x0, x1, y0, y1, z0, z1;
      BlockLattice<T, ADEDESCRIPTOR> *particleLattice;
      BlockLattice<T, ADEDESCRIPTOR> *soluteLattice;
      BlockLattice<T, ADEDESCRIPTOR> *loadingLattice;
      AdsorptionReaction<T, ADEDESCRIPTOR> *adsorptionReaction;
  };

/**
 * Generates post processor AdsorptionCouplingPostProcessor3D
 *
 * @tparam T
 * @tparam NSDESCRIPTOR
 * @tparam ADEDESCRIPTOR
 *
 * @see AdsorptionCouplingPostProcessor3D
 */
  template<typename T, typename NSDESCRIPTOR, typename ADEDESCRIPTOR>
  class AdsorptionCouplingPostProcessorGenerator3D
          : public LatticeCouplingGenerator3D<T, NSDESCRIPTOR> {
  public:
      /** Generate adsorption reaction post processor
       *
       * @param reaction class implementing AdsorptionReaction with the reaction and isotherm data
       */
      explicit AdsorptionCouplingPostProcessorGenerator3D(AdsorptionReaction<T, ADEDESCRIPTOR> *reaction,
                                                          int x0=0, int x1=0, int y0=0, int y1=0, int z0=0, int z1=0);

      PostProcessor3D<T, NSDESCRIPTOR> *generate(std::vector<BlockStructureD<3> *> partners) const;

      LatticeCouplingGenerator3D<T, NSDESCRIPTOR> *clone() const override;

  private:
      AdsorptionReaction<T, ADEDESCRIPTOR> *reaction;
  };

/// Coupler for solute that is only coupled to the fluid velocity.
///
/// \tparam T
/// \tparam NSDESCRIPTOR
/// \tparam CADDESCRIPTOR
template<typename T, typename NSDESCRIPTOR, typename CADDESCRIPTOR>
class PassiveSoluteCouplingPostProcessor3D : public LocalPostProcessor3D<T, NSDESCRIPTOR> {
private:
    bool tick = true;

public:
    PassiveSoluteCouplingPostProcessor3D(
            int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
            std::vector<BlockStructureD<3> *> partners_);

    int extent() const override {
      return 0;
    }

    int extent(int whichDirection) const override {
      return 0;
    }

    void process(BlockLattice<T, NSDESCRIPTOR> &blockLattice) override;

    void processSubDomain(BlockLattice<T, NSDESCRIPTOR> &blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);

private:
    int x0, x1, y0, y1, z0, z1;
    BlockLattice<T, CADDESCRIPTOR> *partner;
};

template<typename T, typename NSDESCRIPTOR, typename CADDESCRIPTOR>
class PassiveSoluteCouplingPostProcessorGenerator3D
        : public LatticeCouplingGenerator3D<T, NSDESCRIPTOR> {
public:
    PassiveSoluteCouplingPostProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);

    PostProcessor3D<T, NSDESCRIPTOR> *generate(std::vector<BlockStructureD<3> *> partners) const override;

    LatticeCouplingGenerator3D<T, NSDESCRIPTOR> *clone() const override;

private:
    T k;
};


template<typename T, typename NSDESCRIPTOR, typename ADEDESCRIPTOR>
class AdsorptionFullCouplingPostProcessor3D : public LocalPostProcessor3D<T, NSDESCRIPTOR> {
public:
      AdsorptionFullCouplingPostProcessor3D(
              int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int iC,
              std::vector<BlockStructureD<3> *> partners_,
              AdsorptionReaction<T, ADEDESCRIPTOR> *adsorptionReaction_,
              std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T, NSDESCRIPTOR,ADEDESCRIPTOR> > > forces_);

      int extent() const override {
        return 0;
      }

      int extent(int whichDirection) const override {
        return 0;
      }
      void process(BlockLattice<T, NSDESCRIPTOR> &blockLattice) override;
      void processSubDomain(BlockLattice<T, NSDESCRIPTOR> &blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);

protected:
      std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T, NSDESCRIPTOR, ADEDESCRIPTOR> > > forces;

private:
      int x0, x1, y0, y1, z0, z1, iC;
      BlockLattice<T, ADEDESCRIPTOR> *particleLattice;
      BlockLattice<T, ADEDESCRIPTOR> *soluteLattice;
      BlockLattice<T, ADEDESCRIPTOR> *loadingLattice;
      AdsorptionReaction<T, ADEDESCRIPTOR> *adsorptionReaction;
  };

template<typename T, typename NSDESCRIPTOR, typename ADEDESCRIPTOR>
class AdsorptionFullCouplingPostProcessorGenerator3D
          : public LatticeCouplingGenerator3D<T, NSDESCRIPTOR> {
public:
      explicit AdsorptionFullCouplingPostProcessorGenerator3D(AdsorptionReaction<T, ADEDESCRIPTOR> *reaction,
                                                          int x0=0, int x1=0, int y0=0, int y1=0, int z0=0, int z1=0);

      PostProcessor3D<T, NSDESCRIPTOR> *generate(std::vector<BlockStructureD<3> *> partners) const;

      LatticeCouplingGenerator3D<T, NSDESCRIPTOR> *clone() const override;
      void addForce(AdvectionDiffusionForce3D<T, NSDESCRIPTOR, ADEDESCRIPTOR> &force);

private:
      AdsorptionReaction<T, ADEDESCRIPTOR> *reaction;

protected:
  std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T, NSDESCRIPTOR, ADEDESCRIPTOR> > > ADforces;
  };
}
#endif //OLB_ADSORPTIONCOUPLING3D_H
