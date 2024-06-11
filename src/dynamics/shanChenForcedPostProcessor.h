/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert, Luiz Czelusniak,
 *                     Maximilian Schecher, Adrian Kummerlaender
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

#ifndef SHAN_CHEN_FORCED_POST_PROCESSOR_H
#define SHAN_CHEN_FORCED_POST_PROCESSOR_H

#define THIRD_COMPONENT
//#define FOURTH_COMPONENT

#include "core/blockStructure.h"
#include "core/postProcessing.h"

namespace olb {

/**
* Multiphysics class for coupling between different lattices.
*/

struct RhoStatistics  {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform
  {
    auto statistic = cell.template getField<descriptors::STATISTIC>();
    statistic[0] = cell.computeRho();
    cell.template setField<descriptors::STATISTIC>(statistic);
  }
};

// =========================================================================//
// ===========Shan Chen coupling without wall interaction===================//
// =========================================================================//

template <typename POTENTIAL>
struct ShanChenForcedPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct RHO0 : public descriptors::FIELD_BASE<2> { };
  struct G    : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_A : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_B : public descriptors::FIELD_BASE<1> { };

  using parameters = typename POTENTIAL::parameters::template decompose_into<
    meta::list<RHO0,G,OMEGA_A,OMEGA_B>::template include
  >;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::A>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::A>::descriptor_t;

    auto& cellA = cells.template get<names::A>();
    auto& cellB = cells.template get<names::B>();

    auto rho0 = parameters.template get<RHO0>();

    Vector<V,2> rhoField{
      cellA.template getFieldComponent<descriptors::STATISTIC>(0) * rho0[0],
      cellB.template getFieldComponent<descriptors::STATISTIC>(0) * rho0[1]
    };
    Vector<V,DESCRIPTOR::d> uA{};
    lbm<DESCRIPTOR>::computeJ(cellA, uA);
    Vector<V,DESCRIPTOR::d> uB{};
    lbm<DESCRIPTOR>::computeJ(cellB, uB);

    V omegaA = parameters.template get<OMEGA_A>();
    V omegaB = parameters.template get<OMEGA_B>();
    // Computation of the common velocity, shared among the two populations
    V rhoTot = rhoField[0]*omegaA
             + rhoField[1]*omegaB;

    Vector<V,DESCRIPTOR::d> uTot = (uA*rho0[0]*omegaA + uB*rho0[1]*omegaB) / rhoTot;

    // Computation of the interaction potential
    Vector<V,DESCRIPTOR::d> rhoBlockContribution{};
    Vector<V,DESCRIPTOR::d> rhoPartnerContribution{};
    V psi1 = POTENTIAL().compute(rhoField[0], parameters);
    V psi2 = POTENTIAL().compute(rhoField[1], parameters);

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      auto nextCellA = cellA.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      auto nextCellB = cellB.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      V rhoA = POTENTIAL().compute(nextCellA.template getFieldComponent<descriptors::STATISTIC>(0)*rho0[0], parameters);
      V rhoB = POTENTIAL().compute(nextCellB.template getFieldComponent<descriptors::STATISTIC>(0)*rho0[1], parameters);
      rhoBlockContribution   += psi2 * rhoA * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop);
      rhoPartnerContribution += psi1 * rhoB * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop);
    }

    // Computation and storage of the final velocity, consisting
    //   of u and the momentum difference due to interaction
    //   potential plus external force
    auto externalBlockForce   = cellA.template getField<descriptors::EXTERNAL_FORCE>();
    auto externalPartnerForce = cellB.template getField<descriptors::EXTERNAL_FORCE>();

    auto g = parameters.template get<G>();

    cellA.template setField<descriptors::VELOCITY>(uTot);
    cellB.template setField<descriptors::VELOCITY>(uTot);
    cellA.template setField<descriptors::FORCE>(externalBlockForce
        - g*rhoPartnerContribution/rhoField[0]);
    cellB.template setField<descriptors::FORCE>(externalPartnerForce
        - g*rhoBlockContribution/rhoField[1]);
  }

};


// =========================================================================//
// ==================Multi-component-multi-phase coupling===================//
// =========================================================================//

/** Multi-Component-Multi-Phase Shan-Chen force with thermodynamic equation
 * of state based on
 *
 * 1. Czelusniak L E, Mapelli V P, Guzella M S, Cabezas-Gómez L, Wagner A J (2020)
 * Force approach for the pseudopotential lattice Boltzmann method. Physical
 * Review E, 102(3), 033307.
 * DOI: 10.1103/PhysRevE.102.033307
 * 2. Peng C, Ayala L F, Ayala O M (2021) A thermodynamically consistent pseudo-
 * potential lattice Boltzmann model for multi-component, multiphase, partially
 * miscible mixtures. Journal of Computational Physics, 429, 110018.
 * DOI: 10.1016/j.jcp.2020.110018
 **/
template <typename POTENTIAL, unsigned N_COMPONENTS>
struct MCMPForcedPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct CHI         : public descriptors::FIELD_BASE<N_COMPONENTS> { };
  struct TEMPERATURE : public descriptors::FIELD_BASE<1> { };
  struct G           : public descriptors::FIELD_BASE<1> { };
  struct SIGMA       : public descriptors::FIELD_BASE<1> { };
  struct EPSILON     : public descriptors::FIELD_BASE<1> { };
  struct K           : public descriptors::FIELD_BASE<1> { };
  struct A           : public descriptors::FIELD_BASE<N_COMPONENTS> { };
  struct B           : public descriptors::FIELD_BASE<N_COMPONENTS> { };
  struct MM          : public descriptors::FIELD_BASE<N_COMPONENTS> { };
  struct TCRIT       : public descriptors::FIELD_BASE<N_COMPONENTS> { };
  struct DEVI        : public descriptors::FIELD_BASE<N_COMPONENTS> { };
  struct ALPHA       : public descriptors::FIELD_BASE<N_COMPONENTS*N_COMPONENTS> { };
  struct GI          : public descriptors::FIELD_BASE<N_COMPONENTS*N_COMPONENTS> { };
  struct GII         : public descriptors::FIELD_BASE<N_COMPONENTS*N_COMPONENTS> { };

  using parameters = meta::list<CHI,TEMPERATURE,G,SIGMA,EPSILON,K,A,B,MM,TCRIT,DEVI,ALPHA,GI,GII>;

  int getPriority() const {
    return 0;
  }

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::Component1>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::Component1>::descriptor_t;

    auto chi = parameters.template get<CHI>();
    auto t = parameters.template get<TEMPERATURE>();
    auto g = parameters.template get<G>();
    auto sig = parameters.template get<SIGMA>();
    auto eps = parameters.template get<EPSILON>();
    auto k = parameters.template get<K>();
    auto a = parameters.template get<A>();
    auto b = parameters.template get<B>();
    auto M = parameters.template get<MM>();
    auto T_c = parameters.template get<TCRIT>();
    auto m = parameters.template get<DEVI>();
    auto alpha = parameters.template get<ALPHA>();
    auto g_I = parameters.template get<GI>();
    auto g_II = parameters.template get<GII>();

    Vector<V,10> rhoField{};
    V rhoTot = 0;
    Vector<V,DESCRIPTOR::d> u_sum{};
    auto& cell1 = cells.template get<names::Component1>();
    rhoField[0] = cell1.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> u1{};
    lbm<DESCRIPTOR>::computeJ(cell1, u1);
    u_sum += u1;

    auto& cell2 = cells.template get<names::Component2>();
    rhoField[1] = cell2.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> u2{};
    lbm<DESCRIPTOR>::computeJ(cell2, u2);
    u_sum += u2;

#ifdef THIRD_COMPONENT
    auto& cell3 = cells.template get<names::Component3>();
    rhoField[2] = cell3.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> u3{};
    lbm<DESCRIPTOR>::computeJ(cell3, u3);
    u_sum += u3;
#endif
#ifdef FOURTH_COMPONENT
    auto& cell4 = cells.template get<names::Component4>();
    rhoField[3] = cell4.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> u4{};
    lbm<DESCRIPTOR>::computeJ(cell4, u4);
    u_sum += u4;
#endif
    for (unsigned n = 0; n < N_COMPONENTS; ++n) {
      rhoTot += rhoField[n];
    }

    // Computation of the interaction potential
    Vector<V,DESCRIPTOR::d> totalForce{};
    Vector<V,DESCRIPTOR::d> A_ij{};
    //V psi = POTENTIAL().compute(rhoField, t, k, a, b, T_c, m, alpha, g_I, g_II, M);
    V p = POTENTIAL().computeP(rhoField, t, a, b, T_c, m, alpha, g_I, g_II, M);
    V psi = util::sqrt(util::abs(6.*(k*p - rhoTot/3.)));
    //V psi = util::sqrt(6.*p);

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      Vector<V,10> rhoField_i{};
      auto nextCell1 = cell1.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      rhoField_i[0] = nextCell1.template getFieldComponent<descriptors::STATISTIC>(0);
      auto nextCell2 = cell2.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      rhoField_i[1] = nextCell2.template getFieldComponent<descriptors::STATISTIC>(0);
#ifdef THIRD_COMPONENT
      auto nextCell3 = cell3.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      rhoField_i[2] = nextCell3.template getFieldComponent<descriptors::STATISTIC>(0);
#endif
#ifdef FOURTH_COMPONENT
      auto nextCell4 = cell4.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      rhoField_i[3] = nextCell4.template getFieldComponent<descriptors::STATISTIC>(0);
#endif
      V rhoTot_i = 0;
      for (unsigned n = 0; n < N_COMPONENTS; ++n) {
        rhoTot_i += rhoField_i[n];
      }
      V p_i = POTENTIAL().computeP(rhoField_i, t, a, b, T_c, m, alpha, g_I, g_II, M);
      //V psi_i = POTENTIAL().compute(rhoField_i, t, k, a, b, T_c, m, alpha, g_I, g_II, M);
      V phi_i = 6.*(k*p_i - rhoTot_i/3.)/g;
      //V phi_i = 6.*p_i;
      V psi_i = 0;
      if (phi_i >= 0){
        psi_i = util::sqrt(util::abs(phi_i));
      }
      else {
        psi_i = -1*util::sqrt(util::abs(phi_i));
      }
      totalForce += -g * psi * psi_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop);
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        Vector<V,10> rhoField_j{};
        auto nextCell1 = cell1.neighbor(descriptors::c<DESCRIPTOR>(jPop));
        rhoField_j[0] = nextCell1.template getFieldComponent<descriptors::STATISTIC>(0);
        auto nextCell2 = cell2.neighbor(descriptors::c<DESCRIPTOR>(jPop));
        rhoField_j[1] = nextCell2.template getFieldComponent<descriptors::STATISTIC>(0);
#ifdef THIRD_COMPONENT
        auto nextCell3 = cell3.neighbor(descriptors::c<DESCRIPTOR>(jPop));
        rhoField_j[2] = nextCell3.template getFieldComponent<descriptors::STATISTIC>(0);
#endif
#ifdef FOURTH_COMPONENT
        auto nextCell4 = cell4.neighbor(descriptors::c<DESCRIPTOR>(jPop));
        rhoField_j[3] = nextCell4.template getFieldComponent<descriptors::STATISTIC>(0);
#endif
        V rhoTot_j = 0;
        for (unsigned n = 0; n < N_COMPONENTS; ++n) {
          rhoTot_j += rhoField_j[n];
        }
        V p_j = POTENTIAL().computeP(rhoField_j, t, a, b, T_c, m, alpha, g_I, g_II, M);
        //V psi_j = POTENTIAL().compute(rhoField_j, t, k, a, b, T_c, m, alpha, g_I, g_II, M);
        V phi_j = 6.*(k*p_j - (rhoTot_j)/3.)/g;
        //V phi_j = 6.*p_j;
        V psi_j = 0;
        if (phi_j >= 0){
          psi_j = util::sqrt(util::abs(phi_j));
        }
        else {
          psi_j = -1*util::sqrt(util::abs(phi_j));
        }
        Vector<V,DESCRIPTOR::d> c1{};
        Vector<V,DESCRIPTOR::d> c2{};
        for (int m = 0; m < DESCRIPTOR::d; ++m) {
          Vector<V,DESCRIPTOR::d> C1{};
          V C2 = 0.;
          for (int n = 0; n < DESCRIPTOR::d; ++n) {
            C1[n] += descriptors::c<DESCRIPTOR>(jPop)[m] * descriptors::c<DESCRIPTOR>(jPop)[n];
            C2 += descriptors::c<DESCRIPTOR>(jPop)[n] * descriptors::c<DESCRIPTOR>(jPop)[n] - 1./descriptors::invCs2<V,DESCRIPTOR>();
            if (n==m){
              C1[n] -= 1./descriptors::invCs2<V,DESCRIPTOR>();
            }
          }
          for (int n = 0; n < DESCRIPTOR::d; ++n) {
            c1[m] += descriptors::c<DESCRIPTOR>(iPop)[n] * C1[n];
          }
          c2[m] += descriptors::c<DESCRIPTOR>(iPop)[m] * C2;
        }
        A_ij = 0.5*g*descriptors::t<V,DESCRIPTOR>(iPop)*descriptors::t<V,DESCRIPTOR>(jPop)*descriptors::invCs2<V,DESCRIPTOR>()*
                 ((3./2.*eps-(sig-1))*c1 + (sig-1)*c2);
        totalForce += A_ij * psi_i * psi_j;
      }

    }
    // Computation of the common velocity, shared among the two populations, consisting
    // of u and the momentum difference due to interaction potential plus external force
    auto externalForce1 = cell1.template getField<descriptors::EXTERNAL_FORCE>();
    Vector<V,DESCRIPTOR::d> uTot = (u_sum + totalForce/2.) / rhoTot + externalForce1/2.;
    cell1.template setField<descriptors::FORCE>(externalForce1 + chi[0]*totalForce/rhoField[0]);
    cell1.template setField<descriptors::VELOCITY>(uTot);

    auto externalForce2 = cell2.template getField<descriptors::EXTERNAL_FORCE>();
    cell2.template setField<descriptors::FORCE>(externalForce2 + chi[1]*totalForce/rhoField[1]);
    cell2.template setField<descriptors::VELOCITY>(uTot);

#ifdef THIRD_COMPONENT
    auto externalForce3 = cell3.template getField<descriptors::EXTERNAL_FORCE>();
    cell3.template setField<descriptors::FORCE>(externalForce3 + chi[2]*totalForce/rhoField[2]);
    cell3.template setField<descriptors::VELOCITY>(uTot);
#endif
#ifdef FOURTH_COMPONENT
    auto externalForce4 = cell4.template getField<descriptors::EXTERNAL_FORCE>();
    cell4.template setField<descriptors::FORCE>(externalForce4 + chi[3]*totalForce/rhoField[3]);
    cell4.template setField<descriptors::VELOCITY>(uTot);
#endif

    cell1.template setField<descriptors::SCALAR>(p);
    cell2.template setField<descriptors::SCALAR>(psi);
  }

};

}

#endif
