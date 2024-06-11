/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_H
#define NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_H

#include "core/operator.h"

namespace olb {

/// TotalEnthalpyPhaseChange between a Navier-Stokes and an Advection-Diffusion lattice
struct TotalEnthalpyPhaseChangeCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };
  struct T0 : public descriptors::FIELD_BASE<1> { };

  struct CP_S : public descriptors::FIELD_BASE<1> { };
  struct CP_L : public descriptors::FIELD_BASE<1> { };
  struct T_S  : public descriptors::FIELD_BASE<1> { };
  struct T_L  : public descriptors::FIELD_BASE<1> { };
  struct L    : public descriptors::FIELD_BASE<1> { };
  // V H_s  = cp_s * T_s;
  // V H_l  = cp_l * T_l + l;


  using parameters = meta::list<FORCE_PREFACTOR,T0, CP_S, CP_L, T_S, T_L, L>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Temperature>();

    V enthalpy = cellADE.computeRho();
    auto& dynamics = cellADE.getDynamics();
    /*auto& dynParams = static_cast<ParametersOfDynamicsD<DYNAMICS>&>(
      tPartner->template getData<OperatorParameters<DYNAMICS>>());*/

    //cellNSE.template setField<descriptors::POROSITY>(
      //dynamics->template computeLiquidFraction<T>(parameters, enthalpy)
    //);

    auto temperature = cellADE.template getFieldPointer<descriptors::TEMPERATURE>();
    //temperature[0] = dynamics->computeTemperature<V, parameters, enthalpy>();

    // computation of the bousinessq force
    auto force = cellNSE.template getFieldPointer<descriptors::FORCE>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();
    V temperatureDifference = cellADE.computeRho() - parameters.template get<T0>();
    for (unsigned iD = 0; iD < DESCRIPTOR::d; ++iD) {
      force[iD] = forcePrefactor[iD] * temperatureDifference;
    }
    // Velocity coupling
    V u[DESCRIPTOR::d] { };
    cellNSE.computeU(u);
    cellADE.template setField<descriptors::VELOCITY>(u);
  }

};



/// Coupling between a Navier-Stokes and an Advection-Diffusion lattice
struct NavierStokesAdvectionDiffusionCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };
  struct T0 : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<FORCE_PREFACTOR,T0>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Temperature>();

    // computation of the bousinessq force
    auto force = cellNSE.template getFieldPointer<descriptors::FORCE>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();
    V temperatureDifference = cellADE.computeRho() - parameters.template get<T0>();
    for (unsigned iD = 0; iD < DESCRIPTOR::d; ++iD) {
      force[iD] = forcePrefactor[iD] * temperatureDifference;
    }

    // Velocity coupling
    V u[DESCRIPTOR::d] { };
    cellNSE.computeU(u);
    cellADE.template setField<descriptors::VELOCITY>(u);
  }
};

/// Velocity coupling between Navier-Stokes and an Advection-Diffusion lattice
struct NavierStokesAdvectionDiffusionVelocityCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Concentration0>();

    // Velocity coupling
    V u[DESCRIPTOR::d] { };
    cellNSE.computeU(u);
    cellADE.template setField<descriptors::VELOCITY>(u);
  }
};

/// AD coupling with Boussinesq bouancy for Smagorinsky-LES
struct SmagorinskyBoussinesqCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };
  struct SMAGORINSKY_PREFACTOR : public descriptors::FIELD_BASE<1> { };
  struct PR_TURB : public descriptors::FIELD_BASE<1> { };
  struct T0 : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_NSE : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_ADE : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<FORCE_PREFACTOR,SMAGORINSKY_PREFACTOR,PR_TURB,T0,OMEGA_NSE,OMEGA_ADE>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    using DESCRIPTOR_ADE = typename CELLS::template value_t<names::Temperature>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Temperature>();

    // computation of the bousinessq force
    auto force = cellNSE.template getFieldPointer<descriptors::FORCE>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();
    V temperatureDifference = cellADE.computeRho() - parameters.template get<T0>();

    for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
      force[iD] = forcePrefactor[iD] * temperatureDifference;
    }

    // Velocity coupling
    auto u = cellADE.template getField<descriptors::VELOCITY>();
    // tau coupling
    auto tauNS = cellNSE.template getFieldPointer<descriptors::TAU_EFF>();
    auto tauAD = cellADE.template getFieldPointer<descriptors::TAU_EFF>();

    V rho, pi[util::TensorVal<DESCRIPTOR>::n] { };
    cellNSE.computeAllMomenta(rho, u.data(), pi);
    cellADE.template setField<descriptors::VELOCITY>(u);
    V PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
    }
    V PiNeqNorm    = util::sqrt(PiNeqNormSqr);
    /// Molecular realaxation time
    V tau_mol_NS = V{1} / parameters.template get<OMEGA_NSE>();
    V tau_mol_AD = V{1} / parameters.template get<OMEGA_ADE>();
    /// Turbulent realaxation time
    V smagoPrefactor = parameters.template get<SMAGORINSKY_PREFACTOR>();
    V tau_turb_NS = V{0.5}*(util::sqrt(tau_mol_NS*tau_mol_NS + smagoPrefactor/rho*PiNeqNorm) - tau_mol_NS);
    /// Effective realaxation time
    tauNS[0] = tau_mol_NS+tau_turb_NS;

    V prTurb = parameters.template get<PR_TURB>();
    V tauTurbADPrefactor = descriptors::invCs2<V,DESCRIPTOR_ADE>() / descriptors::invCs2<V,DESCRIPTOR>() / prTurb;
    V tau_turb_AD = tau_turb_NS * tauTurbADPrefactor;
    tauAD[0] = tau_mol_AD+tau_turb_AD;
  }

};


/// Reaction Coupling for the In-Bulk appraoch of lognitudinalMixing3d example
template<typename T>
struct LongitudinalMixingReactionCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct REACTION_CONSTANT : public descriptors::FIELD_BASE<1> { };
  struct EQUILIBRIUM : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<REACTION_CONSTANT, EQUILIBRIUM>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
        T forwardConstant = parameters.template get<REACTION_CONSTANT>();
        T Ceq = parameters.template get<EQUILIBRIUM>();

        auto allow_source = cells.template get<names::Concentration0>().template getField<descriptors::GAMMA>(); //allow source only at surface
        T concC = cells.template get<names::Concentration0>().computeRho(); //current concentration
        T source = -forwardConstant*(concC - Ceq); //difference between current Conc and equilibrium concentration determines source
        cells.template get<names::Concentration0>().template setField<descriptors::SOURCE>(source*allow_source);
  }
};



/// LES-ADE coupling for multiple reactions
template<typename T, int numComp>
struct LESReactionCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct LATTICE_REACTION_COEFF : public descriptors::FIELD_BASE<1> { };
  struct STOCH_COEFF : public descriptors::FIELD_BASE<numComp> { };
  struct REACTION_ORDER : public descriptors::FIELD_BASE<numComp> { };
  struct SMAGORINSKY_PREFACTOR : public descriptors::FIELD_BASE<1> { };
  struct SCHMIDT : descriptors::FIELD_BASE<numComp> { };
  struct OMEGA_NSE : public descriptors::FIELD_BASE<1> { };
  struct OMEGAS_ADE : public descriptors::FIELD_BASE<numComp> { };

  using parameters = meta::list<LATTICE_REACTION_COEFF, STOCH_COEFF, REACTION_ORDER, SMAGORINSKY_PREFACTOR, SCHMIDT, OMEGA_NSE, OMEGAS_ADE>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    using DESCRIPTOR_ADE = typename CELLS::template value_t<names::Concentration0>::descriptor_t;
    // Velocity coupling
    auto u = cells.template get<names::Concentration0>().template getField<descriptors::VELOCITY>();
    T rho, pi[util::TensorVal<DESCRIPTOR>::n] { };
    cells.template get<names::NavierStokes>().computeAllMomenta(rho, u.data(), pi);
    cells.template get<names::Concentration0>().template setField<descriptors::VELOCITY>(u);
    cells.template get<names::Concentration1>().template setField<descriptors::VELOCITY>(u);
    cells.template get<names::Concentration2>().template setField<descriptors::VELOCITY>(u);
    // Stress tensor
    T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
    }
    T PiNeqNorm    = util::sqrt(PiNeqNormSqr);
    /// Molecular realaxation time
    T tau_mol_NS = T{1} / parameters.template get<OMEGA_NSE>();
    auto omegasAD = parameters.template get<OMEGAS_ADE>();
    T tau_mol_AD0 = T{1} / omegasAD[0];
    T tau_mol_AD1 = T{1} / omegasAD[1];
    T tau_mol_AD2 = T{1} / omegasAD[2];
    /// Turbulent realaxation time
    T smagoPrefactor = parameters.template get<SMAGORINSKY_PREFACTOR>();
    T tau_turb_NS = T{0.5}*(util::sqrt(tau_mol_NS*tau_mol_NS + smagoPrefactor/rho*PiNeqNorm) - tau_mol_NS);
    // Schmidt number stabilization
    auto Sc = parameters.template get<SCHMIDT>();
    T tauTurbADPrefactor0 = descriptors::invCs2<T,DESCRIPTOR_ADE>() / descriptors::invCs2<T,DESCRIPTOR>() / Sc[0];
    T tauTurbADPrefactor1 = descriptors::invCs2<T,DESCRIPTOR_ADE>() / descriptors::invCs2<T,DESCRIPTOR>() / Sc[1];
    T tauTurbADPrefactor2 = descriptors::invCs2<T,DESCRIPTOR_ADE>() / descriptors::invCs2<T,DESCRIPTOR>() / Sc[2];
    T tau_turb_AD0 = tau_turb_NS * tauTurbADPrefactor0;
    T tau_turb_AD1 = tau_turb_NS * tauTurbADPrefactor1;
    T tau_turb_AD2 = tau_turb_NS * tauTurbADPrefactor2;
    cells.template get<names::Concentration0>().template setField<descriptors::OMEGA>(T{1} / (tau_mol_AD0 + tau_turb_AD0));
    cells.template get<names::Concentration1>().template setField<descriptors::OMEGA>(T{1} / (tau_mol_AD1 + tau_turb_AD1));
    cells.template get<names::Concentration2>().template setField<descriptors::OMEGA>(T{1} / (tau_mol_AD2 + tau_turb_AD2));

    auto stochCoeff = parameters.template get<STOCH_COEFF>();
    auto reactionOrder = parameters.template get<REACTION_ORDER>();
    T source = parameters.template get<LATTICE_REACTION_COEFF>();
    {
      T conc = cells.template get<names::Concentration0>().computeRho();
      source *= util::pow(conc, reactionOrder[0]);
    }
    {
      T conc = cells.template get<names::Concentration1>().computeRho();
      source *= util::pow(conc, reactionOrder[1]);
    }
    {
      T conc = cells.template get<names::Concentration2>().computeRho();
      source *= util::pow(conc, reactionOrder[2]);
    }
    {
      cells.template get<names::Concentration0>().template setField<descriptors::SOURCE>(stochCoeff[0]*source);
      cells.template get<names::Concentration1>().template setField<descriptors::SOURCE>(stochCoeff[1]*source);
      cells.template get<names::Concentration2>().template setField<descriptors::SOURCE>(stochCoeff[2]*source);
    }
  }

};

/// LES-ADE coupling with Schmidt number stabilization
template<typename T>
struct LESADECoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct SMAGORINSKY_PREFACTOR : public descriptors::FIELD_BASE<1> { };
  struct SCHMIDT : descriptors::FIELD_BASE<1> { };
  struct OMEGA_NSE : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_ADE : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<SMAGORINSKY_PREFACTOR, SCHMIDT, OMEGA_NSE, OMEGA_ADE>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    using DESCRIPTOR_ADE = typename CELLS::template value_t<names::Concentration0>::descriptor_t;
    // Velocity coupling
    auto u = cells.template get<names::Concentration0>().template getField<descriptors::VELOCITY>();
    T rho, pi[util::TensorVal<DESCRIPTOR>::n] { };
    cells.template get<names::NavierStokes>().computeAllMomenta(rho, u.data(), pi);
    cells.template get<names::Concentration0>().template setField<descriptors::VELOCITY>(u);
    // Stress tensor
    T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
    }
    T PiNeqNorm    = util::sqrt(PiNeqNormSqr);
    /// Molecular realaxation time
    T tau_mol_NS = T{1} / parameters.template get<OMEGA_NSE>();
    T tau_mol_AD = T{1} / parameters.template get<OMEGA_ADE>();
    /// Turbulent realaxation time
    T smagoPrefactor = parameters.template get<SMAGORINSKY_PREFACTOR>();
    T tau_turb_NS = T{0.5}*(util::sqrt(tau_mol_NS*tau_mol_NS + smagoPrefactor/rho*PiNeqNorm) - tau_mol_NS);
    // Schmidt number stabilization
    T tauTurbADPrefactor = descriptors::invCs2<T,DESCRIPTOR_ADE>() / descriptors::invCs2<T,DESCRIPTOR>() / parameters.template get<SCHMIDT>();
    T tau_turb_AD = tau_turb_NS * tauTurbADPrefactor;
    cells.template get<names::Concentration0>().template setField<descriptors::OMEGA>(T{1} / (tau_mol_AD + tau_turb_AD));
  }

};

/// VANS-ADE coupling
template<typename T>
struct VANSADECoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct PARTICLE_DIAMETER : public descriptors::FIELD_BASE<1> { };
  struct VISCOSITY : public descriptors::FIELD_BASE<1> { };
  struct CONV_VEL : public descriptors::FIELD_BASE<1> { };
  struct CONV_DENS : public descriptors::FIELD_BASE<1> { };
  struct PART_DENS : public descriptors::FIELD_BASE<1> { };
  struct DT : public descriptors::FIELD_BASE<1> { };
  struct CONV_FORCE : public descriptors::FIELD_BASE<1> { };
  struct CONV_MASS : public descriptors::FIELD_BASE<1> { };
  struct EARTH_ACC : public descriptors::FIELD_BASE<3> { };

  using parameters = meta::list<PARTICLE_DIAMETER,VISCOSITY,CONV_VEL,CONV_DENS,PART_DENS,DT,CONV_FORCE,CONV_MASS,EARTH_ACC>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
    T particleDiam = parameters.template get<PARTICLE_DIAMETER>();
    T visc = parameters.template get<VISCOSITY>();
    T convVel = parameters.template get<CONV_VEL>();
    T convDens = parameters.template get<CONV_DENS>();
    T convForce = parameters.template get<CONV_FORCE>();
    T convMass = parameters.template get<CONV_MASS>();
    T partDens = parameters.template get<PART_DENS>();
    T dt = parameters.template get<DT>();
    auto earthAcc = parameters.template get<EARTH_ACC>();
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Concentration0>();

    auto u_p = cellADE.template getField<descriptors::VELOCITY>();
    auto u_pXp = cellADE.neighbor({1,0,0}).template getField<descriptors::VELOCITY2>();
    auto u_pXm = cellADE.neighbor({-1,0,0}).template getField<descriptors::VELOCITY2>();
    auto u_pYp = cellADE.neighbor({0,1,0}).template getField<descriptors::VELOCITY2>();
    auto u_pYm = cellADE.neighbor({0,-1,0}).template getField<descriptors::VELOCITY2>();
    auto u_pZp = cellADE.neighbor({0,0,1}).template getField<descriptors::VELOCITY2>();
    auto u_pZm = cellADE.neighbor({0,0,-1}).template getField<descriptors::VELOCITY2>();
    auto u_f = cellADE.template getField<descriptors::VELOCITY>();
    T por = (1. - cellADE.computeRho());
    T rho { };
    cellNSE.computeRhoU(rho, u_f.data());

    T relU[3] = {0.};
    relU[0] = (u_f[0] - u_p[0]) * convVel;
    relU[1] = (u_f[1] - u_p[1]) * convVel;
    relU[2] = (u_f[2] - u_p[2]) * convVel;
    T norm = util::sqrt(relU[0]*relU[0] + relU[1]*relU[1] + relU[2]*relU[2]);
    T Re_p = norm * particleDiam / visc;
    T C_D = 0.;
    if(Re_p != 0.){
    if(Re_p <= 1000.){
      C_D = 24.*(1. + 0.15*util::pow(Re_p,0.687))/Re_p;
    }else{
      C_D = 0.44;
    }}
    T F_D[3] = {0.};
    //F_D[0] = C_D * 3./4. * rho * convDens * norm * relU[0] / particleDiam / partDens;
    //F_D[1] = C_D * 3./4. * rho * convDens * norm * relU[1] / particleDiam / partDens;
    //F_D[2] = C_D * 3./4. * rho * convDens * norm * relU[2] / particleDiam / partDens;
    //dragCoeff = (9.*converter_.getPhysViscosity()*converter_.getPhysDensity()*converter_.getConversionFactorTime()) / (2.*pRho_*pRadius_*pRadius_);
    F_D[0] = 9.*visc*convDens/partDens/particleDiam*relU[0];
    F_D[1] = 9.*visc*convDens/partDens/particleDiam*relU[1];
    F_D[2] = 9.*visc*convDens/partDens/particleDiam*relU[2];

    T mass = partDens * 3.14/4. * particleDiam * particleDiam;
    if(cellADE.computeRho() > 0.001){
    u_p[0] += dt * (F_D[0]/*mass*/ + earthAcc[0]*(partDens - convDens)/partDens) / convVel - u_p * 0.5 * (u_pXp - u_pXm);
    u_p[1] += dt * (F_D[1]/*mass*/ + earthAcc[1]*(partDens - convDens)/partDens) / convVel - u_p * 0.5 * (u_pYp - u_pYm);
    u_p[2] += dt * (F_D[2]/*mass*/ + earthAcc[2]*(partDens - convDens)/partDens) / convVel - u_p * 0.5 * (u_pZp - u_pZm);
    cellADE.template setField<descriptors::VELOCITY>(u_p);
    }

    T porPlusX = (1. - cellADE.neighbor({1,0,0}).computeRho());
    T porMinusX = (1. - cellADE.neighbor({-1,0,0}).computeRho());
    T porPlusY = (1. - cellADE.neighbor({0,1,0}).computeRho());
    T porMinusY = (1. - cellADE.neighbor({0,-1,0}).computeRho());
    T porPlusZ = (1. - cellADE.neighbor({0,0,1}).computeRho());
    T porMinusZ = (1. - cellADE.neighbor({0,0,-1}).computeRho());

    T porDX2 = porPlusX -2.*por + porMinusX;
    T porDY2 = porPlusY -2.*por + porMinusY;
    T porDZ2 = porPlusZ -2.*por + porMinusZ;

    int choice = 0;
    if( porDX2 != 0.) choice += 1;
    if( porDY2 != 0.) choice += 1;
    if( porDZ2 != 0.) choice += 1;

    T coeff = (choice == 1) *1./4. + (choice == 2) * 1./6. + (choice == 3) * 5./36.;
    T pressCorr = por + coeff * (porDX2 + porDY2 + porDZ2);
    if(pressCorr < 0.001){ pressCorr = 1.; }
    cellNSE.template setField<descriptors::PRESSCORR>(pressCorr);
    T press = rho / pressCorr / descriptors::invCs2<T,DESCRIPTOR>();
    T pressCorrForce[3] = {0.};
    pressCorrForce[0] = press * 0.5 * (porPlusX - porMinusX);
    pressCorrForce[1] = press * 0.5 * (porPlusY - porMinusY);
    pressCorrForce[2] = press * 0.5 * (porPlusZ - porMinusZ);

    T force[3] = {0.};
    if(cellADE.computeRho() > 0.001){
    force[0] = pressCorrForce[0] / rho - (1. - por) * por * F_D[0] / convForce * convMass * partDens / rho + por * earthAcc[0] / convVel * dt;
    force[1] = pressCorrForce[1] / rho - (1. - por) * por * F_D[1] / convForce * convMass * partDens / rho + por * earthAcc[1] / convVel * dt;
    force[2] = pressCorrForce[2] / rho - (1. - por) * por * F_D[2] / convForce * convMass * partDens / rho + por * earthAcc[2] / convVel * dt;
    }
    cellNSE.template setField<descriptors::FORCE>(force);
  }

};

/// granular flow
template<typename T>
struct GranularCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct FRICTION_ANGLE : public descriptors::FIELD_BASE<1> { };
  struct RHO0 : public descriptors::FIELD_BASE<1> { };
  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };

  using parameters = meta::list<FRICTION_ANGLE,RHO0,FORCE_PREFACTOR>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR>::n] { };
    cells.template get<names::NavierStokes>().computeAllMomenta(rho, u, pi);
    // Stress tensor
    T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
    }
    T PiNeqNorm = util::sqrt(PiNeqNormSqr);

    if(PiNeqNorm != T{0.}){
      auto omega = cells.template get<names::NavierStokes>().template getFieldPointer<descriptors::OMEGA>();
      PiNeqNorm /= (-rho/omega[0]/descriptors::invCs2<T,DESCRIPTOR>());
      T angle = parameters.template get<FRICTION_ANGLE>();
      T ps = ( T{0.5}*omega[0] - T{1} ) * T{1./3.} * ( pi[0] + pi[3] + pi[5] );
      T granularTau = ps/rho * util::sin(T{3.14}*angle/T{180}) / PiNeqNorm * descriptors::invCs2<T,DESCRIPTOR>() + T{0.5};
      omega[0] = T{1} / granularTau;
    }

    // computation of the bousinessq force
    auto force = cells.template get<names::NavierStokes>().template getFieldPointer<descriptors::FORCE>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();
    T densityDifference = (rho - parameters.template get<RHO0>()) / rho;
    for (unsigned iD = 0; iD < DESCRIPTOR::d; ++iD) {
      force[iD] = forcePrefactor[iD] * densityDifference;
    }
  }

};


}

#endif
