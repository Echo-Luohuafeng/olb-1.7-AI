/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Mathias J. Krause
 *                2024 Julius Jessberger, Shota Ito, Adrian Kummerlaender
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

/** \file
 * The description of optimization algorthims -- header file.
 */
/** \file
 * A collection of dynamics classes for dual LB methods
 * (e.g. dual BGK) with which a Cell object can be
 * instantiated -- header file.
 */

#ifndef DUAL_DYNAMICS_H
#define DUAL_DYNAMICS_H

#include "dualLbHelpers.h"

#include "utilities/vectorHelpers.h"

// All OpenLB code is contained in this namespace.
namespace olb {

namespace collision {

struct DualPorousBGK {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "DualPorousBGK";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;

    template <CONCEPT(MinimalCell) CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      // Revert for backwards-in-time propagation
      //cell.revert();
      for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
        V cell_iPop = cell[iPop];
        cell[iPop] = cell[descriptors::opposite<DESCRIPTOR>(iPop)];
        cell[descriptors::opposite<DESCRIPTOR>(iPop)] = cell_iPop;
      }

      // Preparation
      const auto pop_f =   cell.template getField<descriptors::F>();
      const auto dJdF  =   cell.template getFieldPointer<descriptors::DJDF>();
      const V d      =  cell.template getField<descriptors::POROSITY>();
      const V omega = parameters.template get<descriptors::OMEGA>();

      // Forward density and velocity
      V rho_f { };
      V u_f[DESCRIPTOR::d] { };
      //this->computeRhoU( pop_f.data(), rho_f, u_f);
      rho_f = V();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        u_f[iD] = V();
      }
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        rho_f += pop_f[iPop];
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          u_f[iD] += pop_f[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
        }
      }
      rho_f += (V)1;
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        u_f[iD] /= rho_f;
      }

      // Adjoint equilibrium
      V pheq[DESCRIPTOR::q] { };
      for (int i=0; i < DESCRIPTOR::q; ++i) {
        pheq[i] = V{0};
        for (int j=0; j < DESCRIPTOR::q; ++j) {
          V feq_j = equilibrium<DESCRIPTOR>::secondOrder(j, rho_f, u_f)
              + descriptors::t<V,DESCRIPTOR>(j);
          V dot_ij = V();
          for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
            dot_ij += ( descriptors::c<DESCRIPTOR>(j,iD) - d*u_f[iD] )
              * ( descriptors::c<DESCRIPTOR>(i,iD) - u_f[iD] );
          }
          pheq[i] += cell[j]*feq_j*( V{1} + descriptors::invCs2<V,DESCRIPTOR>()*d*dot_ij );
        }
        pheq[i] /= rho_f;
      }

      // Collision
      for (int i=0; i < DESCRIPTOR::q; ++i) {
        cell[i] = cell[i] - omega*( cell[i] - pheq[i] ) + dJdF[i];
      }

      // Statistics
      V rho_phi, u_phi[DESCRIPTOR::d];
      MomentaF().computeRhoU( cell, rho_phi, u_phi );
      V uSqr_phi = util::normSqr<V,DESCRIPTOR::d>( u_phi );

      // Undo revert
      //cell.revert();
      for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
        V cell_iPop = cell[iPop];
        cell[iPop] = cell[descriptors::opposite<DESCRIPTOR>(iPop)];
        cell[descriptors::opposite<DESCRIPTOR>(iPop)] = cell_iPop;
      }
      return {rho_phi, uSqr_phi};
    };
  };
};

}

template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using DualPorousBGKDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::DualPorousBGK
>;

/// All optimization code is contained in this namespace.
namespace opti {

template<typename T> class Controller;
template<typename T, typename DESCRIPTOR> class DualController;

/// Implementation of the dual BGK collision step with external force
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
class DualForcedBGKdynamics : public legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA> {
public:
  /// Constructor
  DualForcedBGKdynamics(T omega_);
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const;
  /// Collision step
  virtual CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell);
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
  /// Does nothing
  virtual void defineRho(Cell<T,DESCRIPTOR>& cell, T rho);
private:
  T omega;  ///< relaxation parameter
};

////////////////////// Class ForcedBGKdynamics /////////////////////////

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA>
DualForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::DualForcedBGKdynamics (
  T omega_)
  : legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA>(),
    omega(omega_)
{
  this->getName() = "DualForcedBGKdynamics";
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void DualForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  // The meaning of this function is not clear. Please contact Julius.
  throw std::bad_function_call();
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T DualForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const
{
  return equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> DualForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::collide (
  Cell<T,DESCRIPTOR>& cell)
{
  cell.revert();

  // Preparation
  auto dJdF = cell.template getFieldPointer<descriptors::DJDF>();
  auto f = cell.template getField<descriptors::F>();
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  T rho_f;
  T u_f[DESCRIPTOR::d];

  dualLbDynamicsHelpers<T,DESCRIPTOR>::computeRhoU(f.data(), rho_f, u_f);
  T eq_phi[DESCRIPTOR::q];
  T force_phi[DESCRIPTOR::q];

  // Force
  T F1_phi[DESCRIPTOR::q][DESCRIPTOR::q];
  T F2_phi[DESCRIPTOR::q][DESCRIPTOR::q];
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    T f_c = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      f_c += force[iD]*descriptors::c<DESCRIPTOR>(iPop,iD);
    }
    for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
      T sum = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        sum += (f_c*descriptors::c<DESCRIPTOR>(iPop,iD)-force[iD]/(T)descriptors::invCs2<T,DESCRIPTOR>())*descriptors::c<DESCRIPTOR>(jPop,iD);
      }
      F1_phi[iPop][jPop] = descriptors::t<T,DESCRIPTOR>(iPop)*descriptors::invCs2<T,DESCRIPTOR>()*f_c;
      F2_phi[iPop][jPop] = descriptors::t<T,DESCRIPTOR>(iPop)*descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*sum;
    }
  }

  // Collision preperation
  for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
    eq_phi[jPop] = T();
    force_phi[jPop] = T();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      eq_phi[jPop] += cell[iPop]*dualLbHelpers<T,DESCRIPTOR>::equilibrium(iPop, jPop, rho_f, u_f);
      force_phi[jPop] += cell[iPop]*(F1_phi[iPop][jPop] + (1.-.5*omega)*F2_phi[iPop][jPop]);
    }
  }

  // Collision
  for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
    cell[jPop] = cell[jPop] - omega*(cell[jPop]-eq_phi[jPop]) + force_phi[jPop] + dJdF[jPop];
  }

  // Incrementing statistic values for convergence
  T phi2 = T();
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    phi2 += cell[iPop]*cell[iPop];
  }
  cell.revert();

  return {T(1) + cell[0], phi2};
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T DualForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void DualForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::setOmega(T omega_)
{
  omega = omega_;
}


///////////// Porous //////////////////////////

/// Implementation of the dual BGK collision step with external force
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
class DualPorousBGKdynamics : public legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA> {
public:
  using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;

  /// Constructor
  DualPorousBGKdynamics(T omega_);
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const;
  /// Compute moments
  virtual void computeRhoU(const T fPop[DESCRIPTOR::q], T& rho, T u[DESCRIPTOR::d]) const;
  /// Collision step
  virtual CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell);
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
  /// Does nothing
  virtual void defineRho(Cell<T,DESCRIPTOR>& cell, T rho);
private:
  T omega;  ///< relaxation parameter
};

////////////////////// Class DualPorousBGKdynamics /////////////////////////

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA>
DualPorousBGKdynamics<T,DESCRIPTOR,MOMENTA>::DualPorousBGKdynamics (
  T omega_)
  : legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA>(),
    omega(omega_)
{
  this->getName() = "DualPorousBGKdynamics";
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::POROSITY>() );
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void DualPorousBGKdynamics<T,DESCRIPTOR,MOMENTA>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  // The meaning of this function is not clear. Please contact Julius.
  throw std::bad_function_call();
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T DualPorousBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeEquilibrium(
  int iPop, T rho, const T u[DESCRIPTOR::d]) const
{
  // Compute equilibrium minus weights
  // used in iniEquilibrium
  return equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u);
}

// ComputeRhoU for array instead of cell
template<typename T, typename DESCRIPTOR, typename MOMENTA>
void DualPorousBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeRhoU(
  const T fPop[DESCRIPTOR::q], T& rho_f, T u_f[DESCRIPTOR::d]) const
{
  rho_f = T();
  for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
    u_f[iD] = T();
  }
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    rho_f += fPop[iPop];
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u_f[iD] += fPop[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
    }
  }
  rho_f += (T)1;
  for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
    u_f[iD] /= rho_f;
  }
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> DualPorousBGKdynamics<T,DESCRIPTOR,MOMENTA>::collide (
  Cell<T,DESCRIPTOR>& cell)
{
  // Revert for backwards-in-time propagation
  cell.revert();

  // Preparation
  auto pop_f =   cell.template getField<descriptors::F>();
  auto dJdF  =   cell.template getFieldPointer<descriptors::DJDF>();
  T d      =  cell.template getField<descriptors::POROSITY>();

  // Forward density and velocity
  T rho_f, u_f[DESCRIPTOR::d];
  this->computeRhoU( pop_f.data(), rho_f, u_f);

  // Adjoint equilibrium
  T pheq[DESCRIPTOR::q];
  for (int i=0; i < DESCRIPTOR::q; ++i) {
    pheq[i] = T{0};
    for (int j=0; j < DESCRIPTOR::q; ++j) {
      T feq_j = computeEquilibrium(j, rho_f, u_f)
                + descriptors::t<T,DESCRIPTOR>(j);
      T dot_ij = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        dot_ij += ( descriptors::c<DESCRIPTOR>(j,iD) - d*u_f[iD] )
                  * ( descriptors::c<DESCRIPTOR>(i,iD) - u_f[iD] );
      }
      pheq[i] += cell[j]*feq_j*( T{1} + descriptors::invCs2<T,DESCRIPTOR>()*d*dot_ij );
    }
    pheq[i] /= rho_f;
  }

  // Collision
  for (int i=0; i < DESCRIPTOR::q; ++i) {
    cell[i] = cell[i] - omega*( cell[i] - pheq[i] ) + dJdF[i];
  }

  // Statistics
  T rho_phi, u_phi[DESCRIPTOR::d];
  MomentaF().computeRhoU( cell, rho_phi, u_phi );
  T uSqr_phi = util::normSqr<T,DESCRIPTOR::d>( u_phi );

  // Undo revert
  cell.revert();
  return {rho_phi, uSqr_phi};
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T DualPorousBGKdynamics<T,DESCRIPTOR,MOMENTA>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void DualPorousBGKdynamics<T,DESCRIPTOR,MOMENTA>::setOmega(T omega_)
{
  omega = omega_;
}

} // namespace opti

} // namespace olb

#endif

#include "dualDynamics.cse.h"
