/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef DUAL_LB_HELPERS_H
#define DUAL_LB_HELPERS_H

#include "dualLatticeDescriptors.h"
#include "dynamics/lbm.h"
//#include "core/cell.h"
//#include "core/util.h"
//#include "utilities/vectorHelpers.h"


namespace olb {

/// All optimization code is contained in this namespace.
namespace opti {

// Forward declarations
template<typename T, typename DESCRIPTOR> struct dualLbDynamicsHelpers;
template<typename T, typename DESCRIPTOR> struct dualLbExternalHelpers;
template<typename T, typename DESCRIPTOR> struct dualLbLatticeHelpers;

/// This structure forwards the calls to the appropriate helper class
template<typename T, typename DESCRIPTOR>
struct dualLbHelpers {

  static T equilibrium(int iPop, int jPop, T rho, const T u[DESCRIPTOR::d])
  {
    return dualLbDynamicsHelpers<T,DESCRIPTOR>
           ::equilibrium(iPop, jPop, rho, u);
  }

  static T equilibrium2(int iPop, int jPop, T rho, const T u[DESCRIPTOR::d])
  {
    return dualLbDynamicsHelpers<T,DESCRIPTOR>
           ::equilibrium2(iPop, jPop, rho, u);
  }
  /*
    static T incEquilibrium(int iPop, const T j[DESCRIPTOR::d], const T jSqr, const T pressure) {
      return lbm<DESCRIPTOR>
             ::incEquilibrium(iPop, j, jSqr, pressure);
    }

    static void computeFneq ( ConstCell<T,DESCRIPTOR>& cell,
                              T fNeq[DESCRIPTOR::q], T rho, const T u[DESCRIPTOR::d] )
    {
      lbm<DESCRIPTOR>
      ::computeFneq(&cell[0], fNeq, rho, u);
    }

    static T bgkCollision(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T omega)
    {
      return lbm<DESCRIPTOR>
             ::bgkCollision(&cell[0], rho, u, omega);
    }

    static T incBgkCollision(Cell<T,DESCRIPTOR>& cell, T pressure, const T j[DESCRIPTOR::d], T omega)
    {
      return lbm<DESCRIPTOR>
             ::incBgkCollision(&cell[0], pressure, j, omega);
    }

    static T constRhoBgkCollision(Cell<T,DESCRIPTOR>& cell,
                                  T rho, const T u[DESCRIPTOR::d], T ratioRho, T omega)
    {
      return lbm<DESCRIPTOR>
             ::constRhoBgkCollision(&cell[0], rho, u, ratioRho, omega);
    }

    static T computeRho(ConstCell<T,DESCRIPTOR>& cell) {
      return lbm<DESCRIPTOR>
             ::computeRho(&cell[0]);
    }

    static void computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d] ) {
      lbm<DESCRIPTOR>
      ::computeJ(&cell[0], j);
    }
  */
  static void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d])
  {
    dualLbDynamicsHelpers<T,DESCRIPTOR>
    ::computeRhoU(&cell[0], rho, u);
  }

  static void computeRhoJ(ConstCell<T,DESCRIPTOR>& cell, T& rho, T j[DESCRIPTOR::d])
  {
    dualLbDynamicsHelpers<T,DESCRIPTOR>
    ::computeRhoJ(&cell[0], rho, j);
  }
  /*
    static void computeStress(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d],
                              T pi[util::TensorVal<DESCRIPTOR >::n] )
    {
      lbm<DESCRIPTOR>
      ::computeStress(&cell[0], rho, u, pi);
    }

    static void computeAllMomenta(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d],
                                  T pi[util::TensorVal<DESCRIPTOR >::n] )
    {
      lbm<DESCRIPTOR>
      ::computeAllMomenta(&cell[0], rho, u, pi);
    }

    static void modifyVelocity(Cell<T,DESCRIPTOR>& cell, const T newU[DESCRIPTOR::d]) {
      lbm<DESCRIPTOR>
      ::modifyVelocity(&cell[0], newU);
    }

    static void addExternalForce(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d], T omega, T amplitude=(T)1)
    {
      lbExternalHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, omega, amplitude);
    }
  */
};  // struct lbHelpers


/// All helper functions are inside this structure
template<typename T, typename DESCRIPTOR>
struct dualLbDynamicsHelpers {
  /// Computation of adjoint equilibrium distribution operator dEq_i/dF_j -> Mathias J. Krause
  static T equilibrium(int iPop, int jPop, T rho, const T u[DESCRIPTOR::d])
  {
    T eq = olb::equilibrium<DESCRIPTOR>::template secondOrder(iPop, rho, u) + descriptors::t<T,DESCRIPTOR>(iPop);
    T sum = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      sum += (u[iD] - descriptors::c<DESCRIPTOR>(jPop,iD))*(u[iD] - descriptors::c<DESCRIPTOR>(iPop,iD));
    } //std::cout<<u[0]<<std::endl;
    return (descriptors::invCs2<T,DESCRIPTOR>()*sum + T(1))/rho*eq;
  }

  /// Computation of adjoint equilibrium distribution operator dEq_i/dF_j (discrete approach) -> Geng Liu
  static T equilibrium2(int iPop, int jPop, T rho, const T u[DESCRIPTOR::d])
  {
    T u_c = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u_c += u[iD]*descriptors::c<DESCRIPTOR>(iPop,iD);
    }
    T sum = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      sum += (u[iD]*u[iD] + 2.*descriptors::c<DESCRIPTOR>(jPop,iD)*(descriptors::c<DESCRIPTOR>(iPop,iD)-u[iD]) )
            * descriptors::invCs2<T,DESCRIPTOR>()*0.5
         + (u_c*descriptors::c<DESCRIPTOR>(iPop,iD)*(2.*descriptors::c<DESCRIPTOR>(jPop,iD)-u[iD]))
            * descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*0.5;
    }
    return descriptors::t<T,DESCRIPTOR>(iPop)*(1.+sum);
  }

  /*
    static T incEquilibrium( int iPop, const T j[DESCRIPTOR::d],
                             const T jSqr, const T pressure )
    {
      T c_j = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        c_j += descriptors::c<DESCRIPTOR>(iPop,iD)*j[iD];
      }
      T rho = (T)1 + pressure*descriptors::invCs2<T,DESCRIPTOR>();
      return descriptors::t<T,DESCRIPTOR>(iPop) * ( rho +
                                     descriptors::invCs2<T,DESCRIPTOR>() * c_j +
                                     descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>()/(T)2 * c_j*c_j -
                                     descriptors::invCs2<T,DESCRIPTOR>()/(T)2 * jSqr
                                   ) - descriptors::t<T,DESCRIPTOR>(iPop);
    }

    static void computeFneq(T const* cell, T fNeq[DESCRIPTOR::q], T rho, const T u[DESCRIPTOR::d]) {
      const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fNeq[iPop] = cell[iPop] - equilibrium(iPop, rho, u, uSqr);
      }
    }

    /// BGK collision step
    static T bgkCollision(T* cell, T rho, const T u[DESCRIPTOR::d], T omega) {
      const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] *= (T)1-omega;
        cell[iPop] += omega * lbm<DESCRIPTOR>::equilibrium (
                        iPop, rho, u, uSqr );
      }
      return uSqr;
    }

    /// Incompressible BGK collision step
    static T incBgkCollision(T* cell, T pressure, const T j[DESCRIPTOR::d], T omega) {
      const T jSqr = util::normSqr<T,DESCRIPTOR::d>(j);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] *= (T)1-omega;
        cell[iPop] += omega * lbm<DESCRIPTOR>::incEquilibrium (
                        iPop, j, jSqr, pressure );
      }
      return jSqr;
    }

    /// BGK collision step with density correction
    static T constRhoBgkCollision(T* cell, T rho, const T u[DESCRIPTOR::d], T ratioRho, T omega) {
      const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        T feq = lbm<DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr );
        cell[iPop] =
          ratioRho*(feq+descriptors::t<T,DESCRIPTOR>(iPop))-descriptors::t<T,DESCRIPTOR>(iPop) +
          ((T)1-omega)*(cell[iPop]-feq);
      }
      return uSqr;
    }

    /// Computation of density
    static T computeRho(T const* cell) {
      T rho = T();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        rho += cell[iPop];
      }
      rho += (T)1;
      return rho;
    }

    /// Computation of momentum
    static void computeJ(T const* cell, T j[DESCRIPTOR::d]) {
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        j[iD] = T();
      }
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          j[iD] += cell[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
        }
      }
    }
  */
  /// Computation of hydrodynamic variables
  static void computeRhoU(T const* cell, T& rho, T u[DESCRIPTOR::d])
  {

    rho = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] = T();
    }
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      rho += cell[iPop];
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        u[iD] += cell[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
      }
    }
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] /= rho;
    }
  }

  /// Computation of hydrodynamic variables
  static void computeRhoJ(T const* cell, T& rho, T j[DESCRIPTOR::d])
  {

    lbm<DESCRIPTOR>::computeRhoJ(&cell[0], rho, j);
    rho-=T(1);
  }
  /*
    /// Computation of stress tensor
    static void computeStress(T const* cell, T rho, const T u[DESCRIPTOR::d],
                              T pi[util::TensorVal<DESCRIPTOR>::n] )
    {
      int iPi = 0;
      for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
        for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
          pi[iPi] = T();
          for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
            pi[iPi] += descriptors::c<DESCRIPTOR>(iPop,iAlpha) *
                       descriptors::c<DESCRIPTOR>(iPop,iBeta) * cell[iPop];
          }
          // stripe off equilibrium contribution
          pi[iPi] -= rho*u[iAlpha]*u[iBeta];
          if (iAlpha==iBeta) {
            pi[iPi] -= 1./descriptors::invCs2<T,DESCRIPTOR>()*(rho-(T)1);
          }
          ++iPi;
        }
      }
    }

    /// Computation of all hydrodynamic variables
    static void computeAllMomenta(T const* cell, T& rho, T u[DESCRIPTOR::d],
                                  T pi[util::TensorVal<DESCRIPTOR>::n] )
    {
      computeRhoU(cell, rho, u);
      computeStress(cell, rho, u, pi);
    }

    static void modifyVelocity(T* cell, const T newU[DESCRIPTOR::d]) {
      T rho, oldU[DESCRIPTOR::d];
      computeRhoU(cell, rho, oldU);
      const T oldUSqr = util::normSqr<T,DESCRIPTOR::d>(oldU);
      const T newUSqr = util::normSqr<T,DESCRIPTOR::d>(newU);
      for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
        cell[iPop] = cell[iPop]
                     - equilibrium(iPop, rho, oldU, oldUSqr)
                     + equilibrium(iPop, rho, newU, newUSqr);
      }
    }*/

};  // struct lbHelpers

/// Helper functions for dynamics that access external field
template<typename T, typename DESCRIPTOR>
struct dualLbExternalHelpers {
  /*
    /// Add a force term after BGK collision
    static void addExternalForce(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d], T omega, T amplitude) {
      static const int forceBeginsAt = DESCRIPTOR::template index<descriptors::FORCE>();
      T* force = cell.template getFieldPointer<descriptors::FORCE>();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        T c_u = T();
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
        }
        c_u *= descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>();
        T forceTerm = T();
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          forceTerm +=
            (   (descriptors::c<DESCRIPTOR>(iPop,iD)-u[iD]) * descriptors::invCs2<T,DESCRIPTOR>()
                + c_u * descriptors::c<DESCRIPTOR>(iPop,iD)
            )
            * force[iD];
        }
        forceTerm *= descriptors::t<T,DESCRIPTOR>(iPop);
        forceTerm *= 1-omega/(T)2;
        forceTerm *= amplitude;
        cell[iPop] += forceTerm;
      }
    }*/
};  // struct externalFieldHelpers

/// Helper functions with full-lattice access
template<typename T, typename DESCRIPTOR>
struct dualLbLatticeHelpers {
};

/// All boundary helper functions are inside this structure
template<typename T, typename DESCRIPTOR, int direction, int orientation>
struct DualBoundaryHelpers {
  /*
    static void computeStress (
      ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d],
      T pi[util::TensorVal<DESCRIPTOR >::n] )
    {
      typedef DESCRIPTOR L;
      const T uSqr = util::normSqr<T,L::d>(u);

      std::vector<int> const& onWallIndices = util::subIndex<L, direction, 0>();
      std::vector<int> const& normalIndices = util::subIndex<L, direction, orientation>();

      T fNeq[DESCRIPTOR::q];
      for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
        int iPop = onWallIndices[fIndex];
        fNeq[iPop] =
          cell[iPop] -
          lbm<DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
      }
      for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
        int iPop = normalIndices[fIndex];
        if (iPop == 0) {
          fNeq[iPop] = T();  // fNeq[0] will not be used anyway
        }
        else {
          fNeq[iPop] =
            cell[iPop] -
            lbm<DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
        }
      }

      int iPi = 0;
      for (int iAlpha=0; iAlpha<L::d; ++iAlpha) {
        for (int iBeta=iAlpha; iBeta<L::d; ++iBeta) {
          pi[iPi] = T();
          for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex)
          {
            const int iPop = onWallIndices[fIndex];
            pi[iPi] +=
              descriptors::c<L>(iPop)[iAlpha]*descriptors::c<L>(iPop)[iBeta]*fNeq[iPop];
          }
          for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex)
          {
            const int iPop = normalIndices[fIndex];
            pi[iPi] += (T)2 * descriptors::c<L>(iPop)[iAlpha]*descriptors::c<L>(iPop)[iBeta]*
                       fNeq[iPop];
          }
          ++iPi;
        }
      }
    }
  */
};  // struct boundaryHelpers


} // namespace opti

} // namespace olb

// The specialized code is directly included. That is because we never want
// it to be precompiled so that in both the precompiled and the
// "include-everything" version, the compiler can apply all the
// optimizations it wants.
//#include "lbHelpersD2Q9.h"
//#include "lbHelpersD3Q9.h"

#endif
