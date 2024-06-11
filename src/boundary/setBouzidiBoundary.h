/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Michael Crocoll, Adrian Kummerlaender, Shota Ito
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

#ifndef SET_BOUZIDI_BOUNDARY_H
#define SET_BOUZIDI_BOUNDARY_H

// M. Bouzidi, M. Firdaouss, and P. Lallemand.
// Momentum transfer of a Boltzmann-lattice fluid with boundaries.
// Physics of Fluids 13, 3452 (2001).
// DOI: 10.1063/1.1399290

// D. Yu, R. Mei, and W. Shyy.
// A Unified Boundary Treatment in Lattice Boltzmann Method.
// Session: FD-27: CFD Methodology III (2012)
// DOI: 10.2514/6.2003-953

namespace olb {

namespace descriptors {

/// Interpolated Bounce Back (Bouzidi) distance field
struct BOUZIDI_DISTANCE : public descriptors::FIELD_BASE<0,0,1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<BOUZIDI_DISTANCE>()>(-1);
  }
};

/// Interpolated Bounce Back (Bouzidi) velocity coefficient field
struct BOUZIDI_VELOCITY : public descriptors::FIELD_BASE<0,0,1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<BOUZIDI_VELOCITY>()>(0);
  }
};

/// Interpolated Bounce Back (Bouzidi) for ADE Dirichlet field
struct BOUZIDI_ADE_DIRICHLET : public descriptors::FIELD_BASE<0,0,1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<BOUZIDI_ADE_DIRICHLET>()>(0);
  }
};
}

/// Post processor for the zero-velocity Bouzidi boundary
class BouzidiPostProcessor {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return -1;
  }

  template <typename CELL, typename V = typename CELL::value_t>
  void apply(CELL& x_b) any_platform {
    using DESCRIPTOR = typename CELL::descriptor_t;
    const auto q = x_b.template getFieldPointer<descriptors::BOUZIDI_DISTANCE>();
    for (int iPop = 1; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
      // update missing population if valid bouzidi distance
      if (q[iPop] > V{0}) {
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
        auto x_s = x_b.neighbor(c);                                             // solid side neighbor
        auto x_f = x_b.neighbor(descriptors::c<DESCRIPTOR>(iPop_opposite));     // fluid side neighbor opposite to the missing population

        x_b[iPop_opposite] = (q[iPop] <= V{0.5}) // cut is closer to the fluid cell
                           * (V{2} * q[iPop] * x_s[iPop] + (V{1} - V{2} * q[iPop]) * x_b[iPop])
                           + (q[iPop] >  V{0.5}) // cut is closer to the solid cell
                           * (V{0.5} / q[iPop] * x_s[iPop] + V{0.5} * (V{2} * q[iPop] - V{1}) / q[iPop] * x_f[iPop_opposite]);
      }
      // if intersection point is on the cell then fall back to full-way bounce back
      else if (q[iPop] == V{0}) {
        const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
        x_b[iPop_opposite] = x_b[iPop];
      }
    }
  }
};

class BouzidiAdeDirichletPostProcessor {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return -1;
  }

  template <typename CELL, typename V = typename CELL::value_t>
  void apply(CELL& x_b) any_platform {
    using DESCRIPTOR = typename CELL::descriptor_t;
    const auto q = x_b.template getFieldPointer<descriptors::BOUZIDI_DISTANCE>();
    const auto phi_d = x_b.template getFieldPointer<descriptors::BOUZIDI_ADE_DIRICHLET>();
    V f = V{1/3};       // D2Q5
    if (descriptors::q<DESCRIPTOR>() == 7) {
      f = V{0.25};      // D3Q7
    }
    for (int iPop = 1; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
      // update missing population if valid bouzidi distance
      if (q[iPop] > V{0}) {
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
        auto x_s = x_b.neighbor(c);                                             // solid side neighbor
        auto x_f = x_b.neighbor(descriptors::c<DESCRIPTOR>(iPop_opposite));     // fluid side neighbor opposite to the missing population
        auto source_d = phi_d[iPop] * f;                                        // source term set by the dirichlet condition
        auto t_i = descriptors::t<V,DESCRIPTOR>(iPop);
        auto t_iopp = descriptors::t<V,DESCRIPTOR>(iPop_opposite);

        x_b[iPop_opposite] = (q[iPop] <= V{0.5}) // cut is closer to the fluid cell
                           * (V{-2} * q[iPop] * (x_s[iPop] + t_i) + (V{2} * q[iPop] - V{1}) * (x_b[iPop] + t_i) + source_d)
                           + (q[iPop] >  V{0.5}) // cut is closer to the solid cell
                           * (V{-1} / (V{2} * q[iPop]) * (x_s[iPop] + t_i) + (V{1} - V{1} / (V{2} * q[iPop])) * (x_f[iPop_opposite] + t_iopp) + (V{1} / (V{2} * q[iPop])) * source_d)
                           - t_iopp;
      }
      // if intersection point is on the cell then fall back to full-way bounce back
      else if (q[iPop] == V{0}) {
        const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
        auto source_d = phi_d[iPop] * f;
        auto t_i = descriptors::t<V,DESCRIPTOR>(iPop);
        auto t_iopp = descriptors::t<V,DESCRIPTOR>(iPop_opposite);
        x_b[iPop_opposite] = -(x_b[iPop] + t_i) + source_d - t_iopp;
      }
    }
  }
};

/// Post processor for the velocity Bouzidi boundary
class BouzidiVelocityPostProcessor {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return -1;
  }

  template <typename CELL, typename V = typename CELL::value_t>
  void apply(CELL& x_b) any_platform {
    using DESCRIPTOR = typename CELL::descriptor_t;
    const auto q = x_b.template getFieldPointer<descriptors::BOUZIDI_DISTANCE>();
    const auto veloCoeff = x_b.template getFieldPointer<descriptors::BOUZIDI_VELOCITY>();
    for (int iPop=1; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
      // update missing population if valid bouzidi distance
      if (q[iPop] > 0) {
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
        auto x_s = x_b.neighbor(c);                                             // solid side neighbor
        auto x_f = x_b.neighbor(descriptors::c<DESCRIPTOR>(iPop_opposite));     // fluid side neighbor opposite to the missing population
        auto veloTerm = veloCoeff[iPop] * (descriptors::t<V,DESCRIPTOR>(iPop)) * (descriptors::invCs2<V,DESCRIPTOR>());

        x_b[iPop_opposite] = (q[iPop] <= V{0.5}) // cut is closer to the fluid cell
                           * (V{2} * q[iPop] * x_s[iPop] + (V{1} - V{2} * q[iPop]) * x_b[iPop] - V{2} * veloTerm)
                           + (q[iPop] >  V{0.5}) // cut is closer to the solid cell
                           * (V{0.5} / q[iPop] * x_s[iPop] + V{0.5} * (V{2} * q[iPop] - V{1}) / q[iPop] * x_f[iPop_opposite] - V{1}/q[iPop] * veloTerm);
      }
      // if intersection point is on the cell then fall back to full-way bounce back
      else if (q[iPop] == V{0}) {
        const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
        auto veloTerm = veloCoeff[iPop] * (descriptors::t<V,DESCRIPTOR>(iPop)) * (descriptors::invCs2<V,DESCRIPTOR>());
        x_b[iPop_opposite] = x_b[iPop] - V{2} * veloTerm;
      }
    }
  }
};

// Post processor for the zero-velocity Yu IBB scheme
class YuPostProcessor {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return -1;
  }

  template <typename CELL, typename V = typename CELL::value_t>
  void apply(CELL& x_b) any_platform {
    using DESCRIPTOR = typename CELL::descriptor_t;
    const auto q = x_b.template getFieldPointer<descriptors::BOUZIDI_DISTANCE>();
    for (int iPop=1; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
      if (q[iPop] >= 0) {
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
        auto x_s = x_b.neighbor(c);                                         // solid cell inside obstacle material
        auto x_f = x_b.neighbor(descriptors::c<DESCRIPTOR>(iPop_opposite)); // fluid boundary cell
        auto f_tmp = x_b[iPop] + q[iPop]*(x_s[iPop] - x_b[iPop]);           // population at fictitious ghost particle
        x_b[iPop_opposite] = f_tmp + q[iPop]/(V{1}+q[iPop]) * (x_f[iPop_opposite] - f_tmp);
      }
    }
  }
};

/// Set Bouzidi boundary on indicated cells of sLattice
template<typename T, typename DESCRIPTOR, typename OPERATOR = BouzidiPostProcessor>
void setBouzidiBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,
                        FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& boundaryIndicator,
                        FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& bulkIndicator,
                        IndicatorF<T,DESCRIPTOR::d>&                   indicatorAnalyticalBoundary)
{
  int _overlap = 1;
  OstreamManager clout(std::cout, "BouzidiBoundarySetter");
  auto& load = sLattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    setBouzidiBoundary<T,DESCRIPTOR,OPERATOR>(sLattice.getBlock(iC),
                          (bulkIndicator->getBlockIndicatorF(iC)).getBlockGeometry(),
                          boundaryIndicator->getBlockIndicatorF(iC),
                          bulkIndicator->getBlockIndicatorF(iC),
                          indicatorAnalyticalBoundary);
  }
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(boundaryIndicator)>(boundaryIndicator), _overlap);
  clout.setMultiOutput(false);
}

/// Set Bouzidi boundary on material cells of sLattice
template<typename T, typename DESCRIPTOR, typename OPERATOR = BouzidiPostProcessor>
void setBouzidiBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,
                        SuperGeometry<T,DESCRIPTOR::d>& superGeometry,
                        int materialOfSolidObstacle,
                        IndicatorF<T,DESCRIPTOR::d>& indicatorAnalyticalBoundary,
                        std::vector<int> bulkMaterials = std::vector<int>(1,1))
{
  //Getting the indicators by material numbers and calling the superLattice method via the indicators:
  setBouzidiBoundary<T,DESCRIPTOR,OPERATOR>(sLattice,
                      FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>(superGeometry.getMaterialIndicator(materialOfSolidObstacle)),
                      FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>(superGeometry.getMaterialIndicator(std::move(bulkMaterials))),
                      indicatorAnalyticalBoundary);
}


/// Set Bouzidi boundary on indicated cells of block lattice
template<typename T, typename DESCRIPTOR, typename OPERATOR = BouzidiPostProcessor>
void setBouzidiBoundary(BlockLattice<T,DESCRIPTOR>& block,
                        BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
                        BlockIndicatorF<T,DESCRIPTOR::d>& boundaryIndicator,
                        BlockIndicatorF<T,DESCRIPTOR::d>& bulkIndicator,
                        IndicatorF<T,DESCRIPTOR::d>& indicatorAnalyticalBoundary,
                        bool verbose = false)
{
  OstreamManager clout(std::cout, "BouzidiBoundarySetter");
  clout.setMultiOutput(true);

  const T deltaR = blockGeometry.getDeltaR();
  // for each solid cell: all of its fluid neighbors need population updates
  block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> solidLatticeR) {
    // Check if cell is solid cell
    if (boundaryIndicator(solidLatticeR)) {
      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        Vector<T,DESCRIPTOR::d> boundaryLatticeR(solidLatticeR + descriptors::c<DESCRIPTOR>(iPop));
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        const auto iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);

        if (blockGeometry.isInside(boundaryLatticeR)) {
          T boundaryPhysR[DESCRIPTOR::d] { };
          blockGeometry.getPhysR(boundaryPhysR,boundaryLatticeR);
          // check if neighbor is fluid cell
          if (bulkIndicator(boundaryLatticeR)) {
            T dist = -1;    // distance to boundary
            T qIpop = -1;   // normed distance (Bouzidi distance) to boundary
            const T norm = deltaR * util::norm<DESCRIPTOR::d>(descriptors::c<DESCRIPTOR>(iPop));
            auto direction = -deltaR * c; // vector pointing from the boundary cell to the solid cell

            // Check if distance calculation was performed correctly
            if (indicatorAnalyticalBoundary.distance(dist, boundaryPhysR, direction, blockGeometry.getIcGlob())) {
              qIpop = dist / norm;

              // if distance function returned a dist. not suitable for Bouzidi -> fall-back
              if ((qIpop < 0) || (qIpop > 1)) {
                if(verbose) {
                  clout << "Error, non suitable dist. at lattice: (" << boundaryLatticeR << "), physical: (" << blockGeometry.getPhysR(boundaryLatticeR) << "), direction " << iPop << ". Fall-back to bounce-back." << std::endl;
                }

                // fall-back: half-way bounce back
                qIpop = 0.5;
              }
            }
            // if distance function couldn't compute any distance -> fall-back
            else {
              if(verbose) {
                clout << "Error, no boundary found at lattice:(" << boundaryLatticeR << "), physical: (" << blockGeometry.getPhysR(boundaryLatticeR) << "), direction: " << iPop << ".Fall-back to bounce-back." << std::endl;
              }

              // fall-back: half-way bounce back
              qIpop = 0.5;
            }

            // double check
            if (qIpop >= 0) {
              // Bouzidi require the fluid side neighbor of the boundary cell also to be fluid
              if (bulkIndicator(boundaryLatticeR + descriptors::c<DESCRIPTOR>(iPop))) {
                // Standard case, c.f. Bouzidi paper, setting Bouzidi-distance
                block.get(boundaryLatticeR).template setFieldComponent<descriptors::BOUZIDI_DISTANCE>(iPop_opposite, qIpop);
              }
              else {
                // If no fluid cell found: fall-back to bounce-back
                block.get(boundaryLatticeR).template setFieldComponent<descriptors::BOUZIDI_DISTANCE>(iPop_opposite, T{0.5});
              }
              // Initialize velocity coefficients if necessary
              if constexpr (std::is_same_v<OPERATOR, BouzidiVelocityPostProcessor>) {
                block.get(boundaryLatticeR).template setFieldComponent<descriptors::BOUZIDI_VELOCITY>(iPop_opposite, 0);
              }
              // Initialize ade dirichlet coefficients if necessary
              if constexpr (std::is_same_v<OPERATOR, BouzidiAdeDirichletPostProcessor>) {
                block.get(boundaryLatticeR).template setFieldComponent<descriptors::BOUZIDI_ADE_DIRICHLET>(iPop_opposite, 0);
              }
              // Setting up the post processor, if this cell does not have one yet.
              if (!block.isPadding(boundaryLatticeR)) {
                block.addPostProcessor(typeid(stage::PostStream),
                                      boundaryLatticeR,
                                      meta::id<OPERATOR>{});
              }
            }
          }
          // if neigbour cell is not fluid
          else {
            // check if neighbor cell is not solid
            if (blockGeometry.getMaterial(boundaryLatticeR) != 0) {
              // fall-back to half-way bounce-back
              block.get(boundaryLatticeR).template setFieldComponent<descriptors::BOUZIDI_DISTANCE>(iPop_opposite, T{0.5});
              // Initialize velocity coefficients if necessary
              if constexpr (std::is_same_v<OPERATOR, BouzidiVelocityPostProcessor>) {
                block.get(boundaryLatticeR).template setFieldComponent<descriptors::BOUZIDI_VELOCITY>(iPop_opposite, 0);
              }
              // Initialize velocity coefficients if necessary
              if constexpr (std::is_same_v<OPERATOR, BouzidiAdeDirichletPostProcessor>) {
                block.get(boundaryLatticeR).template setFieldComponent<descriptors::BOUZIDI_ADE_DIRICHLET>(iPop_opposite, 0);
              }
              if (!block.isPadding(boundaryLatticeR)) {
                block.addPostProcessor(typeid(stage::PostStream),
                                      boundaryLatticeR,
                                      meta::id<OPERATOR>{});
              }
            }
          }
        }
      }
    }
  });
}

/// Set Bouzidi velocity boundary on material cells of sLattice
template<typename T, typename DESCRIPTOR>
void setBouzidiVelocity(SuperLattice<T,DESCRIPTOR>& sLattice,
                        SuperGeometry<T,DESCRIPTOR::d>& superGeometry, int material,
                        AnalyticalF<DESCRIPTOR::d,T,T>& u,
                        std::vector<int> bulkMaterials = std::vector<int>(1,1))
{
  setBouzidiVelocity<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material), u, bulkMaterials);
}

/// Set Bouzidi velocity boundary on indicated cells of sLattice
template<typename T, typename DESCRIPTOR>
void setBouzidiVelocity(SuperLattice<T,DESCRIPTOR>& sLattice,
                        FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& boundaryIndicator,
                        AnalyticalF<DESCRIPTOR::d,T,T>& u,
                        std::vector<int> bulkMaterials = std::vector<int>(1,1))
{
  setBouzidiVelocity<T,DESCRIPTOR>(sLattice, std::forward<decltype(boundaryIndicator)>(boundaryIndicator),
                               boundaryIndicator->getSuperGeometry().getMaterialIndicator(std::move(bulkMaterials)),
                               u);
}

/// Set Bouzidi velocity boundary on indicated cells of sLattice
template<typename T, typename DESCRIPTOR>
void setBouzidiVelocity(SuperLattice<T,DESCRIPTOR>& sLattice,
                        FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& boundaryIndicator,
                        FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& bulkIndicator,
                        AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  auto& load = sLattice.getLoadBalancer();
  auto& cuboidGeometry = boundaryIndicator->getSuperGeometry().getCuboidGeometry();
  for (int iCloc = 0; iCloc < load.size(); ++iCloc) {
    auto& cuboid = cuboidGeometry.get(load.glob(iCloc));
    setBouzidiVelocity<T,DESCRIPTOR>(sLattice.getBlock(iCloc),
                                     boundaryIndicator->getBlockIndicatorF(iCloc),
                                     bulkIndicator->getBlockIndicatorF(iCloc),
                                     u,
                                     cuboid);
  }
}

/// Set Bouzidi velocity boundary on indicated cells of block lattice
template<typename T, typename DESCRIPTOR>
void setBouzidiVelocity(BlockLattice<T,DESCRIPTOR>& block,
                        BlockIndicatorF<T,DESCRIPTOR::d>& boundaryIndicator,
                        BlockIndicatorF<T,DESCRIPTOR::d>& bulkIndicator,
                        AnalyticalF<DESCRIPTOR::d,T,T>& u,
                        Cuboid<T,DESCRIPTOR::d>& cuboid)
{
  const T deltaR = cuboid.getDeltaR();
  block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> solidLatticeR) {
    if (boundaryIndicator(solidLatticeR)) {
      for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
        Vector<T,DESCRIPTOR::d> boundaryLatticeR(solidLatticeR + descriptors::c<DESCRIPTOR>(iPop));
        if (block.isInside(boundaryLatticeR)) {
          const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
          auto x_b = block.get(boundaryLatticeR);
          const auto opp_bouzidi_dist = x_b.template getFieldComponent<descriptors::BOUZIDI_DISTANCE>(iPop_opposite);

          // check if distance from the fluid cell to the solid cell is a valid bouzidi distance
          if (opp_bouzidi_dist >= 0) {
            T wallVelocity[DESCRIPTOR::d] = { };
            T physicalIntersection[DESCRIPTOR::d] = { };
            T boundaryPhysR[DESCRIPTOR::d] { };
            cuboid.getPhysR(boundaryPhysR, boundaryLatticeR);

            // calculating the intersection of the boundary with the missing link in physical coordinates
            for ( int i = 0; i< DESCRIPTOR::d; ++i) {
              physicalIntersection[i] = boundaryPhysR[i] + opp_bouzidi_dist * deltaR * descriptors::c<DESCRIPTOR>(iPop_opposite,i);
            }

            //Calculating the velocity at the wall intersection
            u(wallVelocity, physicalIntersection);
            const auto c = descriptors::c<DESCRIPTOR>(iPop_opposite);
            T vel_coeff = c * Vector<T,DESCRIPTOR::d>(wallVelocity);

            // set computed velocity into the bouzidi velocity field
            block.get(boundaryLatticeR).template setFieldComponent<descriptors::BOUZIDI_VELOCITY>(iPop_opposite, vel_coeff);
          }
        }
      }
    }
  });
}

template<typename T, typename DESCRIPTOR>
void setBouzidiAdeDirichlet(SuperLattice<T,DESCRIPTOR>& sLattice,
                            SuperGeometry<T,DESCRIPTOR::d>& superGeometry, int material,
                            T phi_d,
                            std::vector<int> bulkMaterials = std::vector<int>(1,1))
{
  setBouzidiAdeDirichlet<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material), phi_d, bulkMaterials);
}

template<typename T, typename DESCRIPTOR>
void setBouzidiAdeDirichlet(SuperLattice<T,DESCRIPTOR>& sLattice,
                            SuperGeometry<T,DESCRIPTOR::d>& superGeometry, int material,
                            AnalyticalF<DESCRIPTOR::d,T,T>& phi_d,
                            std::vector<int> bulkMaterials = std::vector<int>(1,1))
{
  setBouzidiAdeDirichlet<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material), phi_d, bulkMaterials);
}

template<typename T, typename DESCRIPTOR>
void setBouzidiAdeDirichlet(SuperLattice<T,DESCRIPTOR>& sLattice,
                            FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& boundaryIndicator,
                            FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& bulkIndicator,
                            T phi_d)
{
  auto& load = sLattice.getLoadBalancer();
  auto& cuboidGeometry = boundaryIndicator->getSuperGeometry().getCuboidGeometry();
  for (int iCloc = 0; iCloc < load.size(); ++iCloc) {
    auto& cuboid = cuboidGeometry.get(load.glob(iCloc));
    setBouzidiAdeDirichlet<T,DESCRIPTOR>(sLattice.getBlock(iCloc),
                                     boundaryIndicator->getBlockIndicatorF(iCloc),
                                     bulkIndicator->getBlockIndicatorF(iCloc),
                                     phi_d,
                                     cuboid);
  }
}

template<typename T, typename DESCRIPTOR>
void setBouzidiAdeDirichlet(SuperLattice<T,DESCRIPTOR>& sLattice,
                            FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& boundaryIndicator,
                            FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& bulkIndicator,
                            AnalyticalF<DESCRIPTOR::d,T,T>& phi_d)
{
  auto& load = sLattice.getLoadBalancer();
  auto& cuboidGeometry = boundaryIndicator->getSuperGeometry().getCuboidGeometry();
  for (int iCloc = 0; iCloc < load.size(); ++iCloc) {
    auto& cuboid = cuboidGeometry.get(load.glob(iCloc));
    setBouzidiAdeDirichlet<T,DESCRIPTOR>(sLattice.getBlock(iCloc),
                                     boundaryIndicator->getBlockIndicatorF(iCloc),
                                     bulkIndicator->getBlockIndicatorF(iCloc),
                                     phi_d,
                                     cuboid);
  }
}

template<typename T, typename DESCRIPTOR>
void setBouzidiAdeDirichlet(SuperLattice<T,DESCRIPTOR>& sLattice,
                            FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& boundaryIndicator,
                            T phi_d,
                            std::vector<int> bulkMaterials = std::vector<int>(1,1))
{
  setBouzidiAdeDirichlet<T,DESCRIPTOR>(sLattice, std::forward<decltype(boundaryIndicator)>(boundaryIndicator),
                               boundaryIndicator->getSuperGeometry().getMaterialIndicator(std::move(bulkMaterials)),
                               phi_d);
}

template<typename T, typename DESCRIPTOR>
void setBouzidiAdeDirichlet(SuperLattice<T,DESCRIPTOR>& sLattice,
                            FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& boundaryIndicator,
                            AnalyticalF<DESCRIPTOR::d,T,T>& phi_d,
                            std::vector<int> bulkMaterials = std::vector<int>(1,1))
{
  setBouzidiAdeDirichlet<T,DESCRIPTOR>(sLattice, std::forward<decltype(boundaryIndicator)>(boundaryIndicator),
                               boundaryIndicator->getSuperGeometry().getMaterialIndicator(std::move(bulkMaterials)),
                               phi_d);
}

template<typename T, typename DESCRIPTOR>
void setBouzidiAdeDirichlet(BlockLattice<T,DESCRIPTOR>& block,
                            BlockIndicatorF<T,DESCRIPTOR::d>& boundaryIndicator,
                            BlockIndicatorF<T,DESCRIPTOR::d>& bulkIndicator,
                            T phi_d,
                            Cuboid<T,DESCRIPTOR::d>& cuboid)
{
  block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> solidLatticeR) {
    if (boundaryIndicator(solidLatticeR)) {
      for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
        Vector<T,DESCRIPTOR::d> boundaryLatticeR(solidLatticeR + descriptors::c<DESCRIPTOR>(iPop));
        if (block.isInside(boundaryLatticeR)) {
          const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
          auto x_b = block.get(boundaryLatticeR);
          const auto opp_bouzidi_dist = x_b.template getFieldComponent<descriptors::BOUZIDI_DISTANCE>(iPop_opposite);

          // check if distance from the fluid cell to the solid cell is a valid bouzidi distance
          if (opp_bouzidi_dist >= 0) {
            // set computed velocity into the bouzidi velocity field
            block.get(boundaryLatticeR).template setFieldComponent<descriptors::BOUZIDI_ADE_DIRICHLET>(iPop_opposite, phi_d);
          }
        }
      }
    }
  });
}

template<typename T, typename DESCRIPTOR>
void setBouzidiAdeDirichlet(BlockLattice<T,DESCRIPTOR>& block,
                            BlockIndicatorF<T,DESCRIPTOR::d>& boundaryIndicator,
                            BlockIndicatorF<T,DESCRIPTOR::d>& bulkIndicator,
                            AnalyticalF<DESCRIPTOR::d,T,T>& phi_d,
                            Cuboid<T,DESCRIPTOR::d>& cuboid)
{
  block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> solidLatticeR) {
    if (boundaryIndicator(solidLatticeR)) {
      for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
        Vector<T,DESCRIPTOR::d> boundaryLatticeR(solidLatticeR + descriptors::c<DESCRIPTOR>(iPop));
        if (block.isInside(boundaryLatticeR)) {
          const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
          auto x_b = block.get(boundaryLatticeR);
          const auto opp_bouzidi_dist = x_b.template getFieldComponent<descriptors::BOUZIDI_DISTANCE>(iPop_opposite);

          // check if distance from the fluid cell to the solid cell is a valid bouzidi distance
          if (opp_bouzidi_dist >= 0) {
            T boundaryPhysR[DESCRIPTOR::d] { };
            T phi_at_cell[1] { };
            cuboid.getPhysR(boundaryPhysR, boundaryLatticeR);

            phi_d(phi_at_cell, boundaryPhysR);

            // set computed velocity into the bouzidi velocity field
            block.get(boundaryLatticeR).template setFieldComponent<descriptors::BOUZIDI_ADE_DIRICHLET>(iPop_opposite, phi_at_cell[0]);
          }
        }
      }
    }
  });
}
}

#endif
