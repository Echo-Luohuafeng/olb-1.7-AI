/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2021 Davide Dapelo
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

/* advectionDiffusionReaction2d:
 * simulating a simple domain with no fluid motion and homogeneous
 * species concentration. Reproducing the chemical reactions:
 * (1) |a|A -> |b|B
 * (2) |c|B -> |d|D
 * (3) |e|E -> |f|F + |g|G
 * The (simplified) reaction rates are respectively:
 * nu1 = [A]/t0
 * nu2 - [C]/t0
 * nuE - [C]/t0
 * where t0 is a time conversion factor.
 * The initial conditions are respectively:
 * [A](t=0)=1; [B](t=0)=0
 * [D](t=0)=0
 * [E](t=0)=1; [F](t=0)=0; [G](t=0)=0
 * Analytical solution:
 * [A](t) = util::exp(-|a|*t/t0)
 * [B](t) = |b|/(|c|-|a|)*( util::exp(-|a|*t/t0) - util::exp(-|c|*t/t0) )
 * [D](t) = |bd|/(|c|-|a|)*( 1/|a|*(1 - util::exp(-|a|*t/t0)) - 1/|c|*(1 - util::exp(-|c|*t/t0)) )
 * [E](t) = util::exp(-|e|*t/t0 )
 * [F](t) = |f/e|*( 1 - util::exp(-|e|*t/t0) )
 * [G](t) = |g/e|*( 1 - util::exp(-|e|*t/t0) )
 */

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;
using namespace olb::graphics;

namespace olb {
namespace descriptors {
struct A_FIELD  : public FIELD_BASE<2,  0, 0> { };
struct A_SOURCE : public FIELD_BASE<1,  0, 0> { };
struct B_FIELD  : public FIELD_BASE<2,  0, 0> { };
struct B_SOURCE : public FIELD_BASE<1,  0, 0> { };
struct D_FIELD  : public FIELD_BASE<2,  0, 0> { };
struct D_SOURCE : public FIELD_BASE<1,  0, 0> { };
struct E_FIELD  : public FIELD_BASE<2,  0, 0> { };
struct E_SOURCE : public FIELD_BASE<1,  0, 0> { };
struct F_FIELD  : public FIELD_BASE<2,  0, 0> { };
struct F_SOURCE : public FIELD_BASE<1,  0, 0> { };
struct G_FIELD  : public FIELD_BASE<2,  0, 0> { };
struct G_SOURCE : public FIELD_BASE<1,  0, 0> { };
}
namespace names {
struct reactImplem1 : public descriptors::DESCRIPTOR_TAG { };
}
}

using T = FLOATING_POINT_TYPE;
typedef D2Q9<A_FIELD,A_SOURCE,B_FIELD,B_SOURCE,D_FIELD,D_SOURCE,E_FIELD,E_SOURCE,F_FIELD,F_SOURCE,G_FIELD,G_SOURCE,NORMAL_X,NORMAL_Y> DESCRIPTOR;

using FdParams = meta::list<
  fd::fdParams::Timestep,
  fd::fdParams::Diffusivity
>;

typedef fd::tag::UPWIND  SCHEME_ADV;
typedef fd::tag::CENTRAL SCHEME_DIFF;


///////////////////////////////////////////////////////////////////////////////////////////////////
// Functional for the analytical solution of single reactions of the type:
// aA -> bB + cC
// with rate A/t0
template <unsigned D, typename T, typename S>
class Sol_OneDecoupledEq final: public AnalyticalF<D,T,S> {
private:
  T _t, _t0, _a, _b, _c;
public:
  Sol_OneDecoupledEq(T t, T t0, T a, T b, T c=0) : AnalyticalF<D,T,S>(2), _t(t), _t0(t0), _a(a), _b(b), _c(c)
  {
    this->getName() = "sol_OneDecoupledEq ";
    if (t0 <= 0) throw std::invalid_argument("Time constant t0 must be strictly positive.");
    if (a  >= 0) throw std::invalid_argument("Stoichiometric coefficient a must be strictly negative.");
    if (b  <= 0) throw std::invalid_argument("Stoichiometric coefficient b must be strictly positive.");
    if (b  <  0) throw std::invalid_argument("Stoichiometric coefficient c must be positive or zero." );
  }
  bool operator() (T output[], const S x[]) override
  {
    output[0] = util::exp(_a*_t/_t0);
    output[1] = -_b/_a*(1. - util::exp(_a*_t/_t0));
    output[2] = -_c/_a*(1. - util::exp(_a*_t/_t0));
    return true;
  }
};


///////////////////////////////////////////////////////////////////////////////////////////////////
// Functional for the analytical solution of two coupled reactions of the type:
// aA -> bB
// cB -> qQ
// with rates A/t0 and B/t0
template <unsigned D, typename T, typename S>
class Sol_TwoCoupledEq final: public AnalyticalF<D,T,S> {
private:
  T _t, _t0, _a, _b, _c, _q;
public:
  Sol_TwoCoupledEq(T t, T t0, T a, T b, T c, T q) : AnalyticalF<D,T,S>(2), _t(t), _t0(t0), _a(a), _b(b), _c(c), _q(q)
  {
    this->getName() = "sol_TwoCoupledEq ";
    if (t0 <= 0) throw std::invalid_argument("Time constant t0 must be strictly positive.");
    if (a  >= 0) throw std::invalid_argument("Stoichiometric coefficient a must be strictly negative.");
    if (b  <= 0) throw std::invalid_argument("Stoichiometric coefficient b must be strictly positive.");
    if (c  >= 0) throw std::invalid_argument("Stoichiometric coefficient c must be strictly negative.");
    if (q  <= 0) throw std::invalid_argument("Stoichiometric coefficient q must be strictly positive.");
  }
  bool operator() (T output[], const S x[]) override
  {
    output[0] = util::exp(_a*_t/_t0);
    output[1] = _b/(_a - _c)*( util::exp(_a*_t/_t0) - util::exp(_c*_t/_t0) );
    output[2] = _b*_q/(_a - _c)*( 1/_c*(1 - util::exp(_c*_t/_t0)) - 1/_a*(1 - util::exp(_a*_t/_t0)) );
    return true;
  }
};


///////////////////////////////////////////////////////////////////////////////////////////////////
std::size_t iT   = 0; // global timestep

// (dimensionless) Parameters for the simulation setup
int nx   = 50;  // Number of internal lattice points along the x direction
int ny   = 50;   // Number of internal lattice points along the y direction
std::size_t tMax = 2000; // Total number of lattice updates
std::size_t tVtm = 100;  // Number of timesteps before producing output
std::size_t tGnu = 50; // Number of timesteps before producing Gnuplot output
T  stoichio_a    = -3.;  // A specie's stoichiometric coefficient, reaction 1. Negative because it is the reagent.
T  stoichio_b    =  2.;  // B specie's stoichiometric coefficient, reaction 1. Positive because it is the product.
T  stoichio_c    = -7.;  // B specie's stoichiometric coefficient, reaction 2. Negative because it is the reagent.
T  stoichio_d    =  5.;  // D specie's stoichiometric coefficient, reaction 2. Positive because it is the product.
T  stoichio_e    = -11.; // E specie's stoichiometric coefficient, reaction 3. Negative because it is the reagent.
T  stoichio_f    =  23.; // F specie's stoichiometric coefficient, reaction 3. Positive because it is the product.
T  stoichio_g    =  29.; // G specie's stoichiometric coefficient, reaction 3. Positive because it is the product.

// Dimensional parameters
T deltaX = 0.1;  // Lattice spacing (m)
T deltaT = 0.01; // Timestep (s)
T t0     = 10;  // time conversion factor for reaction rates


///////////////////////////////////////////////////////////////////////////////////////////////////
// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry<T,2>& superGeometry)
{
  /* MAT NUM | GEOMETRY     | FINITE-DIFFERENCE | lATTICE-BOLTZMANN
   * 1       | Bulk         | Bulk              | Bulk
   * 2       | Wall         | No-penetration    | Bulk
   */
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,2);
  superGeometry.rename(2,1,{1,1});

  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,2>& superGeometry )
{
  /* MAT NUM | GEOMETRY     | FINITE-DIFFERENCE | lATTICE-BOLTZMANN
   * 1       | Bulk         | Bulk              | Bulk
   * 2       | Wall         | No-penetration    | Bulk
   */
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  using AD_MODEL = FdAdvectionDiffusionModel<T, fd::AdvectionScheme<2,T,fd::tag::UPWIND>,
                                                fd::DiffusionScheme<2,T,fd::tag::CENTRAL>>;

  using REACTIONS = std::tuple<
    std::tuple<FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,A_SOURCE,A_FIELD>,
               FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,B_SOURCE,B_FIELD>>,
    std::tuple<FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,B_SOURCE,B_FIELD>,
               FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,D_SOURCE,D_FIELD>>,
    std::tuple<FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,E_SOURCE,E_FIELD>,
               FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,F_SOURCE,F_FIELD>,
               FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,G_SOURCE,G_FIELD>>
  >;

  ReactionGenerator2D<T,DESCRIPTOR,REACTIONS> reactionCoupling (
    {
      std::make_shared<ExpOn1stSpecieRate<T>>(converter.getLatticeTime(t0)),
      std::make_shared<ExpOn1stSpecieRate<T>>(converter.getLatticeTime(t0)),
      std::make_shared<ExpOn1stSpecieRate<T>>(converter.getLatticeTime(t0))
    },
    REACTIONS{ {{stoichio_a, iT}, {stoichio_b, iT}},
               {{stoichio_c, iT}, {stoichio_d, iT}},
               {{stoichio_e, iT}, {stoichio_f, iT}, {stoichio_g, iT}} }
  );
  sLattice.addLatticeCoupling(reactionCoupling, sLattice, sLattice, sLattice, sLattice, sLattice, sLattice, sLattice);
  setFdPostProcessor2D<T,DESCRIPTOR,AD_MODEL,FdParams,A_FIELD,A_SOURCE>(sLattice, superGeometry, 1);
  setFdPostProcessor2D<T,DESCRIPTOR,AD_MODEL,FdParams,B_FIELD,B_SOURCE>(sLattice, superGeometry, 1);
  setFdPostProcessor2D<T,DESCRIPTOR,AD_MODEL,FdParams,D_FIELD,D_SOURCE>(sLattice, superGeometry, 1);
  setFdPostProcessor2D<T,DESCRIPTOR,AD_MODEL,FdParams,E_FIELD,E_SOURCE>(sLattice, superGeometry, 1);
  setFdPostProcessor2D<T,DESCRIPTOR,AD_MODEL,FdParams,F_FIELD,F_SOURCE>(sLattice, superGeometry, 1);
  setFdPostProcessor2D<T,DESCRIPTOR,AD_MODEL,FdParams,G_FIELD,G_SOURCE>(sLattice, superGeometry, 1);

  setFdBoundary2D<T,DESCRIPTOR,AD_MODEL,fd::AdNeumannZeroBoundaryScheme<2,T,fd::tag::UPWIND>,FdParams,A_FIELD,A_SOURCE>(sLattice, superGeometry, 2);
  setFdBoundary2D<T,DESCRIPTOR,AD_MODEL,fd::AdNeumannZeroBoundaryScheme<2,T,fd::tag::UPWIND>,FdParams,B_FIELD,B_SOURCE>(sLattice, superGeometry, 2);
  setFdBoundary2D<T,DESCRIPTOR,AD_MODEL,fd::AdNeumannZeroBoundaryScheme<2,T,fd::tag::UPWIND>,FdParams,D_FIELD,D_SOURCE>(sLattice, superGeometry, 2);
  setFdBoundary2D<T,DESCRIPTOR,AD_MODEL,fd::AdNeumannZeroBoundaryScheme<2,T,fd::tag::UPWIND>,FdParams,E_FIELD,E_SOURCE>(sLattice, superGeometry, 2);
  setFdBoundary2D<T,DESCRIPTOR,AD_MODEL,fd::AdNeumannZeroBoundaryScheme<2,T,fd::tag::UPWIND>,FdParams,F_FIELD,F_SOURCE>(sLattice, superGeometry, 2);
  setFdBoundary2D<T,DESCRIPTOR,AD_MODEL,fd::AdNeumannZeroBoundaryScheme<2,T,fd::tag::UPWIND>,FdParams,G_FIELD,G_SOURCE>(sLattice, superGeometry, 2);

  sLattice.defineDynamics<BGKdynamics>(superGeometry.getMaterialIndicator({1,2}));
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  sLattice.setParameter<fd::fdParams::Timestep>(iT);
  sLattice.setParameter<fd::fdParams::Diffusivity>(T(0));

  auto& commFields = sLattice.getCommunicator(stage::PostPostProcess());
  commFields.requestField<A_FIELD>();
  commFields.requestField<A_SOURCE>();
  commFields.requestField<B_FIELD>();
  commFields.requestField<B_SOURCE>();
  commFields.requestField<D_FIELD>();
  commFields.requestField<D_SOURCE>();
  commFields.requestField<E_FIELD>();
  commFields.requestField<E_SOURCE>();
  commFields.requestField<F_FIELD>();
  commFields.requestField<F_SOURCE>();
  commFields.requestField<G_FIELD>();
  commFields.requestField<G_SOURCE>();
  commFields.requestOverlap(sLattice.getOverlap());
  commFields.exchangeRequests();
  clout << "Prepare Lattice ... OK" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Sets values at the boundary and external velocity field following linear Couette at iT=0
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR>& converter,
                        SuperGeometry<T,2>& superGeometry )
{
  /* MAT NUM | GEOMETRY     | FINITE-DIFFERENCE | lATTICE-BOLTZMANN
   * 1       | Bulk         | Bulk              | Bulk
   * 2       | Wall         | No-penetration    | Bulk
   */
  OstreamManager clout( std::cout,"setBoundaryValues" );

  AnalyticalConst<2,T,T> one1  {1.};
  AnalyticalConst<2,T,T> zero2 {0., 0.};
  AnalyticalConst<2,T,T> one2  {1., 1.};
  AnalyticalConst<2,T,T> u0 {T(), T()};

  sLattice.defineRhoU( superGeometry, 0, one1, u0 );
  sLattice.template defineField<A_FIELD> ( superGeometry.getMaterialIndicator({1,2}), one2  );
  sLattice.template defineField<A_SOURCE>( superGeometry.getMaterialIndicator({1,2}), zero2 );
  sLattice.template defineField<B_FIELD> ( superGeometry.getMaterialIndicator({1,2}), zero2 );
  sLattice.template defineField<B_SOURCE>( superGeometry.getMaterialIndicator({1,2}), zero2 );
  sLattice.template defineField<D_FIELD> ( superGeometry.getMaterialIndicator({1,2}), zero2 );
  sLattice.template defineField<D_SOURCE>( superGeometry.getMaterialIndicator({1,2}), zero2 );
  sLattice.template defineField<E_FIELD> ( superGeometry.getMaterialIndicator({1,2}), one2  );
  sLattice.template defineField<E_SOURCE>( superGeometry.getMaterialIndicator({1,2}), zero2 );
  sLattice.template defineField<F_FIELD> ( superGeometry.getMaterialIndicator({1,2}), zero2 );
  sLattice.template defineField<F_SOURCE>( superGeometry.getMaterialIndicator({1,2}), zero2 );
  sLattice.template defineField<G_FIELD> ( superGeometry.getMaterialIndicator({1,2}), zero2 );
  sLattice.template defineField<G_SOURCE>( superGeometry.getMaterialIndicator({1,2}), zero2 );

  sLattice.initialize();
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the results
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR>& converter,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer )
{
  OstreamManager clout( std::cout,"getResults" );

  static Gnuplot<T> gplot1_2("concentrations1_2");
  static Gnuplot<T> gplot3  ("concentrations3"  );

  Sol_TwoCoupledEq<2,T,int> sol1_2 (converter.getPhysTime(iT), t0, stoichio_a, stoichio_b, stoichio_c, stoichio_d);
  Sol_OneDecoupledEq<2,T,int> sol3 (converter.getPhysTime(iT), t0, stoichio_e, stoichio_f, stoichio_g);
  SuperLatticeExternal2D<T,DESCRIPTOR,A_FIELD> slA ( sLattice, iT );
  SuperLatticeExternal2D<T,DESCRIPTOR,B_FIELD> slB ( sLattice, iT );
  SuperLatticeExternal2D<T,DESCRIPTOR,D_FIELD> slD ( sLattice, iT );
  SuperLatticeExternal2D<T,DESCRIPTOR,E_FIELD> slE ( sLattice, iT );
  SuperLatticeExternal2D<T,DESCRIPTOR,F_FIELD> slF ( sLattice, iT );
  SuperLatticeExternal2D<T,DESCRIPTOR,G_FIELD> slG ( sLattice, iT );
  AnalyticalFfromSuperF2D<T> A( slA, true, 1 );
  AnalyticalFfromSuperF2D<T> B( slB, true, 1 );
  AnalyticalFfromSuperF2D<T> D( slD, true, 1 );
  AnalyticalFfromSuperF2D<T> E( slE, true, 1 );
  AnalyticalFfromSuperF2D<T> F( slF, true, 1 );
  AnalyticalFfromSuperF2D<T> G( slG, true, 1 );

  int point[] { (int)(nx/2)+1, (int)(ny/2)+1 };
  T pointP[] { point[0]*deltaX, point[1]*deltaX };
  T analytical1_2[] { T(), T(), T() };
  T analytical3[]   { T(), T(), T() };
  T numericalA[] { T() };
  T numericalB[] { T() };
  T numericalD[] { T() };
  T numericalE[] { T() };
  T numericalF[] { T() };
  T numericalG[] { T() };

  sol1_2(analytical1_2,point);
  sol3  (analytical3,  point);
  A(numericalA,pointP);
  B(numericalB,pointP);
  D(numericalD,pointP);
  E(numericalE,pointP);
  F(numericalF,pointP);
  G(numericalG,pointP);

  if (iT % tGnu ==0) {
    gplot1_2.setData(converter.getPhysTime(iT),
              { analytical1_2[0], numericalA[0], analytical1_2[1], numericalB[0], analytical1_2[2], numericalD[0] },
              { "A analytical", "A numerical", "B analytical", "B numerical", "D analytical", "D numerical" },
              "bottom right", { 'l', 'p', 'l', 'p', 'l', 'p' });
    gplot3.setData(converter.getPhysTime(iT),
              { analytical3[0], numericalE[0], analytical3[1], numericalF[0], analytical3[2], numericalG[0] },
              { "E analytical", "E numerical", "F analytical", "F numerical", "G analytical", "G numerical" },
              "bottom right", { 'l', 'p', 'l', 'p', 'l', 'p' });
    gplot1_2.writePNG();
    gplot3.writePNG();
  }

  if (iT % tVtm == 0) {
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime( iT ));
    clout << "A:  analytical=" << analytical1_2[0] << "  numerical=" << numericalA[0]<< "  error=" << util::abs(analytical1_2[0]-numericalA[0])/numericalA[0] << std::endl
          << "B:  analytical=" << analytical1_2[1] << "  numerical=" << numericalB[0]<< "  error=" << util::abs(analytical1_2[1]-numericalB[0])/numericalB[0] << std::endl
          << "D:  analytical=" << analytical1_2[2] << "  numerical=" << numericalD[0]<< "  error=" << util::abs(analytical1_2[2]-numericalD[0])/numericalD[0] << std::endl
          << "E:  analytical=" <<   analytical3[0] << "  numerical=" << numericalE[0]<< "  error=" << util::abs(  analytical3[0]-numericalE[0])/numericalE[0] << std::endl
          << "F:  analytical=" <<   analytical3[1] << "  numerical=" << numericalF[0]<< "  error=" << util::abs(  analytical3[1]-numericalF[0])/numericalF[0] << std::endl
          << "G:  analytical=" <<   analytical3[2] << "  numerical=" << numericalG[0]<< "  error=" << util::abs(  analytical3[2]-numericalG[0])/numericalG[0] << std::endl
          << std::endl;
  }

  if (iT == tMax) {
    gplot1_2.writePDF("concentrations1_2");
    gplot3.writePDF("concentrations3");
  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );

  OstreamManager clout( std::cout,"main" );

  UnitConverter<T,DESCRIPTOR> converter (
    (T) deltaX,    // physDeltaX
    (T) deltaT,    // physDeltaT
    (T) nx*deltaX, // charPhysLength
    (T) 0.,        // charPhysVelocity
    (T) 0.1,       // physViscosity
    (T) 1.         // physDensity
  );
  clout << "---------- Input data: ------------" << std::endl
        << "nx         = " << nx << std::endl
        << "ny         = " << ny << std::endl
        << "tMax       = " << tMax << std::endl
        << "tVtm       = " << tVtm << std::endl
        << "deltaX     = " << deltaX << std::endl
        << "deltaT     = " << deltaT << std::endl
        << "stoichio_a = " << stoichio_a << std::endl
        << "stoichio_b = " << stoichio_b << std::endl
        << "stoichio_c = " << stoichio_c << std::endl
        << "stoichio_d = " << stoichio_d << std::endl
        << "stoichio_e = " << stoichio_e << std::endl
        << "stoichio_f = " << stoichio_f << std::endl
        << "stoichio_g = " << stoichio_g << std::endl
        << "t0         = " << t0 << std::endl
        << "------------------------------------" << std::endl;
  converter.print();

  /// === 2rd Step: Prepare Geometry ===
  std::vector<T> origin { 0., 0. };
  std::vector<T> extend { (nx-1)*deltaX, (ny-1)*deltaX };
  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 3*singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry2D<T> cuboidGeometry(cuboid, deltaX, noOfCuboids);
  cuboidGeometry.setPeriodicity(false, false);
  cuboidGeometry.print();

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice( superGeometry );
  prepareLattice( sLattice, converter, superGeometry );

  // === 4th Step: Definition of Initial and Boundary Conditions ===
  setBoundaryValues( sLattice, converter, superGeometry );

  // === 5th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( tMax, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( iT = 0; iT <= tMax; ++iT ) {
    sLattice.setParameter<fd::fdParams::Timestep>(iT);

    // === 6th Step: Collide and Stream Execution ===
    sLattice.executeCoupling();
    sLattice.collideAndStream();

    // === 8th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, superGeometry, timer );
  }

  timer.stop();
  timer.printSummary();
}
