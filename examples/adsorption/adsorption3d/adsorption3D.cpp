/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2022 Florian Raichle, Fedor Bukreev
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

/**
 * adsorption3D
 * Simulating adsorption in a batch reactor using an Euler-Euler approach.
 * The model is based on the linear driving force model and uses advection diffusion reaction lattices for particles,
 * solute and particle loading.
 *
 * Different isotherms and mass transfer models can be used.
 * An analytical solution is implemented when using the linear isotherm and surface diffusion.
 */

#include "olb3D.h"
#include "olb3D.hh"   // use only generic version!
#include <vector>
#include <iostream>

#include "../isotherms.h"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;

typedef double T;
typedef D3Q19<VELOCITY> NSDESCRIPTOR;
typedef D3Q7<VELOCITY, SOURCE> ADEDESCRIPTOR;

using bulkDynamics = BGKdynamics<T, NSDESCRIPTOR>;
using bulkDynamicsAD = SourcedAdvectionDiffusionBGKdynamics<T, ADEDESCRIPTOR>;

const int dim = 3;
T cubeLength = 0.1;

const int N = 21;          // resolution of the model
const T Re = 17;           // Reynolds number
const T Sc = 1.5;
const T particleRadius = 1.5e-04;  // particles radius
const T partRho = 0.940;   // particles density
const T c_0 = 1;
const T tau_ads = 0.6952380952380952;
const T fluidViscosity = 0.00011764705882352942;
const T maxPhysT = 15;

const T isoConstA = 45;
const T isoConstB = 0.5;
const T k_f = 5.37E-5;  //5.*7E-11/particleRadius;
const T D_s = 5E-11;

template <typename T>
class BatchSolution: public AnalyticalF3D<T, T> {
protected:
    T t;
    T D_b_;

public:
    BatchSolution(T time, T D_b, T ks) : AnalyticalF3D<T, T>(1), t(/*time * ks*/time*3.*k_f/particleRadius/D_b), D_b_(D_b) {}
    bool operator () (T output[], const T input[]) override {
      output[0] = T(1)/(D_b_+T(1)) + D_b_ / (D_b_+T(1)) * util::exp(-(D_b_+T(1))*t);
      return true;
    }
};

enum material {
  buffer = 0,
  fluid = 1,
  boundary = 2
};

void prepareGeometry(UnitConverter<T, NSDESCRIPTOR> const &converter, SuperGeometry<T, dim> &superGeometry) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(buffer, fluid);

  superGeometry.communicate();

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLatticeNS(SuperLattice<T, NSDESCRIPTOR> &NSLattice,
                      UnitConverter<T, NSDESCRIPTOR> const &converter,
                      SuperGeometry<T, dim> &superGeometry,
                      T omega) {
  OstreamManager clout(std::cout, "prepareLatticeNS");
  clout << "Prepare NSE Lattice ..." << std::endl;

  // Material=1 --> bulk dynamics
  NSLattice.defineDynamics<bulkDynamics>(superGeometry, 1);

  // Initial conditions
  AnalyticalConst3D<T, T> rhoFluid(1.);
  AnalyticalConst3D<T, T> u0(0.0, 0.0, 0.0);

  auto bulkIndicator = superGeometry.getMaterialIndicator({0, 1});

  // Initialize all values of distribution functions to their local equilibrium
  NSLattice.defineRhoU(bulkIndicator, rhoFluid, u0);
  NSLattice.iniEquilibrium(bulkIndicator, rhoFluid, u0);

  // Lattice initialize
  NSLattice.setParameter<descriptors::OMEGA>(omega);
  NSLattice.initialize();

  {
    auto &communicator = NSLattice.getCommunicator(stage::Full());
    communicator.requestField<descriptors::VELOCITY>();
    communicator.requestOverlap(NSLattice.getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Prepare NSE Lattice ... OK" << std::endl;
}

void prepareLatticeAD(SuperLattice<T, ADEDESCRIPTOR>*& ADLattice,
                      UnitConverter<T, NSDESCRIPTOR> const &converter,
                      SuperGeometry<T, dim> &superGeometry,
                      T rho_,
                      T omegaAD) {
  OstreamManager clout(std::cout, "prepareLatticeAD");
  clout << "Prepare ADE Lattice ..." << std::endl;

  // Material=1 --> bulk dynamics
  ADLattice->defineDynamics<bulkDynamicsAD>(superGeometry, 1);

  // Initial conditions
  AnalyticalConst3D<T, T> rho(rho_);
  AnalyticalConst3D<T, T> u0(0.0, 0.0, 0.0);

  auto bulkIndicator = superGeometry.getMaterialIndicator({0, 1});

  // Initialize all values of distribution functions to their local equilibrium
  ADLattice->defineRho(bulkIndicator, rho);
  ADLattice->iniEquilibrium(bulkIndicator, rho, u0);

  // Lattice initialize
  ADLattice->setParameter<descriptors::OMEGA>(omegaAD);
  ADLattice->initialize();

  {
    auto& communicator = ADLattice->getCommunicator(stage::Full());
    communicator.requestField<descriptors::VELOCITY>();
    communicator.requestField<descriptors::SOURCE>();
    communicator.requestOverlap(ADLattice->getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Prepare ADE Lattice ... OK" << std::endl;
}

void getResults(SuperLattice<T, NSDESCRIPTOR> &NSLattice,
                SuperLattice<T, ADEDESCRIPTOR> &PADLattice,
                SuperLattice<T, ADEDESCRIPTOR> &QADLattice,
                SuperLattice<T, ADEDESCRIPTOR> &CADLattice,
                UnitConverter<T, NSDESCRIPTOR> const &converter,
                UnitConverter<T, ADEDESCRIPTOR> const &adsConverter,
                int iT,
                SuperGeometry<T, dim> &superGeometry,
                Timer<T> &timer,
                T k_s,
                T D_b,
                Isotherm<T> const &isotherm) {
  OstreamManager clout(std::cout, "getResults");

  BatchSolution<T> concentrationSol(adsConverter.getPhysTime(iT), D_b, k_s);
  SuperLatticeFfromAnalyticalF3D<T, ADEDESCRIPTOR> solution(concentrationSol, QADLattice);

  SuperVTMwriter3D<T> vtmWriter("adsorption3D");
  SuperLatticeGeometry3D<T, NSDESCRIPTOR> materials(NSLattice, superGeometry);
  SuperLatticeDensity3D<T, ADEDESCRIPTOR> particles(PADLattice);
  SuperLatticeDensity3D<T, ADEDESCRIPTOR> loading(QADLattice);
  SuperLatticeDensity3D<T, ADEDESCRIPTOR> phosphateConcentration(CADLattice);
  SuperLatticePhysField3D<T, ADEDESCRIPTOR, descriptors::VELOCITY> extFieldParticles(PADLattice, converter.getConversionFactorVelocity());
  SuperLatticePhysField3D<T, ADEDESCRIPTOR, descriptors::VELOCITY> extFieldPhosphate(CADLattice, converter.getConversionFactorVelocity());
  SuperLatticePhysField3D<T, ADEDESCRIPTOR, descriptors::SOURCE> sourcePhosphate(CADLattice, 1);
  AnalyticalFfromSuperF3D<T> concentrationInterpolation(phosphateConcentration, true, true);
  AnalyticalFfromSuperF3D<T> loadingInterpolation(loading, true, true);

  vtmWriter.addFunctor(materials);
  vtmWriter.addFunctor(particles, "particle density");
  vtmWriter.addFunctor(extFieldParticles, "particle velocity");
  vtmWriter.addFunctor(loading, "particle loading");
  vtmWriter.addFunctor(phosphateConcentration, "solute concentration");
  vtmWriter.addFunctor(sourcePhosphate, "solute source");
  vtmWriter.addFunctor(solution, "solution");

  const int vtkIter = adsConverter.getLatticeTime(maxPhysT/T(100))+1;

  if (iT == 0) {
    SuperLatticeGeometry3D<T, NSDESCRIPTOR> geometry(NSLattice, superGeometry);
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid(NSLattice);
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank(NSLattice);
    vtmWriter.write(geometry);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();

    // Print some output of the chosen simulation setup
    clout << "N=" << N << "; maxTimeSteps=" << converter.getLatticeTime(maxPhysT)
          << "; noOfCuboid=" << superGeometry.getCuboidGeometry().getNc() << "; Re=" << Re
          << "; St="
          << (T(2) * partRho * particleRadius * particleRadius * converter.getCharPhysVelocity())
              / (T(9) * converter.getPhysViscosity() * converter.getPhysDensity() * converter.getCharPhysLength())
          << std::endl;
  }

  if (iT % vtkIter == 0) {
    // Writes the vtk files
    vtmWriter.write(iT);

    auto indicatorF = superGeometry.getMaterialIndicator(1);
    SuperRelativeErrorL2Norm3D<T> relConcentrationError(phosphateConcentration, concentrationSol, indicatorF);

    // Writes output on the console
    timer.update(iT);
    timer.printStep();
    NSLattice.getStatistics().print(iT, converter.getPhysTime(iT)/*iT * converter.getCharLatticeVelocity() / T(converter.getResolution())*/);
  }
}

int main(int argc, char *argv[]) {
  // === 1st Step: Initialization ===
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager
  clout(std::cout, "main");

 AdsorptionConverterFromSchmidtNumberAndRelaxation<T, ADEDESCRIPTOR> const adsConverter(
   (T) Sc,
   (T) Re,
   (T) tau_ads,
   (T) cubeLength,
   (T) cubeLength,
   (T) N,
   (T) 0.02,
   (T) 1,
   (T) partRho
 );

 UnitConverter<T, NSDESCRIPTOR> const converter(
   adsConverter.getPhysDeltaX(),
   adsConverter.getPhysDeltaT(),
   cubeLength,
   Re * fluidViscosity / (cubeLength),
   fluidViscosity,
   1.225
   );

  time_t now = time(nullptr);
  char buffer [80];
  strftime (buffer,80,"%F %H-%M", localtime(&now));

  // Prints the converter log as console output
  converter.print();
  adsConverter.print();
  // Writes the converter log in a file
  converter.write("adsorption3D");
  adsConverter.write("phosphate");

  // compute relaxation parameter to solve the advection-diffusion equation in the lattice Boltzmann context
  T omegaAD = adsConverter.getLatticeRelaxationFrequency();
  T omega = converter.getLatticeRelaxationFrequency();

  // === 2rd Step: Prepare Geometry ===
  Vector<T, dim> extend(cubeLength, cubeLength, cubeLength);
  Vector<T, dim> origin
      (-converter.getCharPhysLength() / 2, -converter.getCharPhysLength() / 2, -converter.getCharPhysLength() / 2);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry3D<T> cuboidGeometry(cuboid, converter.getConversionFactorLength(), noOfCuboids);
  cuboidGeometry.setPeriodicity(true, true, true);

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  // Instantiation of a superGeometry
  SuperGeometry<T, dim> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(converter, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, NSDESCRIPTOR> NSLattice(superGeometry);
  SuperLattice<T, ADEDESCRIPTOR> PADLattice(superGeometry);
  SuperLattice<T, ADEDESCRIPTOR> QADLattice(superGeometry);
  SuperLattice<T, ADEDESCRIPTOR> CADLattice(superGeometry);

  //prepareLattice and setBoundaryConditions
  T rho[3] = {0.};
  rho[0] = 1; //particle density at the start
  rho[1] = c_0; //concentration at the start
  rho[2] = 0.; //adsorbed load at the start

  std::vector<SuperLattice < T, ADEDESCRIPTOR>*> partners;
  partners.emplace_back(&PADLattice);
  partners.emplace_back(&CADLattice);
  partners.emplace_back(&QADLattice);

  prepareLatticeNS(NSLattice,
                   converter,
                   superGeometry,
                   omega);

  for(int i = 0; i<3; i++){
    prepareLatticeAD(partners[i],
                   converter,
                   superGeometry,
                   rho[i],
                   omegaAD);
    }

  LinearIsotherm<T> const isotherm(isoConstA, isoConstB);
  T V_l = cubeLength * cubeLength * cubeLength;
  T D_b = V_l*partRho * isotherm.getLoading(c_0)/(V_l*c_0);
  AdsorptionReaction<T, ADEDESCRIPTOR> reaction(adsConverter, &isotherm, k_f, D_s, c_0, particleRadius);
  reaction.write("adsorption3D");
  reaction.print(clout);
  AdsorptionFullCouplingPostProcessorGenerator3D<T, NSDESCRIPTOR, ADEDESCRIPTOR> adsorptionCoupling(&reaction);

  AdvDiffDragForce3D<T, NSDESCRIPTOR, ADEDESCRIPTOR> dragForce(converter, particleRadius, partRho);
  adsorptionCoupling.addForce(dragForce);
  NSLattice.addLatticeCoupling(adsorptionCoupling, partners);

  // === 4th Step: Main Loop with Timer ===
  Timer<T> timer(adsConverter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel());
  timer.start();
  for (std::size_t iT = 0; iT < adsConverter.getLatticeTime(maxPhysT); ++iT) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    SuperLatticeVelocity3D<T,NSDESCRIPTOR> velocity(NSLattice);
    NSLattice.defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator({1, 2}), velocity);
    // === 6th Step: Computation and Output of the Results ===
    getResults(NSLattice, PADLattice, QADLattice, CADLattice, converter, adsConverter, iT, superGeometry, timer, reaction.getPhysSurfaceTransferConstant(), D_b, isotherm);

    // === 7th Step: Collide and Stream Execution ===
    NSLattice.executeCoupling();
    for(int i = 0; i<3; i++) {
      partners[i]->executeCoupling();
      }

    for(int i = 0; i<3; i++) {
      partners[i]->collideAndStream();
      }
    NSLattice.collideAndStream();
  }

  timer.stop();
  timer.printSummary();
}
