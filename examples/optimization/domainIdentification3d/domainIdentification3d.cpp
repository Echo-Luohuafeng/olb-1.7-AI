/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020-22 Fabian Klemens, Mathias J. Krause, Benjamin
 *  Förster, Julius Jeßberger
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
 * In this example, a domain identification problem is solved with the help of
 * adjoint LBM. Cf. https://doi.org/10.1016/j.camwa.2018.07.010 for details.
 */

#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;
using namespace olb::opti;

using Descriptor = DualPorousD3Q19Descriptor;
using Lattices = meta::map<
                 NavierStokes, Descriptor
                 >;


/// Collect parameters which are relevant for optimization.
// Inherit/ extend DistributedOptiSimulation which has the adjoint-optimization-
// specific data.
template <typename T, SolverMode MODE>
struct TestDomainOptiParameters
 : public parameters::DistributedOptiSimulation<T,Lattices,MODE>
{
  std::shared_ptr<IndicatorF3D<T>>      designDomain;
  std::shared_ptr<IndicatorF3D<T>>      objectiveDomain;
};

/// Enable reading the parameters in DistributedOptiSimulation from xml
template<typename T, SolverMode MODE, typename TAG>
struct parameters::Reader<TestDomainOptiParameters<T,MODE>, TAG>
 : public parameters::ReaderBase<TestDomainOptiParameters<T,MODE>>
{
  using ReaderBase<TestDomainOptiParameters<T,MODE>>::ReaderBase;

  void read(XMLreader const& xml)
  {
    using namespace parameters;
    Reader<DistributedOptiSimulation<T,Lattices,MODE>, TAG>(this->params).read(xml);
    this->params->objectiveDomain = createIndicatorF3D<T> (xml["ObjectiveDomain"], false);
    this->params->designDomain = createIndicatorF3D<T> (xml["DesignDomain"]);
  }
};

/// Parameters which are relevant for simulation
template <typename T>
struct TestDomainSimulationParameters
 : public parameters::XmlSimulation<T,Lattices>
{
  std::shared_ptr<IndicatorF3D<T>> simulationObject;
};

/// XML interface for TestDomainSimulationParameters
template<typename T, typename TAG>
struct parameters::Reader<TestDomainSimulationParameters<T>, TAG>
 : public parameters::ReaderBase<TestDomainSimulationParameters<T>>
{
  using ReaderBase<TestDomainSimulationParameters<T>>::ReaderBase;

  void read(XMLreader const& xml) {
    using namespace parameters;
    Reader<XmlSimulation<T,Lattices>, TAG>(this->params).read(xml);
    this->params->simulationObject = createIndicatorF3D<T> (xml["SimulationObject"]);
  }
};



/// Pack all parameter structs together in a map
template<typename T, opti::SolverMode MODE>
using Parameters = meta::map<
  Opti,             TestDomainOptiParameters<T,MODE>,
  Simulation,       TestDomainSimulationParameters<T>,
  Stationarity,     parameters::Stationarity<T>,
  OutputOpti,       parameters::OptiOutput<T,MODE>,
  Results,          parameters::DistributedOptiSimulationResults<T,Lattices,MODE>,
  Output,           parameters::OutputGeneral<T>,
  VisualizationVTK, parameters::OutputPlot<T>
>;

/// Solver class implements the simulation
// Inherit from AdjointLbSolverBase since most of the adjoint-specific details
// are implemented there
template<typename T, SolverMode MODE>
class DomainIdSolver : public AdjointLbSolverBase<T,Parameters<T,MODE>,Lattices,MODE>
{
private:
  mutable OstreamManager                clout {std::cout, "DomainIdSolver"};

public:
  using LATTICES = Lattices;
  using descriptor = typename Lattices::values_t::template get<0>;

  DomainIdSolver(utilities::TypeIndexedSharedPtrTuple<Parameters<T,MODE>> params)
   : DomainIdSolver::AdjointLbSolverBase(params)
  {
    const std::string name = (MODE == SolverMode::Reference) ? "Reference_Solution"
      : ((MODE == SolverMode::Primal) ? "Flow_Simulations"
      : "Dual_Simulations");
    this->parameters(VisualizationVTK()).filename = name;
  }

protected:
  /// Set up a cubic geometry and mark the design domain inside
  void prepareGeometry() override
  {
    auto& xml = *(this->parameters(Simulation()).xml);
    // Read mesh geometry from XML file
    std::shared_ptr<IndicatorF3D<T>> cube
                                  = createIndicatorF3D<T>(xml["Application"]["Mesh"], false);
    // Create CuboidGeometry3D from mesh indicator
    this->_cGeometry = std::make_shared<CuboidGeometry3D<T>>(
                         *cube,
                         this->converter().getPhysDeltaX(),
                         singleton::mpi().getSize() * this->parameters(Simulation()).noC);
    // Create HeuristicLoadBalancer from CuboidGeometry3D
    this->_loadBalancer = std::make_shared<HeuristicLoadBalancer<T>>(*this->_cGeometry);
    // Create SuperGeometry
    this->_sGeometry = std::make_shared<SuperGeometry<T,3>> (
                         *this->_cGeometry,
                         *this->_loadBalancer,
                         this->parameters(Simulation()).overlap);

    // Material 2 everywhere
    this->geometry().rename(0,2);
    // Material 1 inside of border
    this->geometry().rename(2,1, {1,1,1});

    // Outer boundaries
    Vector<T,3> origin = {0.,0.,0.};
    origin[0] += this->converter().getPhysDeltaX()/2.;
    origin[1] -= this->converter().getPhysDeltaX()/2.;
    origin[2] += this->converter().getPhysDeltaX()/2.;

    Vector<T,3> extend = {1.,0.,1.};
    extend[0] -= 2.*this->converter().getPhysDeltaX()/2.;
    extend[1] +=    this->converter().getPhysDeltaX()/2.;
    extend[2] -= 2.*this->converter().getPhysDeltaX()/2.;

    IndicatorCuboid3D<T> inflow( extend,origin );
    this->geometry().rename( 2,3,inflow );

    origin = {0.,1.,0.};
    origin[0] += this->converter().getPhysDeltaX()/2.;
    origin[1] -= this->converter().getPhysDeltaX()/2.;
    origin[2] += this->converter().getPhysDeltaX()/2.;

    extend = {1.,1.,1.};
    extend[0] -= 2.*this->converter().getPhysDeltaX()/2.;
    extend[1] +=    this->converter().getPhysDeltaX()/2.;
    extend[2] -= 2.*this->converter().getPhysDeltaX()/2.;

    IndicatorCuboid3D<T> outflow(extend, origin);
    this->geometry().rename( 2,4,outflow );

    // Indicators for material 6 (designDomain)
    this->geometry().rename(1, 6, this->parameters(Opti()).designDomain);

    if constexpr (MODE == SolverMode::Reference) {
      const unsigned nSimulation = this->geometry().getStatistics().getNvoxel();
      clout << "Number of fluid voxels in simulation domain: " << nSimulation << std::endl;
      const unsigned nDesign = this->geometry().getStatistics().getNvoxel(6);
      clout << "Number of voxels in design domain: " << nDesign << std::endl;
    }
  }

  void prepareLattices() override
  {
    auto& lattice = this->lattice(NavierStokes());
    const T omega = this->converter().getLatticeRelaxationFrequency();

    lattice.template defineDynamics<NoDynamics<T,descriptor>>(this->geometry(), 0);

    // PRIMAL
    if constexpr ( (MODE == SolverMode::Reference) || (MODE == SolverMode::Primal) ) {
      auto bulkIndicator = this->geometry().getMaterialIndicator({1, 3, 4, 6});
      lattice.template defineDynamics<PorousBGKdynamics>(bulkIndicator);

      setBounceBackBoundary(lattice, this->geometry(), 2);

      setLocalVelocityBoundary<T,descriptor>(
        lattice,
        omega,
        this->geometry().getMaterialIndicator({3, 4}));

      lattice.template setParameter<descriptors::OMEGA>(omega);
    }
    else {   // DUAL
      auto bulkIndicator = this->geometry().getMaterialIndicator({1, 6});
      lattice.defineDynamics(bulkIndicator,
                             new DualPorousBGKdynamics<T,descriptor> (omega));

      auto bounceBackIndicator = this->geometry().getMaterialIndicator({2, 3, 4});
      setBounceBackBoundary(lattice, bounceBackIndicator);
    }
  }

  void setInitialValues() override
  {
    std::unique_ptr<SuperIndicatorF3D<T>> bulkIndicator;
    if constexpr ( (MODE == SolverMode::Reference) || (MODE == SolverMode::Primal) ) {
      bulkIndicator = this->geometry().getMaterialIndicator({1, 3, 4});
    }
    else {
      bulkIndicator = this->geometry().getMaterialIndicator({1, 2, 3, 4});
    }

    AnalyticalConst3D<T,T> rhoF(1);
    std::vector<T> velocity(3,T());
    AnalyticalConst3D<T,T> uF(velocity);

    this->lattice(NavierStokes()).defineRhoU( bulkIndicator, rhoF, uF );
    this->lattice(NavierStokes()).iniEquilibrium( bulkIndicator, rhoF, uF);


    AnalyticalConst3D<T,T> one(1.);

    for ( int i: {
            0,1,2,3,4
          } ) {
      this->lattice(NavierStokes()).template defineField<POROSITY>(this->geometry(), i, one);
    }

    if constexpr ( MODE == SolverMode::Reference ) {
      // solid cube is the middle
      AnalyticalConst3D<T,T> zero(0.);
      this->lattice(NavierStokes()).template defineField<POROSITY>(
        this->geometry(),
        *this->parameters(Simulation()).simulationObject,
        zero);
    } else {
      this->lattice(NavierStokes()).template defineField<POROSITY>(
        this->geometry().getMaterialIndicator({6}),
        *(this->parameters(Opti()).controlledField));
    }

    if constexpr ( MODE == SolverMode::Dual ) {
      this->loadPrimalPopulations();
    }
  }

  void setBoundaryValues(std::size_t iT) override
  {
    const std::size_t itStart = this->converter().getLatticeTime(this->parameters(Simulation()).startUpTime);

    if (iT <= itStart) {
      PolynomialStartScale<T,std::size_t> StartScale( itStart, T( 1 ) );
      std::size_t iTvec[1] = {iT};
      T frac[1] = {};
      StartScale( frac,iTvec );

      if constexpr ( (MODE == SolverMode::Reference) || (MODE == SolverMode::Primal) ) {
        AnalyticalConst3D<T,T> uF(0., frac[0] * this->converter().getCharLatticeVelocity(), 0.);
        //this->lattice(NavierStokes()).defineU(this->geometry(), 2, uF);
        this->lattice(NavierStokes()).defineU(this->geometry(), 3, uF);
        this->lattice(NavierStokes()).defineU(this->geometry(), 4, uF);
      }
    }
  }

  void prepareVTK() const override
  {
    if constexpr (MODE == SolverMode::Reference) {
      DomainIdSolver::LbSolver::prepareVTK();
    }
  }

  void writeVTK(std::size_t iT) const override
  {
    SuperVTMwriter3D<T> writer(this->parameters(VisualizationVTK()).filename);
    if constexpr (MODE == SolverMode::Reference) {
      // ARTIFICIAL DATA: write solution (unsteady, until convergence in time is obtained)

      SuperLatticeGeometry3D<T,descriptor> geometry(this->lattice(NavierStokes()), this->geometry());
      SuperLatticePhysVelocity3D<T,descriptor> velocity(this->lattice(NavierStokes()), this->converter());
      SuperLatticePhysPressure3D<T,descriptor> pressure(this->lattice(NavierStokes()), this->converter());
      SuperLatticePorosity3D<T,descriptor> porosity(this->lattice(NavierStokes()));
      writer.addFunctor(porosity);
      writer.addFunctor(geometry);
      writer.addFunctor(velocity);
      writer.addFunctor(pressure);
      writer.write(iT);
    }
  }

  /// Compute and store data at the end of the simulation
  void computeResults() override
  {
    auto& lattice = this->lattice(NavierStokes());

    // Store the full lattice for later functor evaluation
    this->parameters(Results()).lattice = std::get<0>(this->_sLattices);

    // Store the full geometry for later functor evaluation
    if constexpr ((MODE == SolverMode::Reference) || (MODE == SolverMode::Primal))
    {
      this->parameters(Results()).geometry = this->_sGeometry;
    }

    if constexpr (MODE == SolverMode::Primal) {
      // compute objective
      int tmp[1];
      T output[1];
      const std::shared_ptr<SuperIndicatorF3D<T>> objectiveDomain
       = std::make_shared<SuperIndicatorFfromIndicatorF3D<T>> (
        this->parameters(Opti()).objectiveDomain,
        this->geometry());

      std::shared_ptr<SuperF3D<T,T>> physLatticeVelocity
                                  = std::make_shared<SuperLatticePhysVelocity3D<T,descriptor>>(
                                      lattice, this->converter());
      RelativeDifferenceObjective3D<T,descriptor> objective(
        lattice,
        physLatticeVelocity,
        this->parameters(Opti()).referenceSolution,
        objectiveDomain);

      objective(output,tmp);
      this->parameters(Results()).objective = output[0];
      this->parameters(Results()).objectiveComputed = true;

      // set population
      this->parameters(Results()).fpop = std::make_shared<SuperLatticeFpop3D<T,descriptor>>(
        lattice);

      // set population derivative
      std::shared_ptr<SuperF3D<T,T>> dPhysVelocityf = std::make_shared<SuperLatticeDphysVelocityDf3D<T,descriptor>> (
        lattice, this->converter());
      this->parameters(Results()).djdf = std::make_shared<DrelativeDifferenceObjectiveDf3D_Lattice<T,descriptor>>(
        lattice,
        physLatticeVelocity,
        dPhysVelocityf,
        this->parameters(Opti()).referenceSolution,
        objectiveDomain);

      // set control derivative
      std::shared_ptr<AnalyticalF<3,T,T>> zero = std::make_shared<AnalyticalConst<3,T,T>> (Vector<T,3>(0));
      this->parameters(Results()).djdalpha = std::make_shared<SuperLatticeFfromAnalyticalF3D<T,descriptor>> (zero, lattice);
    }
  }

public:
  /// Optimization: write stationary flow field once after each simulation
  // This method is called by the optimizer once at each successful optimization step
  // (only for the primal simulation)
  void postProcess() override
  {
    if constexpr (MODE == SolverMode::Primal) {
      if (this->parameters(names::VisualizationVTK()).output != "off") {
        if (this->parameters(OutputOpti()).counterOptiStep == 0) {
          DomainIdSolver::LbSolver::prepareVTK();
        }

        SuperVTMwriter3D<T> writer(this->parameters(VisualizationVTK()).filename);
        SuperLatticeGeometry3D<T,descriptor> geometry(this->lattice(NavierStokes()), this->geometry());
        SuperLatticePhysVelocity3D<T,descriptor> velocity(this->lattice(NavierStokes()), this->converter());
        SuperLatticePhysPressure3D<T,descriptor> pressure(this->lattice(NavierStokes()), this->converter());
        SuperLatticePorosity3D<T,descriptor> porosity(this->lattice(NavierStokes()));
        writer.addFunctor(porosity);
        writer.addFunctor(geometry);
        writer.addFunctor(velocity);
        writer.addFunctor(pressure);

        writer.write(this->parameters(OutputOpti()).counterOptiStep);
      }
    }
  }
};


int main(int argc, char **argv)
{
  olbInit(&argc, &argv);
  XMLreader config("parameter.xml");
  OstreamManager clout("main");

  using T = FLOATING_POINT_TYPE;

  OptiCaseDual<T,DomainIdSolver> optiCase(config);

  auto optimizer = createOptimizerLBFGS<T>(config, optiCase._dimCtrl);

  T startValue;
  config.readOrWarn<T>("Optimization", "StartValue", "", startValue);
  startValue = projection::getInitialControl(startValue, optiCase);
  optimizer->setStartValue(startValue);

  optimizer->optimize(optiCase);
}
