/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2021 Mathias J. Krause, Fabian Klemens,
 *  Julius Jessberger
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

#ifndef TEST_FLOW_OPTI_H
#define TEST_FLOW_OPTI_H


/** \file A fluid flow simulation test case with analytical solution - implementation
 * of parameter identification using an optimization approach,
 * cf. https://publikationen.bibliothek.kit.edu/1000019768.
 */

// There, the standard setup of this simulation is implemented
#include "../../laminar/testFlow3dSolver/testFlow3dSolver.h"

using DualLattices = meta::map<
  NavierStokes, descriptors::DualForcedD3Q19Descriptor
>;

enum TestFlowOptiMode {VELOCITY, DISSIPATION};
enum OptiReferenceSolution {ANALYTICAL, DISCRETE};


/// Adds objective computation to TestFlowBase solver implementation
template<
  typename T,
  typename PARAMETERS,
  typename LATTICES,
  SolverMode MODE=SolverMode::Reference>
class TestFlowOptiBase : public TestFlowBase<T,PARAMETERS,LATTICES,MODE>
{
private:
  mutable OstreamManager            clout {std::cout, "TestFlowOptiBase"};

public:
  TestFlowOptiBase(utilities::TypeIndexedSharedPtrTuple<PARAMETERS> params)
   : TestFlowOptiBase::TestFlowBase(params)
  { }

  using typename TestFlowOptiBase::TestFlowBase::descriptor;

protected:
  /** compute velocity or dissipation L2-error
   * objective = 0.5 * error^2
   */
  T computeObjective() const
  {
    int tmp[1];
    T output[1];

    using BT = BaseType<T>;
    std::shared_ptr<AnalyticalF<3,BT,BT>> help;
    std::shared_ptr<AnalyticalF<3,T,T>> analytical;
    std::shared_ptr<SuperF3D<T,T>>      simulated;

    switch (this->parameters(Opti()).optiMode)
    {
      case DISSIPATION:
      {
        // assign analytical and simulated dissipation for error computation
        switch (this->parameters(Opti()).optiReferenceMode)
        {
          case ANALYTICAL:
          {
            analytical = std::make_shared<DissipationTestFlow3D<T,T,descriptor>>(
              this->converter());
          } break;

          case DISCRETE:
          {
            // user has to guarantee that reference solution is computed
            // cf. node Optimization -> ReferenceSolution in xml file
            help = std::make_shared<AnalyticalFfromSuperF3D<BT>>(
              *(this->parameters(Opti()).referenceSolution));
            analytical = std::make_shared<AnalyticalTypecast<3,BT,T,BT,T>>(*help);
          } break;
        }

        simulated = std::make_shared<SuperLatticePhysDissipation3D<T,descriptor>>(
          this->lattice(), this->converter());
      } break;

      case VELOCITY:
      {
        // assign analytical and simulated velocity field for error computation
        switch (this->parameters(Opti()).optiReferenceMode)
        {
          case ANALYTICAL:
          {
            analytical = this->_velocity;
          } break;

          case DISCRETE:
          {
            // user has to guarantee that reference solution is computed
            // cf. node Optimization -> ReferenceSolution in xml file
            help = std::make_shared<AnalyticalFfromSuperF3D<BT>>(
              *(this->parameters(Opti()).referenceSolution));
            analytical = std::make_shared<AnalyticalTypecast<3,BT,T,BT,T>>(*help);
          } break;
        }

        simulated = std::make_shared<SuperLatticePhysVelocity3D<T,descriptor>>(
          this->lattice(), this->converter());
      } break;
    }

    RelativeDifferenceObjective3D<T,descriptor>(
      this->lattice(),
      simulated,
      analytical,
      this->geometry().getMaterialIndicator(this->parameters(Opti()).objectiveMaterial)
    ).operator()(output, tmp);

    return output[0];
  }
};



// =================== Variant B: optimizable solver ==========================

/// Parameters for optimizable TestFlow: has to inherit interface from OptiSimulation
/** cuboidWise=true: a cuboid decomposition of the domain will be created and the
 * control values are applied as force factors for each cuboid.
 * cuboidWise=false: 3 control values are expected and applied globally as
 * force factors (component-wise).
 */
template<typename T>
struct TfDirectOpti : public parameters::DirectOptiSimulation<T> {

  // cf. above class description
  bool                      cuboidWise        {false};
  /// scales the force field (comes from control variable)
  std::vector<T>            forceFactor       {};
  /// material number, where the objective is computed
  int                       objectiveMaterial {1};
  /// define whether VELOCITY or DISSIPATION error is used for objective computation
  TestFlowOptiMode          optiMode          {VELOCITY};
  /// define whether DISCRETE or ANALYTICAL solution is used to compute the objective
  OptiReferenceSolution     optiReferenceMode {DISCRETE};
  bool                      computeObjective  {true};

  /// here, we can store a discrete reference solution for objective computation
  using BT = BaseType<T>;
  std::shared_ptr<SuperF3D<BT>> referenceSolution;

  /// define how the control variable enters the simulation
  void applyControl(const std::vector<T>& control) override {
    forceFactor = control;
  }
};

/// XML interface for TestFlowOptiSimulation parameters
template<typename T, typename TAG>
struct parameters::Reader<TfDirectOpti<T>, TAG>
 : public parameters::ReaderBase<TfDirectOpti<T>>
{
  using parameters::ReaderBase<TfDirectOpti<T>>::ReaderBase;

  void read(XMLreader const& xml) {
    xml.readOrWarn<bool>("Optimization", "CuboidWiseControl", "",
      this->params->cuboidWise, true, false, true);
    xml.readOrWarn<int>("Optimization", "ObjectiveMaterial", "",
      this->params->objectiveMaterial, true, false, false);

    std::string helper("velocity");
    xml.readOrWarn<std::string>("Optimization", "TestFlowOptiMode", "",
      helper, true, false, true);
    if (helper=="velocity") {
      this->params->optiMode = VELOCITY;
    } else if (helper=="dissipation") {
      this->params->optiMode = DISSIPATION;
    } else {
      throw std::invalid_argument("invalid argument in xml [Optimization][TestFlowOptiMode]");
    }

    helper = "analytical";
    xml.readOrWarn<std::string>("Optimization", "OptiReferenceMode", "",
      helper, true, false, true);
    if (helper=="discrete") {
      this->params->optiReferenceMode = DISCRETE;
    } else if (helper=="analytical") {
      this->params->optiReferenceMode = ANALYTICAL;
    } else {
      throw std::invalid_argument("invalid argument in xml [Optimization][OptiReferenceMode]");
    }
  }
};

/** Optimization-related results of test-flow simulation
 * Additionally save velocity or dissipation functor to allow usage of
 * discrete objective
 */
template<typename T>
struct TfDirectOptiResults : public parameters::DirectOptiResults<T> {

  std::shared_ptr<SuperF3D<T>> solution;
};


template<typename T>
using Params_TfDirectOpti = meta::map<
  Opti,             TfDirectOpti<T>,
  Output,           parameters::OutputGeneral<T>,
  OutputOpti,       parameters::OptiOutput<T,SolverMode::Primal>,
  Results,          TfDirectOptiResults<T>,
  Simulation,       TfSimulationParams<T,Lattices>,
  VisualizationVTK, parameters::OutputPlot<T>
>;


template<typename T>
class TestFlowSolverDirectOpti : public TestFlowOptiBase<T,Params_TfDirectOpti<T>,Lattices> {
private:
  mutable OstreamManager            clout {std::cout, "TestFlowSolverDirectOpti"};
  using typename TestFlowSolverDirectOpti::TestFlowOptiBase::descriptor;

public:
  TestFlowSolverDirectOpti(
    utilities::TypeIndexedSharedPtrTuple<Params_TfDirectOpti<T>> params)
   : TestFlowSolverDirectOpti::TestFlowOptiBase(params)
  { }

protected:
  void prepareLattices() override
  {
    // call prepareLattices from base class
    TestFlowSolverDirectOpti::TestFlowBase::prepareLattices();

    // scale force field according to control
    std::shared_ptr<AnalyticalF<3,T,T>> controlF;
    if (this->parameters(Opti()).cuboidWise) {
      controlF = std::make_shared<AnalyticalCuboidwiseConst<3,T,T>>(
        this->geometry(), this->parameters(Opti()).forceFactor);
    } else {
      OLB_ASSERT((this->parameters(Opti()).forceFactor.size()==3),
        "A control vector of size 3 is required");
      controlF = std::make_shared<AnalyticalConst3D<T,T>>(
        this->parameters(Opti()).forceFactor);
    }
    OLB_ASSERT((controlF->getTargetDim()==3),
      "Dimension of force scaling functor must be 3");
    this->_force = this->_force * controlF;

    this->lattice().template defineField<descriptors::FORCE>
    (this->geometry().getMaterialIndicator({1, 2}), *(this->_force));
  }

  void computeResults() override
  {
    if (this->parameters(Opti()).computeObjective) {
      this->parameters(Results()).objective = this->computeObjective();
    }

    switch (this->parameters(Opti()).optiMode)
    {
      case DISSIPATION:
      {
        this->parameters(Results()).solution
        = std::make_shared<SuperLatticePhysDissipation3D<T,descriptor>>(
          this->lattice(), this->converter());
      } break;

      case VELOCITY:
      {
        this->parameters(Results()).solution
        = std::make_shared<SuperLatticePhysVelocity3D<T,descriptor>>(
          this->lattice(), this->converter());
      } break;
    }
  }
};


// ================ Variant C: Adjoint optimizable solver =====================

/// Parameters for optimizable TestFlow: has to inherit interface from OptiSimulation
template<typename T, opti::SolverMode MODE>
struct TestFlowAdjointOpti
 : public parameters::DistributedOptiSimulation<T,DualLattices,MODE> {

  int                       objectiveMaterial {1};
  TestFlowOptiMode          optiMode          {VELOCITY};
  OptiReferenceSolution     optiReferenceMode {DISCRETE};
  std::shared_ptr<SuperIndicatorF3D<T>> designDomain {};
};

/// XML interface for TestFlowOptiSimulation parameters
template<typename T, opti::SolverMode MODE, typename TAG>
struct parameters::Reader<TestFlowAdjointOpti<T,MODE>, TAG>
 : public parameters::ReaderBase<TestFlowAdjointOpti<T,MODE>>
{
  using parameters::ReaderBase<TestFlowAdjointOpti<T,MODE>>::ReaderBase;

  void read(XMLreader const& xml) {
    Reader<DistributedOptiSimulation<T,DualLattices,MODE>, TAG>(this->params).read(xml);

    xml.readOrWarn<int>("Optimization", "ObjectiveMaterial", "",
      this->params->objectiveMaterial, true, false, false);

    std::string helper ("analytical");
    xml.readOrWarn<std::string>("Optimization", "OptiReferenceMode", "",
      helper, true, false, true);
    if (helper=="discrete") {
      this->params->optiReferenceMode = DISCRETE;
    } else if (helper=="analytical") {
      this->params->optiReferenceMode = ANALYTICAL;
    } else {
      throw std::invalid_argument("invalid argument in xml [Optimization][OptiReferenceMode]");
    }
  }
};


template<typename T, SolverMode MODE>
using Params_TfOptiAdjoint = meta::map<
  Opti,             TestFlowAdjointOpti<T,MODE>,
  OutputOpti,       parameters::OptiOutput<T,MODE>,
  Output,           parameters::OutputGeneral<T>,
  Results,          parameters::DistributedOptiSimulationResults<T,DualLattices,MODE>,
  Simulation,       TfSimulationParams<T,DualLattices>,
  VisualizationVTK, parameters::OutputPlot<T>
>;

template<typename T, SolverMode MODE>
class TestFlowSolverOptiAdjoint : public TestFlowOptiBase<
  T,Params_TfOptiAdjoint<T,MODE>,DualLattices,MODE>
{
private:
  mutable OstreamManager            clout {std::cout, "TestFlowSolverOptiAdjoint"};

public:
  using lattices = DualLattices;
  using typename TestFlowSolverOptiAdjoint::TestFlowOptiBase::descriptor;

  TestFlowSolverOptiAdjoint(
    utilities::TypeIndexedSharedPtrTuple<Params_TfOptiAdjoint<T,MODE>> params)
   : TestFlowSolverOptiAdjoint::TestFlowOptiBase(params)
  { }

  void prepareGeometry() override
  {
    TestFlowSolverOptiAdjoint::TestFlowBase::prepareGeometry();

    this->parameters(Opti()).designDomain =
      std::move(this->geometry().getMaterialIndicator(this->parameters(Opti()).controlMaterial));
  }

  void setInitialValues() override
  {
    TestFlowSolverOptiAdjoint::TestFlowBase::setInitialValues();

    if constexpr ( MODE == SolverMode::Dual ) {
      this->loadPrimalPopulations();
    }
  }

  void prepareVTK() const override { }

  void writeVTK(std::size_t iT) const override
  {
    if constexpr (MODE == SolverMode::Reference) {
      // ARTIFICIAL DATA: write solution (unsteady, until convergence in time is obtained)
      SuperVTMwriter3D<T> writer("Reference_Solution");

      if (iT == 0) {
        // Writes the geometry, cuboid no. and rank no. as vti file for visualization
        SuperLatticeGeometry3D<T,descriptor> geometry(this->lattice(), this->geometry());
        SuperLatticeCuboid3D<T,descriptor> cuboid(this->lattice());
        SuperLatticeRank3D<T,descriptor> rank(this->lattice());

        writer.write(cuboid);
        writer.write(geometry);
        writer.write(rank);
        writer.createMasterFile();
      }

      this->writeFunctorsToVTK(writer, iT);
    }

    else if constexpr (MODE == SolverMode::Primal) {
      // OPTIMISATION: Write stationary flow field
      if (this->_finishedTimeLoop) {
        SuperVTMwriter3D<T> writer("Flow_Simulations");

        if (this->parameters(OutputOpti()).counterOptiStep == 0) {
          writer.createMasterFile();
        }

        this->writeFunctorsToVTK(writer, this->parameters(OutputOpti()).counterOptiStep);
      }
    }

    else {    // Dual (adjoint) mode
      if (this->_finishedTimeLoop) {
        SuperVTMwriter3D<T> writer("Dual_Simulations");

        if (this->parameters(OutputOpti()).counterOptiStep==0) {
          writer.createMasterFile();
        }

        this->writeFunctorsToVTK(writer, this->parameters(OutputOpti()).counterOptiStep);
      }
    }
  }

  /// Store some pointers and fields s.t. they can be passed to the adjoint Solver
  void computeResults() override
  {
    // yield pointer to lattice
    this->parameters(Results()).lattice = std::get<0>(this->_sLattices);

    if constexpr ((MODE == SolverMode::Reference) || (MODE == SolverMode::Primal))
    {
      this->parameters(Results()).geometry = this->_sGeometry;
    }

    if constexpr (MODE == SolverMode::Primal)
    {
      // compute objective
      this->parameters(Results()).objective = this->computeObjective();
      this->parameters(Results()).objectiveComputed = true;

      // set population
      this->parameters(Results()).fpop = std::make_shared<SuperLatticeFpop3D<T,descriptor>>(
        this->lattice());

      // set objective derivative
      this->parameters(Results()).djdf = getObjectiveDerivative();

      // set control derivative
      std::shared_ptr<AnalyticalF<3,T,T>> zero = std::make_shared<AnalyticalConst<3,T,T>> (Vector<T,3>(0));
      this->parameters(Results()).djdalpha = std::make_shared<SuperLatticeFfromAnalyticalF3D<T,descriptor>> (
        zero, this->lattice());
    }
  }

  std::shared_ptr<SuperF3D<T,T>> getObjectiveDerivative() const
  {
    std::shared_ptr<SuperF3D<T,T>>      simulated;
    std::shared_ptr<SuperF3D<T,T>>      simulatedDerivative;
    std::shared_ptr<AnalyticalF<3,T,T>> analytical;

    switch (this->parameters(Opti()).optiMode)
    {
      case DISSIPATION:
      {
        simulated = std::make_shared<SuperLatticePhysDissipation3D<T,descriptor>>(
          this->lattice(), this->converter());
        simulatedDerivative = std::make_shared<SuperLatticeDphysDissipationDf3D<T,descriptor>> (
          this->lattice(), this->converter());

        switch (this->parameters(Opti()).optiReferenceMode)
        {
          case ANALYTICAL:
          {
            analytical = std::make_shared<DissipationTestFlow3D<T,T,descriptor>>(
              this->converter());
          } break;

          case DISCRETE:
          {
            // user has to guarantee that reference solution is computed
            // cf. node Optimization -> ReferenceSolution in xml file
            analytical = std::make_shared<AnalyticalFfromSuperF3D<T>>(
              *(this->parameters(Opti()).referenceSolution));
          } break;
        }
      } break;

      case VELOCITY:
      {
        simulated = std::make_shared<SuperLatticePhysVelocity3D<T,descriptor>>(
          this->lattice(), this->converter());
        simulatedDerivative = std::make_shared<SuperLatticeDphysVelocityDf3D<T,descriptor>> (
          this->lattice(), this->converter());

        switch (this->parameters(Opti()).optiReferenceMode)
        {
          case ANALYTICAL:
          {
            analytical = std::make_shared<VelocityTestFlow3D<T,T,descriptor>>(
              this->converter());
          } break;

          case DISCRETE:
          {
            // user has to guarantee that reference solution is computed
            // cf. node Optimization -> ReferenceSolution in xml file
            analytical = std::make_shared<AnalyticalFfromSuperF3D<T>>(
              *(this->parameters(Opti()).referenceSolution));
          } break;
        }
      } break;
    }

    return std::make_shared<DrelativeDifferenceObjectiveDf3D<T,descriptor>>(
      this->lattice(),
      simulated,
      simulatedDerivative,
      analytical,
      this->geometry().getMaterialIndicator(this->parameters(Opti()).objectiveMaterial));
  }
};

#endif
