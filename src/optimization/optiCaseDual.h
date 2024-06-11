/*  This file is part of the OpenLB library
*
*  Copyright (C) 2012-2021 Mathias J. Krause, Benjamin Förster, Julius Jeßberger
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
* An OptiCase using Adjoint-based Differentiation
*/


#ifndef OPTI_CASE_DUAL_H
#define OPTI_CASE_DUAL_H

#include "io/xmlReader.h"

#include "adjointLbSolver.h"
#include "controller.h"
#include "optiCase.h"
#include "projection.h"
#include "serialization.h"

namespace olb {

namespace opti {

enum ControlType {ForceControl, PorosityControl};

/** This class implements the evaluation of the goal functional
 * and its derivatives by using adjoint LBM.
 * The adjoint equations are problem-specific and have been computed for force
 * and porosity optimization so far.
 *
 * Requirements: S is the arithmetic data type
 * SOLVER implements the (primal/ dual) simulation, inherits from AdjointLbSolverBase
 * An xml file is expected to provide additional information on e.g. simulation
 * and optimization parameters
 */
template<
  typename S,
  template<typename,SolverMode> typename SOLVER,
  typename C = std::vector<S>>
class OptiCaseDual : public OptiCase<S,C> {

private:
  mutable OstreamManager                           clout {std::cout, "OptiCaseDual"};
  using descriptor = typename SOLVER<S,SolverMode::Reference>::AdjointLbSolverBase::DESCRIPTOR;
  static constexpr int dim = descriptor::d;

public:
  bool                                             _verbose {true};
  /// upper limit for the number of control variables (#voxels * field-dimension)
  std::size_t                                      _dimCtrl;

  /// Either force or porosity field
  ControlType                                      _controlType;
  StartValueType                                   _startValueType {Control};
  std::string                                      _projectionName;

  /// Spatial dimension of controlled field
  int                                              _fieldDim;
  /// Material number of design domain
  int                                              _controlMaterial;
  /// Regulatory term in objective functional (so far unused)
  S                                                _regAlpha {0};

  /// Manages the array of control variables
  Controller<S>*                                   _controller {nullptr};
  UnitConverter<S,descriptor>*                     _converter {nullptr};
  std::shared_ptr<SOLVER<S,SolverMode::Primal>>    _primalSolver;
  std::shared_ptr<SOLVER<S,SolverMode::Dual>>      _dualSolver;
  std::shared_ptr<SOLVER<S,SolverMode::Reference>> _referenceSolver;

  bool                                             _computeReference {false};
  std::shared_ptr<SuperGeometry<S,dim>>            _referenceGeometry;
  std::shared_ptr<SuperLattice<S,descriptor>>      _referenceLattice;

  std::shared_ptr<GeometrySerializer<S,dim>>       _serializer;
  std::shared_ptr<projection::Base<S>>             _projection;
  /// maps the control to a lattice functor
  std::shared_ptr<SuperLatticeF<S,descriptor>>     _projectedControl;
  /// derivative of _projectionControl
  std::shared_ptr<SuperLatticeF<S,descriptor>>     _dProjectionDcontrol;

  OptiCaseDual(XMLreader const& xml) {
    readFromXML(xml);
    initialize(xml);
    initializeFields();
  }

  ~OptiCaseDual() {
    free();
  }

  void free() {
    delete _converter;
    _converter = nullptr;
    delete _controller;
    _controller = nullptr;
  }

  /// Solve primal problem and evaluate objective
  S evaluateObjective(const C& control, unsigned optiStep=0) override;

  /// Compute derivatives via adjoint problem
  void computeDerivatives(
    const C& control, C& derivatives, unsigned optiStep=0) override;

  // get the control values as they were computed in the reference solution
  // only ForceControl & without projection is implemented so far
  C getReferenceControl() const;

private:
  void readFromXML(XMLreader const& xml);

  void initialize(XMLreader const& xml);

  void initializeFields();

  void derivativesFromDualSolution(C& derivatives);

};


}

}

#endif
