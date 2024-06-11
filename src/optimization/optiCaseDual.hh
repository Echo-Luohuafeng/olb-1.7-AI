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


#ifndef OPTI_CASE_DUAL_HH
#define OPTI_CASE_DUAL_HH

#include "optiCaseDual.h"
#include "serialization.h"

namespace olb {

namespace opti {

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::readFromXML(XMLreader const& xml)
{
  std::string type ("");
  xml.readOrWarn<std::string>("Optimization", "ControlType", "", type);
  _controlType = ((type == "Force") || (type == "force")) ? ForceControl : PorosityControl;
  xml.readOrWarn<std::string>("Optimization", "Projection", "", _projectionName);
  xml.readOrWarn<int>("Optimization", "ControlMaterial", "", _controlMaterial);
  xml.readOrWarn<S>("Optimization", "RegAlpha", "", _regAlpha);
  type = "";
  xml.readOrWarn<std::string>("Optimization", "StartValueType", "", type);
  if (type == "Porosity") {
    _startValueType = Porosity;
  } else if (type == "Permeability") {
    _startValueType = Permeability;
  } else if (type == "Control") {
    _startValueType = Control;
  } else {
    _startValueType = ProjectedControl;
  }
   xml.readOrWarn<bool>("Optimization", "ReferenceSolution", "", _computeReference);
}

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::initialize(XMLreader const& xml)
{
  _converter = createUnitConverter<S,descriptor>(xml);

  _referenceSolver = createLbSolver <SOLVER<S,SolverMode::Reference>> (xml);
  if (_computeReference) {
    _referenceSolver->solve();
  } else {
    // only build up geometry and lattice
    _referenceSolver->buildAndReturn();
  }

  _referenceGeometry = _referenceSolver->parameters(names::Results()).geometry;
  _referenceLattice = _referenceSolver->parameters(names::Results()).lattice;

  if (_controlType == ForceControl) {
    _fieldDim = dim;
  } else {
    _fieldDim = 1;
  }

  _serializer = std::make_shared<SimpleGeometrySerializer<S,dim>>(*_referenceGeometry);
  _dimCtrl = _serializer->getNoCells() * _fieldDim;

  _controller = new Controller<S>(_dimCtrl);

  _primalSolver = createLbSolver <SOLVER<S,SolverMode::Primal>> (xml);
  _dualSolver = createLbSolver <SOLVER<S,SolverMode::Dual>> (xml);
  this->_postEvaluation = std::bind(&SOLVER<S,SolverMode::Primal>::postProcess, _primalSolver);

  if (_computeReference) {
    _primalSolver->parameters(names::Opti()).referenceSolution
     = std::make_shared<SuperLatticePhysVelocity3D<S,descriptor>>(*_referenceLattice, *_converter);
    _primalSolver->parameters(names::Opti()).referencePorosity
     = std::make_shared<SuperLatticePorosity3D<S,descriptor>>(*_referenceLattice);
  }
}

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::initializeFields()
{
  _projection = projection::construct<S,descriptor>(*_converter, _projectionName);

  _projectedControl = std::make_shared<SuperLatticeSerialDataF<S,descriptor>>(
    *_referenceLattice,
    *_controller,
    _fieldDim,
    _serializer,
    [this](S x) { return _projection->project(x); });
  _dProjectionDcontrol = std::make_shared<SuperLatticeSerialDataF<S,descriptor>>(
    *_referenceLattice,
    *_controller,
    _fieldDim,
    _serializer,
    [this](S x) { return _projection->derivative(x); });
  _primalSolver->parameters(names::Opti()).controlledField = _projectedControl;
  _dualSolver->parameters(names::Opti()).controlledField = _projectedControl;
}


template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
S OptiCaseDual<S,SOLVER,C>::evaluateObjective(
  const C& control, unsigned optiStep)
{
  //project control to enforce volume constraint
#ifdef mathiasProjection

  const int nx = _dualSolver->getGeometry()->getCuboidGeometry().getMotherCuboid().getNx();
  const int ny = _dualSolver->getGeometry()->getCuboidGeometry().getMotherCuboid().getNy();
  const int nz = _dualSolver->getGeometry()->getCuboidGeometry().getMotherCuboid().getNz();

  S pi = 4.0 * util::atan(1.0);

  S totalD = T();
  S wantedTotalD = T();

  S upperBound, lowerBound, volumeRatio;
  xml.readOrWarn<T>("Optimization", "UpperBound", "", upperBound);
  xml.readOrWarn<T>("Optimization", "LowerBound", "", lowerBound);
  xml.readOrWarn<T>("Optimization", "VolumeRatio", "", volumeRatio);

  for (int iX=1; iX<nx-1; iX++) {
    for (int iY=1; iY<ny-1; iY++) {
      for (int iZ=1; iZ<nz-1; iZ++) {
        // TODO: this seems to be fixed to material number = 6
        if (_dualSolver->getGeometry()->getMaterial(iX,iY,iZ)==6) {
          S ctrl = control[(iX-1)*(ny-1)*(nz-1)+(iY-1)*(nz-1)+(iZ-1)];
          totalD += 1./(upperBound-lowerBound)*(ctrl-1.);
          wantedTotalD += 1.0;
        }
      }
    }
  }

  for (int iX=1; iX<nx-1; iX++) {
    for (int iY=1; iY<ny-1; iY++) {
      for (int iZ=1; iZ<nz-1; iZ++) {
        // TODO: this seems to be fixed to material number = 6
        if (_dualSolver->getGeometry()->getMaterial(iX,iY,iZ)==6) {
          control[(iX-1)*(ny-1)*(nz-1)+(iY-1)*(nz-1)+(iZ-1)]
           = lowerBound
             + (control[(iX-1)*(ny-1)*(nz-1)+(iY-1)*(nz-1)+(iZ-1)] - lowerBound)
               * wantedTotalD * (1.-volumeRatio) / totalD;
        }
      }
    }
  }
#endif

  _controller->setControl(control, _dimCtrl);
  _primalSolver->parameters(names::OutputOpti()).counterOptiStep = optiStep;
  _primalSolver->solve();

  return _primalSolver->parameters(names::Results()).objective;
}


template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::computeDerivatives(
  const C& control, C& derivatives, unsigned optiStep)
{
  _controller->setControl(control, _dimCtrl);
  _dualSolver->parameters(names::OutputOpti()).counterOptiStep = optiStep;

  const auto& primalResults = _primalSolver->parameters(names::Results());
  auto& dualParams = _dualSolver->parameters(names::Opti());

  dualParams.fpop = primalResults.fpop;
  dualParams.dObjectiveDf = primalResults.djdf;
  dualParams.dObjectiveDcontrol = primalResults.djdalpha;

  _dualSolver->solve();

  derivativesFromDualSolution(derivatives);
}

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::derivativesFromDualSolution(
  C& derivatives)
{
  const auto& primalGeometry = _primalSolver->parameters(names::Results()).geometry;
  const S omega = _converter->getLatticeRelaxationFrequency();
  const int nC = primalGeometry->getCuboidGeometry().getNc();

  for (int iC = 0; iC < nC; iC++) {

    const Vector<int,3> extend = primalGeometry->getCuboidGeometry().get(iC).getExtent();

    for (int iX = 0; iX < extend[0]; iX++) {
      for (int iY = 0; iY < extend[1]; iY++) {
        for (int iZ = 0; iZ < extend[2]; iZ++) {
          const LatticeR<dim+1> latticeR(iC, iX, iY, iZ);

          if (getMaterialGlobally(*_referenceGeometry, latticeR) == _controlMaterial) {
            C derivativesHelp(_fieldDim, 0);

            if (primalGeometry->getLoadBalancer().rank(iC)  == singleton::mpi().getRank()) {
              const S rho_f = _primalSolver->parameters(names::Results()).lattice->get(latticeR).computeRho();
              S dProjectionDcontrol[_fieldDim];
              (*_dProjectionDcontrol)(dProjectionDcontrol, latticeR.data());

              if ( _controlType == ForceControl ) {
                S dObjectiveDcontrol[_fieldDim];
                (*_primalSolver->parameters(names::Results()).djdalpha)(dObjectiveDcontrol, latticeR.data());
                for (int iDim=0; iDim<_fieldDim; iDim++) {
                  for (int jPop=0; jPop < descriptor::q; ++jPop) {
                    S phi_j = _dualSolver->parameters(names::Results()).lattice->get(latticeR)[jPop];
                    derivativesHelp[iDim] -= rho_f * descriptors::t<S,descriptor>(jPop)
                       * descriptors::invCs2<S,descriptor>() * descriptors::c<descriptor>(jPop,iDim) * phi_j;
                  }
                  // dObjectiveDcontrol is zero if objective does not depend on control
                  // --> dObjectiveDcontrol is 0 and therefore dProjectionDcontrol may be
                  // irrelevant for force and porosity optimisation
                  const int index = _serializer->getSerializedComponentIndex(latticeR, iDim, _fieldDim);
                  derivativesHelp[iDim] += _regAlpha * _controller->getControl(index)
                       + dObjectiveDcontrol[iDim] * dProjectionDcontrol[iDim];
                }
              }
              else if ( _controlType == PorosityControl ) {
                S u_f[3];
                _primalSolver->parameters(names::Results()).lattice->get(latticeR).computeU(u_f);
                const S d = _dualSolver->parameters(names::Results()).lattice->get(latticeR).template getField<descriptors::POROSITY>();

                for (int jPop = 0; jPop < descriptor::q; ++jPop) {
                  S phi_j = _dualSolver->parameters(names::Results()).lattice->get(latticeR)[jPop];
                  S feq_j = equilibrium<descriptor>::secondOrder(jPop, rho_f, u_f) + descriptors::t<S,descriptor>(jPop);
                  for (int iDim = 0; iDim < descriptor::d; iDim++) {
                    derivativesHelp[0] +=  phi_j*feq_j*( descriptors::c<descriptor>(jPop,iDim) - d*u_f[iDim] )*u_f[iDim]*dProjectionDcontrol[0];
                  }
                }
                derivativesHelp[0] *= -omega*descriptors::invCs2<S,descriptor>();
              }
              // correct processing of regularizing term has to be checked in both cases!
            }

#ifdef PARALLEL_MODE_MPI
            singleton::mpi().bCast(&derivativesHelp[0], _fieldDim, primalGeometry->getLoadBalancer().rank(iC));
#endif
            for (int iDim=0; iDim<_fieldDim; ++iDim) {
              const int index = _serializer->getSerializedComponentIndex(latticeR, iDim, _fieldDim);
              derivatives[index] = derivativesHelp[iDim];
            }
          }
        }
      }
    }
  }
}

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
C OptiCaseDual<S,SOLVER,C>::getReferenceControl() const
{
  assert((this->_computeReference) && "Reference solution has not yet been computed\n");
  C result = util::ContainerCreator<C>::create(_dimCtrl);

  if (_controlType == ForceControl) {
      LatticeR<dim+1> latticeR;
      for (int iC=0; iC<_referenceGeometry->getCuboidGeometry().getNc(); iC++) {
      latticeR[0] = iC;
      int nX = _referenceGeometry->getCuboidGeometry().get(iC).getNx();
      int nY = _referenceGeometry->getCuboidGeometry().get(iC).getNy();
      int nZ = _referenceGeometry->getCuboidGeometry().get(iC).getNz();
      for (int iX=0; iX<nX; iX++) {
        latticeR[1] = iX;
        for (int iY=0; iY<nY; iY++) {
          latticeR[2] = iY;
          for (int iZ=0; iZ<nZ; iZ++) {
            latticeR[3] = iZ;

            if (getMaterialGlobally(*_referenceGeometry, latticeR) == _controlMaterial) {
              C force_help(_fieldDim, 0);
              if (_referenceGeometry->getLoadBalancer().rank(iC)  == singleton::mpi().getRank()) {
                for (int iDim=0; iDim<_fieldDim; iDim++) {
                  const auto cell = _referenceLattice->get(latticeR);
                  force_help[iDim]
                   = _projection->inverse(cell.template getFieldComponent<descriptors::FORCE>(iDim));
                }
              }
#ifdef PARALLEL_MODE_MPI
              singleton::mpi().bCast(&force_help[0], _fieldDim, _referenceGeometry->getLoadBalancer().rank(iC));
#endif
              for (int iDim=0; iDim<_fieldDim; iDim++) {
                const int index = _serializer->getSerializedComponentIndex(latticeR, iDim, _fieldDim);
                result[index] = force_help[iDim];
              }
            }
          }
        }
      }
    }
  } else {
    clout << "getReferenceControl() is only implemented for forced optimization" << std::endl;
    exit(1);
  }
  return result;
}

} // namespace opti

} // namespace olb

#endif
