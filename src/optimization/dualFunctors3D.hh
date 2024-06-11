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

#ifndef DUAL_FUNCTORS_3D_HH
#define DUAL_FUNCTORS_3D_HH


#include "optimization/dualFunctors3D.h"


namespace olb {

namespace opti {

template <typename T, typename DESCRIPTOR>
BlockLatticeDphysDissipationDf3D<T,DESCRIPTOR>::BlockLatticeDphysDissipationDf3D(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,DESCRIPTOR::q),
    _overlap(overlap),
    _converter(converter)
{
  this->getName() = "dPhysDissipationDf";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeDphysDissipationDf3D<T,DESCRIPTOR>::operator()(T dDdf[], const int latticeR[])
{
  T omega = _converter.getLatticeRelaxationFrequency();
  T dt = _converter.getConversionFactorTime();

  T rho;
  T u[3];
  T pi[6];

  /*const int localLatticeR[3] = {
    latticeR[0] + _overlap,
    latticeR[1] + _overlap,
    latticeR[2] + _overlap
  };
  this->_blockLattice.get(localLatticeR).computeAllMomenta(rho, u, pi);*/
  this->_blockLattice.get(latticeR).computeAllMomenta(rho, u, pi);

  for (int i=0; i<DESCRIPTOR::q; i++) { //output default 0
    dDdf[i] = T();
  }

  T pi2[DESCRIPTOR::d][DESCRIPTOR::d];

  pi2[0][0]=pi[0];
  pi2[0][1]=pi[1];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    pi2[0][2]=pi[2];
  }
  else {
    pi2[1][1]=pi[2];
  }
  pi2[1][0]=pi[1];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    pi2[1][1]=pi[3];
    pi2[1][2]=pi[4];
    pi2[2][0]=pi[2];
    pi2[2][1]=pi[4];
    pi2[2][2]=pi[5];
  }

  for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
    for (int iDim=0; iDim < DESCRIPTOR::d; ++iDim) {
      for (int jDim=0; jDim < DESCRIPTOR::d; ++jDim) {
        T dPidf2ndTerm = T();
        for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          dPidf2ndTerm += descriptors::c<DESCRIPTOR>(iPop)[iDim] * descriptors::c<DESCRIPTOR>(iPop)[jDim]
                          * dualLbHelpers<T,DESCRIPTOR>::equilibrium(iPop, jPop, rho, u);
        }
        T dpidf = descriptors::c<DESCRIPTOR>(jPop,iDim)*descriptors::c<DESCRIPTOR>(jPop,jDim) - dPidf2ndTerm;
        // rho2 instead of rho vs simon!
        dDdf[jPop] += 2.*(pi2[iDim][jDim]*(dpidf - pi2[iDim][jDim]/2./rho*2.));
      }
    }
    dDdf[jPop] *= util::pow(omega*descriptors::invCs2<T,DESCRIPTOR>()/rho,2)/2.*_converter.getPhysViscosity()/_converter.getLatticeViscosity()/dt/dt;
  }
  return true;
}

template <typename T, typename DESCRIPTOR>
SuperLatticeDphysDissipationDf3D<T,DESCRIPTOR>::SuperLatticeDphysDissipationDf3D(
  SuperLattice<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice, converter, DESCRIPTOR::q*DESCRIPTOR::d)
{
  this->getName() = "dPhysDissipationDf";
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticeDphysDissipationDf3D<T,DESCRIPTOR>(
        sLattice.getBlock(iC),
        sLattice.getOverlap(),
        converter)
    );
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeDphysDissipationDf3D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  auto& load = this->_sLattice.getLoadBalancer();
  if (load.isLocal(input[0])) {
    return this->getBlockF(load.loc(input[0]))(output, &input[1]);
  }
  else {
    return false;
  }
}


template <typename T, typename DESCRIPTOR>
BlockLatticeDphysVelocityDf3D<T,DESCRIPTOR>::BlockLatticeDphysVelocityDf3D(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter,
  int nDim, int extractDim)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,DESCRIPTOR::q*DESCRIPTOR::d),
    _overlap(overlap),
    _converter(converter),
    _nDim(nDim), _extractDim(extractDim)
{
  this->getName() = "dPhysVelocityDf";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeDphysVelocityDf3D<T,DESCRIPTOR>::operator()(T dVelocityDf[], const int latticeR[])
{
  for (int i=0; i<DESCRIPTOR::d*DESCRIPTOR::q; i++) {
    dVelocityDf[i] = T();
  }

  T rho;
  T u[3];
  /*const int localLatticeR[3] = {
    latticeR[0] + _overlap,
    latticeR[1] + _overlap,
    latticeR[2] + _overlap
  };

  this->_blockLattice.get(localLatticeR).computeRhoU(rho, u);*/
  this->_blockLattice.get(latticeR).computeRhoU(rho, u);

  for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
    for (int iDim=0; iDim < DESCRIPTOR::d; ++iDim) {
      if (iDim != _extractDim && _extractDim != -1) {
        dVelocityDf[jPop*DESCRIPTOR::d + iDim] = T(0);
      }
      else {
        dVelocityDf[jPop*DESCRIPTOR::d + iDim] = (descriptors::c<DESCRIPTOR>(jPop,iDim)-u[iDim])/rho;
      }
    }
  }
  return true;
}

template <typename T, typename DESCRIPTOR>
SuperLatticeDphysVelocityDf3D<T,DESCRIPTOR>::SuperLatticeDphysVelocityDf3D(
  SuperLattice<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice, converter, DESCRIPTOR::q)
{
  this->getName() = "dPhysVelocityDf";
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticeDphysVelocityDf3D<T,DESCRIPTOR>(
        sLattice.getBlock(iC),
        sLattice.getOverlap(),
        converter,
        -1, -1)
    );
  }
}

template <typename T, typename DESCRIPTOR>
SuperLatticeDphysVelocityDf3D<T,DESCRIPTOR>::SuperLatticeDphysVelocityDf3D(
  SuperLattice<T,DESCRIPTOR>&      sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperExtractComponentF3D<T,T>&  sExtract)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice, converter, DESCRIPTOR::q)
{
  this->getName() = "dPhysVelocityDf";
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticeDphysVelocityDf3D<T,DESCRIPTOR>(
        sLattice.getBlock(iC),
        sLattice.getOverlap(),
        converter,
        1, sExtract.getExtractDim())
    );
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeDphysVelocityDf3D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  auto& load = this->_sLattice.getLoadBalancer();
  if (load.isLocal(input[0])) {
    return this->getBlockF(load.loc(input[0]))(output, &input[1]);
  }
  else {
    return false;
  }
}


template <typename T, typename DESCRIPTOR>
DifferenceObjective3D<T,DESCRIPTOR>::DifferenceObjective3D(
  SuperLattice<T, DESCRIPTOR>&     sLattice,
  FunctorPtr<SuperF3D<T,T>>&&        f,
  FunctorPtr<AnalyticalF3D<T,T>>&&   wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperIdentity3D<T,T>([&]()
{
  using namespace functor_dsl;

  auto wantedLatticeF = restrictF(wantedF.toShared(), sLattice);
  auto differenceF = norm<2>(f.toShared() - wantedLatticeF, indicatorF.toShared());

  return util::pow(differenceF, 2) / T(2.0);
}())
{
  this->getName() = "differenceObjective";
}

template<typename T, typename DESCRIPTOR>
DifferenceObjective3D<T,DESCRIPTOR>::DifferenceObjective3D(
  FunctorPtr<SuperLatticeF3D<T,DESCRIPTOR>>&& f,
  FunctorPtr<AnalyticalF3D<T,T>>&&            wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&&          indicatorF)
  : DifferenceObjective3D(
      f->getSuperLattice(),
      std::static_pointer_cast<SuperF3D<T,T>>(f.toShared()),
      std::forward<decltype(wantedF)>(wantedF),
      std::forward<decltype(indicatorF)>(indicatorF))
{ }

template<typename T, typename DESCRIPTOR>
DifferenceObjective3D<T,DESCRIPTOR>::DifferenceObjective3D(
  SuperLattice<T, DESCRIPTOR>&   sLattice,
  FunctorPtr<SuperF3D<T,T>>&&      f,
  FunctorPtr<AnalyticalF3D<T,T>>&& wantedF,
  SuperGeometry<T,3>& geometry,
  int                 material)
  : DifferenceObjective3D(
      sLattice,
      std::forward<decltype(f)>(f),
      std::forward<decltype(wantedF)>(wantedF),
      geometry.getMaterialIndicator(material))
{ }


template <typename T, typename DESCRIPTOR>
BlockDdifferenceObjectiveDf3D<T,DESCRIPTOR>::BlockDdifferenceObjectiveDf3D(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  BlockF3D<T>&          f,
  BlockF3D<T>&          dFdF,
  BlockF3D<T>&          wantedF,
  BlockIndicatorF3D<T>& indicatorF,
  T weight)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, f.getTargetDim()),
    _f(f),
    _dFdF(dFdF),
    _wantedF(wantedF),
    _indicatorF(indicatorF),
    _weight(weight)
{
  this->getName() = "BlockDdifferenceObjectiveDf";
}

template <typename T, typename DESCRIPTOR>
bool BlockDdifferenceObjectiveDf3D<T,DESCRIPTOR>::operator()(T dJdF[], const int latticeR[])
{
  for (int i=0; i < DESCRIPTOR::q; i++) {
    dJdF[i] = T();
  }

  if (_indicatorF(latticeR)) {
    const int nDim = _f.getTargetDim();

    std::vector<T> dJdD(nDim,T());

    T f[nDim];
    _f(f, latticeR);

    T wantedF[nDim];
    _wantedF(wantedF, latticeR);

    T dFdF[nDim*DESCRIPTOR::q];
    _dFdF(dFdF,latticeR);

    for (int iDim=0; iDim<nDim; iDim++) {
      dJdD[iDim] = (f[iDim] - wantedF[iDim]) * _weight*_weight;

      for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
        dJdF[jPop] += dJdD[iDim] * dFdF[jPop*nDim + iDim];
      }
    }
  }

  return true;
}

template <typename T, typename DESCRIPTOR>
DdifferenceObjectiveDf3D<T,DESCRIPTOR>::DdifferenceObjectiveDf3D(
  FunctorPtr<SuperLatticePhysF3D<T,DESCRIPTOR>>&& f,
  FunctorPtr<SuperLatticePhysF3D<T,DESCRIPTOR>>&& dFdF,
  FunctorPtr<AnalyticalF3D<T,T>>&&             wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&&           indicatorF)
  : SuperLatticeF3D<T,DESCRIPTOR>(f->getSuperLattice(), f->getTargetDim()),
    _f(std::move(f)),
    _dFdF(std::move(dFdF)),
    _wantedF(std::forward<decltype(wantedF)>(wantedF), f->getSuperLattice()),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "dDifferenceObjectiveDf";

  const T weight = util::pow(_f->getSuperStructure().getCuboidGeometry().getMinDeltaR(), 3);

  for (int iC = 0; iC < this->getSuperLattice().getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockDdifferenceObjectiveDf3D<T,DESCRIPTOR>(
        this->getSuperLattice().getBlock(iC),
        _f->getBlockF(iC),
        _dFdF->getBlockF(iC),
        _wantedF.getBlockF(iC),
        _indicatorF->getBlockIndicatorF(iC),
        weight)
    );
  }
}

template <typename T, typename DESCRIPTOR>
bool DdifferenceObjectiveDf3D<T,DESCRIPTOR>::operator()(T dJdF[], const int latticeR[])
{
  for (int i=0; i < DESCRIPTOR::q; ++i) {
    dJdF[i] = T();
  }

  auto& load = this->getSuperLattice().getLoadBalancer();

  if (load.isLocal(latticeR[0])) {
    return this->getBlockF(load.loc(latticeR[0]))(dJdF, &latticeR[1]);
  }
  else {
    return true;
  }
}


template<typename T, typename DESCRIPTOR>
RelativeDifferenceObjective3D<T,DESCRIPTOR>::RelativeDifferenceObjective3D(
  SuperLattice<T, DESCRIPTOR>&     sLattice,
  FunctorPtr<SuperF3D<T,T>>&&        f,
  FunctorPtr<AnalyticalF3D<T,T>>&&   wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperIdentity3D<T,T>([&]()
{
  using namespace functor_dsl;

  auto wantedLatticeF = restrictF(wantedF.toShared(), sLattice);
  auto relativeDifferenceF = norm<2>(f.toShared() - wantedLatticeF, indicatorF.toShared()) /
                             norm<2>(wantedLatticeF,                indicatorF.toShared());

  return pow(relativeDifferenceF, 2) / T(2.0);
}())
{
  this->getName() = "relativeDifferenceObjective";
}

template<typename T, typename DESCRIPTOR>
RelativeDifferenceObjective3D<T,DESCRIPTOR>::RelativeDifferenceObjective3D(
  FunctorPtr<SuperLatticeF3D<T,DESCRIPTOR>>&& f,
  FunctorPtr<AnalyticalF3D<T,T>>&&            wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&&          indicatorF)
  : RelativeDifferenceObjective3D(
      f->getSuperLattice(),
      std::static_pointer_cast<SuperF3D<T,T>>(f.toShared()),
      std::forward<decltype(wantedF)>(wantedF),
      std::forward<decltype(indicatorF)>(indicatorF))
{ }

template<typename T, typename DESCRIPTOR>
RelativeDifferenceObjective3D<T,DESCRIPTOR>::RelativeDifferenceObjective3D(
  SuperLattice<T, DESCRIPTOR>&   sLattice,
  FunctorPtr<SuperF3D<T,T>>&&      f,
  FunctorPtr<AnalyticalF3D<T,T>>&& wantedF,
  SuperGeometry<T,3>& geometry,
  std::vector<int>    materials)
  : RelativeDifferenceObjective3D(
      sLattice,
      std::forward<decltype(f)>(f),
      std::forward<decltype(wantedF)>(wantedF),
      geometry.getMaterialIndicator(materials))
{ }

template<typename T, typename DESCRIPTOR>
RelativeDifferenceObjective3D<T,DESCRIPTOR>::RelativeDifferenceObjective3D(
  SuperLattice<T, DESCRIPTOR>&   sLattice,
  FunctorPtr<SuperF3D<T,T>>&&      f,
  FunctorPtr<AnalyticalF3D<T,T>>&& wantedF,
  SuperGeometry<T,3>& geometry,
  int                 material)
  : RelativeDifferenceObjective3D(
      sLattice,
      std::forward<decltype(f)>(f),
      std::forward<decltype(wantedF)>(wantedF),
      geometry.getMaterialIndicator(material))
{ }

template<typename T, typename DESCRIPTOR>
RelativeDifferenceObjective3D<T,DESCRIPTOR>::RelativeDifferenceObjective3D(
  SuperLattice<T, DESCRIPTOR>&     sLattice,
  FunctorPtr<SuperF3D<T,T>>&&        f,
  FunctorPtr<SuperF3D<T,T>>&&   wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperIdentity3D<T,T>([&]()
{
  using namespace functor_dsl;

  // auto wantedLatticeF = restrictF(wantedF.toShared(), sLattice);
  auto relativeDifferenceF = norm<2>(f.toShared() - wantedF.toShared(), indicatorF.toShared()) /
                             norm<2>(wantedF.toShared(),                indicatorF.toShared());

  return pow(relativeDifferenceF, 2) / T(2.0);
}())
{
  this->getName() = "relativeDifferenceObjective";
}


template <typename T, typename DESCRIPTOR>
BlockDrelativeDifferenceObjectiveDf3D<T,DESCRIPTOR>::BlockDrelativeDifferenceObjectiveDf3D(
  BlockStructureD<3>& blockStructure,
  BlockF3D<T>&          f,
  BlockF3D<T>&          dFdF,
  BlockF3D<T>&          wantedF,
  BlockIndicatorF3D<T>& indicatorF,
  T globalValue,
  T weight):
  BlockF3D<T>(blockStructure, DESCRIPTOR::q),
  _f(f), _dFdF(dFdF), _wantedF(wantedF), _indicatorF(indicatorF),
  _globalValue(globalValue), _weight(weight)
{
  this->getName() = "BlockDrelativeDifferenceObjectiveDf";
}

template <typename T, typename DESCRIPTOR>
bool BlockDrelativeDifferenceObjectiveDf3D<T,DESCRIPTOR>::operator()(T dJdF[], const int latticeR[])
{
  for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
    dJdF[jPop] = T{};
  }

  if (_indicatorF(latticeR)) {
    const int nDim = DESCRIPTOR::d;

    T dJdD[nDim] {};

    T f[nDim] {};
    _f(f, latticeR);

    T wantedF[nDim] {};
    _wantedF(wantedF, latticeR);

    T dFdF[nDim*DESCRIPTOR::q] {};
    _dFdF(dFdF, latticeR);

    for (int iDim=0; iDim < nDim; ++iDim) {
      dJdD[iDim] = (f[iDim] - wantedF[iDim]) / _globalValue * _weight*_weight;

      for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
        dJdF[jPop] += dJdD[iDim] * dFdF[jPop*nDim + iDim];
      }
    }
  }

  return true;
}

template <typename T, typename DESCRIPTOR>
DrelativeDifferenceObjectiveDf3D<T,DESCRIPTOR>::DrelativeDifferenceObjectiveDf3D(
  SuperLattice<T, DESCRIPTOR>& sLattice,
  FunctorPtr<SuperF3D<T,T>>&&        f,
  FunctorPtr<SuperF3D<T,T>>&&        dFdF,
  FunctorPtr<AnalyticalF3D<T,T>>&&   wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  :  SuperF3D<T,T>(sLattice, DESCRIPTOR::q),
     _f(std::move(f)),
     _dFdF(std::move(dFdF)),
     _wantedF(std::forward<decltype(wantedF)>(wantedF), sLattice),
     _indicatorF(std::move(indicatorF))
{
  this->getName() = "dRelativeDifferenceObjectiveDf";

  SuperL2Norm3D<T> wantedFlp(_wantedF, *_indicatorF);
  T output[1] {};
  wantedFlp(output, nullptr);

  const T globalValue = util::pow(output[0], 2);
  const T weight = util::pow(_f->getSuperStructure().getCuboidGeometry().getMinDeltaR(), 3);

  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockDrelativeDifferenceObjectiveDf3D<T,DESCRIPTOR>(
        sLattice.getBlock(iC),
        _f->getBlockF(iC),
        _dFdF->getBlockF(iC),
        _wantedF.getBlockF(iC),
        _indicatorF->getBlockIndicatorF(iC),
        globalValue,
        weight)
    );
  }
}

template <typename T, typename DESCRIPTOR>
DrelativeDifferenceObjectiveDf3D<T,DESCRIPTOR>::DrelativeDifferenceObjectiveDf3D(
  SuperLattice<T, DESCRIPTOR>& sLattice,
  FunctorPtr<SuperF3D<T,T>>&&      f,
  FunctorPtr<SuperF3D<T,T>>&&      dFdF,
  FunctorPtr<AnalyticalF3D<T,T>>&& wantedF,
  SuperGeometry<T,3>& geometry,
  std::vector<int>    materials)
  : DrelativeDifferenceObjectiveDf3D(sLattice,
                                     std::forward<decltype(f)>(f),
                                     std::forward<decltype(dFdF)>(dFdF),
                                     std::forward<decltype(wantedF)>(wantedF),
                                     geometry.getMaterialIndicator(materials))
{ }

template <typename T, typename DESCRIPTOR>
DrelativeDifferenceObjectiveDf3D<T,DESCRIPTOR>::DrelativeDifferenceObjectiveDf3D(
  SuperLattice<T, DESCRIPTOR>& sLattice,
  FunctorPtr<SuperF3D<T,T>>&&      f,
  FunctorPtr<SuperF3D<T,T>>&&      dFdF,
  FunctorPtr<AnalyticalF3D<T,T>>&& wantedF,
  SuperGeometry<T,3>& geometry,
  int                 material)
  : DrelativeDifferenceObjectiveDf3D(sLattice,
                                     std::forward<decltype(f)>(f),
                                     std::forward<decltype(dFdF)>(dFdF),
                                     std::forward<decltype(wantedF)>(wantedF),
                                     geometry.getMaterialIndicator(material))
{ }

template <typename T, typename DESCRIPTOR>
bool DrelativeDifferenceObjectiveDf3D<T,DESCRIPTOR>::operator()(T dJdF[], const int latticeR[])
{
  for (int i=0; i < DESCRIPTOR::q; ++i) {
    dJdF[i] = T{};
  }

  auto& load = _f->getSuperStructure().getLoadBalancer();

  if (load.isLocal(latticeR[0])) {
    return this->getBlockF(load.loc(latticeR[0]))(dJdF, &latticeR[1]);
  }
  else {
    return true;
  }
}

template <typename T, typename DESCRIPTOR>
DrelativeDifferenceObjectiveDf3D_Lattice<T,DESCRIPTOR>::DrelativeDifferenceObjectiveDf3D_Lattice(
  SuperLattice<T, DESCRIPTOR>& sLattice,
  FunctorPtr<SuperF3D<T,T>>&&        f,
  FunctorPtr<SuperF3D<T,T>>&&        dFdF,
  FunctorPtr<SuperF3D<T,T>>&&        wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  :  SuperF3D<T,T>(sLattice, DESCRIPTOR::q),
     _f(std::move(f)),
     _dFdF(std::move(dFdF)),
     _wantedF(std::move(wantedF)),
     _indicatorF(std::move(indicatorF))
{
  this->getName() = "dRelativeDifferenceObjectiveDf_Lattice";

  SuperL2Norm3D<T> wantedFlp(*_wantedF, *_indicatorF);
  T output[1] {};
  wantedFlp(output, nullptr);

  const T globalValue = util::pow(output[0], 2);
  const T weight = util::pow(sLattice.getCuboidGeometry().getMinDeltaR(), 3);
  // util::pow(_f->getSuperStructure().getCuboidGeometry().getMinDeltaR(), 3);

  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockDrelativeDifferenceObjectiveDf3D<T,DESCRIPTOR>(
        sLattice.getBlock(iC),
        _f->getBlockF(iC),
        _dFdF->getBlockF(iC),
        _wantedF->getBlockF(iC),
        _indicatorF->getBlockIndicatorF(iC),
        globalValue,
        weight)
    );
  }
}

template <typename T, typename DESCRIPTOR>
bool DrelativeDifferenceObjectiveDf3D_Lattice<T,DESCRIPTOR>::operator()(T dJdF[], const int latticeR[])
{
  for (int i=0; i < DESCRIPTOR::q; ++i) {
    dJdF[i] = T{};
  }

  auto& load = _f->getSuperStructure().getLoadBalancer();

  if (load.isLocal(latticeR[0])) {
    return this->getBlockF(load.loc(latticeR[0]))(dJdF, &latticeR[1]);
  }
  else {
    return true;
  }
}


template <typename T, typename DESCRIPTOR>
BlockDrelativeDifferenceObjectiveComponentDf3D<T,DESCRIPTOR>::BlockDrelativeDifferenceObjectiveComponentDf3D(
  BlockStructureD<3>& blockStructure,
  BlockExtractComponentF3D<T>& f,
  BlockF3D<T>&                 dFdF,
  BlockF3D<T>&                 wantedF,
  BlockIndicatorF3D<T>&        indicatorF,
  T globalValue,
  T weight):
  BlockF3D<T>(blockStructure, DESCRIPTOR::q),
  _f(f), _dFdF(dFdF), _wantedF(wantedF), _indicatorF(indicatorF),
  _extractDim(f.getExtractDim()), _globalValue(globalValue), _weight(weight)
{
  this->getName() = "BlockDrelativeDifferenceObjectiveComponentDf";
}

template <typename T, typename DESCRIPTOR>
bool BlockDrelativeDifferenceObjectiveComponentDf3D<T,DESCRIPTOR>::operator()(T dJdF[], const int latticeR[])
{
  for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
    dJdF[jPop] = T{};
  }

  if (_indicatorF(latticeR)) {
    const int nDim = DESCRIPTOR::d;

    T dJdD[nDim] {};

    T f[1] {};
    _f(f, latticeR);

    T wantedF[nDim] {};
    _wantedF(wantedF, latticeR);

    T dFdF[nDim*DESCRIPTOR::q] {};
    _dFdF(dFdF, latticeR);

    for (int iDim=0; iDim < nDim; ++iDim) {
      if (iDim == _extractDim) {
        // Note that by convention wantedF is always one-dimensional.
        // This is why f[0] - wantedF[0] is correct (as compared to f[0] - wantedF[iDim])
        dJdD[iDim] = (f[0] - wantedF[0]) / _globalValue * _weight*_weight;
      }

      for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
        dJdF[jPop] += dJdD[iDim] * dFdF[jPop*nDim + iDim];
      }
    }
  }

  return true;
}

template <typename T, typename DESCRIPTOR>
DrelativeDifferenceObjectiveComponentDf3D<T,DESCRIPTOR>::DrelativeDifferenceObjectiveComponentDf3D(
  SuperLattice<T, DESCRIPTOR>& sLattice,
  FunctorPtr<SuperF3D<T,T>>&&        f, int extractDim,
  FunctorPtr<SuperF3D<T,T>>&&        dFdF,
  FunctorPtr<AnalyticalF3D<T,T>>&&   wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  :  SuperF3D<T,T>(sLattice, DESCRIPTOR::q),
     _f(std::forward<decltype(f)>(f), extractDim),
     _dFdF(std::move(dFdF)),
     _wantedF(std::forward<decltype(wantedF)>(wantedF), sLattice),
     _indicatorF(std::move(indicatorF))
{
  this->getName() = "dRelativeDifferenceObjectiveDf";

  SuperL2Norm3D<T> wantedFlp(_wantedF, *_indicatorF);
  T output[1] {};
  wantedFlp(output, nullptr);

  const T globalValue = util::pow(output[0], 2);
  const T weight = util::pow(_f.getSuperStructure().getCuboidGeometry().getMinDeltaR(), 3);

  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockDrelativeDifferenceObjectiveComponentDf3D<T,DESCRIPTOR>(
        sLattice.getBlock(iC),
        static_cast<BlockExtractComponentF3D<T>&>(_f.getBlockF(iC)),
        _dFdF->getBlockF(iC),
        _wantedF.getBlockF(iC),
        _indicatorF->getBlockIndicatorF(iC),
        globalValue,
        weight)
    );
  }
}

template <typename T, typename DESCRIPTOR>
bool DrelativeDifferenceObjectiveComponentDf3D<T,DESCRIPTOR>::operator()(T dJdF[], const int latticeR[])
{
  for (int i=0; i < DESCRIPTOR::q; ++i) {
    dJdF[i] = T{};
  }

  auto& load = _f.getSuperStructure().getLoadBalancer();

  if (load.isLocal(latticeR[0])) {
    return this->getBlockF(load.loc(latticeR[0]))(dJdF, &latticeR[1]);
  }
  else {
    return true;
  }
}

} // namespace opti

} // namespace olb

#endif
