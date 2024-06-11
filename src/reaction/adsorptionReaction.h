/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2022 Florian Raichle
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

#ifndef OLB_APPS_FLORIAN_ADSORPTION3D_ADSORPTIONREACTION_H_
#define OLB_APPS_FLORIAN_ADSORPTION3D_ADSORPTIONREACTION_H_

#include <vector>
#include <cmath>
#include <iostream>
#include "core/adsorptionConverter.h"
#include "isotherm.h"

namespace olb {

struct SCALAR2 : public descriptors::FIELD_BASE<1,  0, 0> { };

/**
 * @brief Describes adsorption reactions in conjunction with a Isotherm class.
 *
 * @tparam T
 * @tparam DESCRIPTOR
 *
 * @see Isotherm
 */
template<typename T, typename DESCRIPTOR>
class AdsorptionReaction {
 protected:
//  virtual void boundary(int n, T x[], T fvec[], int &iflag, T args[]) = 0;
  bool externalMassTransferEnabled = false;
  const AdsorptionConverter<T, DESCRIPTOR> _unitConverter;
  const Isotherm<T>* isotherm; // needs to be pointer because Isotherm is an abstract class

 public:
  T k_f{}; // film diffusion mass transfer coefficient (m/s)
  T D_s{}; // surface diffusion constant (m^2/s)
  T c_0{}; // initial concentration (mg/L)
  T r_p{}; // particle radius (m)
  T q_0{}; // equilibrium surface loading (mg/g)
  T k_s{}; // surface diffusion mass transfer coefficient (1/s)

//  explicit AdsorptionReaction(void(*fcn)(int, T[], T[], int&, T[])): fcn(fcn){}
  AdsorptionReaction(AdsorptionConverter<T, DESCRIPTOR> const& converter, Isotherm<T> const* isotherm, T k_f, T D_s, T c_0, T r_p):
  _unitConverter(converter),
  isotherm(isotherm)
  {
    // calculate lattice values
    T conversionFactorLength = _unitConverter.getConversionFactorLength();
    this->k_f = k_f * _unitConverter.getConversionFactorTime() * 3 / r_p;
    this->D_s = D_s / (conversionFactorLength*conversionFactorLength) * _unitConverter.getConversionFactorTime();
    this->c_0 = c_0 / _unitConverter.getConversionFactorDensity();
    this->r_p = r_p / conversionFactorLength;
    this->k_s = 15 * this->D_s / (this->r_p * this->r_p);
    this->q_0 = this->isotherm->getLoading(c_0) / _unitConverter.getConversionFactorParticleDensity();
    externalMassTransferEnabled = this->k_f > 0.;
  }

  /**
   * @brief Get surface loading without film diffusion.
   *
   * @param soluteConcentration
   * @return loading in g/m^3
   */
  T getSurfaceLoading(T soluteConcentration) {
    if (soluteConcentration < 0) return 0; // prevent nan
    T load = this->isotherm->getLoading(soluteConcentration);
    if (load < 0) load = 0;
    return load;
  }

  /**
   * @brief  Get surface loading with film diffusion.
   *
   * Solves implicit boundary condition for surface concentration
   *
   * @param soluteConcentration
   * @param particleLoading
   * @param particleConcentration
   * @return T
   */
  T getSurfaceLoading(T soluteConcentration, T particleLoading, T particleConcentration) {
    if (soluteConcentration < 0) return 0; // prevent nan
    // using D = double;  // because fsolve is hardcoded in double type

    // previous surface concentration as initial guess
    T surfaceConcentration = soluteConcentration;
    T args[] = {particleLoading, soluteConcentration, particleConcentration};
    for(int i = 0; i<100; i++){
      //ADf<T,1> y0_ad = iniAD(surfaceConcentration);
      //ADf<T,1> y1_ad = filmDiffBoundary(y0_ad, args);
      //surfaceConcentration -= filmDiffBoundary(surfaceConcentration, args) / y1_ad.d(0);

      T dCs = surfaceConcentration/1E8;
      surfaceConcentration -= 2.*dCs*filmDiffBoundary(surfaceConcentration, args)/(filmDiffBoundary(surfaceConcentration + dCs, args) - filmDiffBoundary(surfaceConcentration - dCs, args));
    }
    T load = isotherm->getLoading(surfaceConcentration);
    if (load < 0) load = 0;
    return load;
  }

  /// Calculates reaction rate in lattice units
  /// \param soluteConcentration
  /// \param particleLoading
  /// \param particleConcentration
  /// \return change in solute concentration
  Vector<T, 2> getReactionRate(T soluteConcentration, T particleLoading, T particleConcentration) {
    T surfaceLoad;
    particleConcentration *= _unitConverter.getConversionFactorParticleDensity(); // convert to kg/m^3
    soluteConcentration *= _unitConverter.getConversionFactorDensity();           // convert to g/m^3
    particleLoading *= _unitConverter.getConversionFactorParticleDensity();       // convert to g/m^3
    if (externalMassTransferEnabled) {
      surfaceLoad = this->getSurfaceLoading(soluteConcentration, particleLoading, particleConcentration);
    } else {
      surfaceLoad = this->getSurfaceLoading(soluteConcentration);
    }

    T reactionRate = this->k_s * (surfaceLoad * particleConcentration - particleLoading);
    Vector<T, 2> reactionRates(reactionRate / _unitConverter.getConversionFactorDensity(),          // solute
                               reactionRate / _unitConverter.getConversionFactorParticleDensity()); // loading
    return reactionRates;
  }

  T getPhysFilmTransferConstant() {
    return this->k_f / _unitConverter.getConversionFactorTime();
  }
  T getPhysSurfaceTransferConstant() {
    return this->k_s / _unitConverter.getConversionFactorTime();
  }

  T getPhysSurfaceDiffusionConstant() {
    return this->D_s * _unitConverter.getConversionFactorLength() * _unitConverter.getConversionFactorLength() / _unitConverter.getConversionFactorTime();
  }

  void print(std::ostream& clout) {
    clout << "----------------- Reaction information -----------------" << std::endl;
    clout << "-- Parameters:" << std::endl;
    clout << "Particle diameter(m):                   d_p=    " << this->r_p*2 * _unitConverter.getConversionFactorLength() << std::endl;
    clout << "Film diffusion constant(m/s):           k_f*=   " << getPhysFilmTransferConstant() << std::endl;
    clout << "Surface diffusion constant(m^2/s):      D_s=    " << getPhysSurfaceDiffusionConstant() << std::endl;
    clout << "Initial solute concentration(mg/mL):    c_0=    " << this->c_0 * _unitConverter.getConversionFactorDensity() << std::endl;
    clout << "Equilibrium surface loading(mg/g):      q_0=    " << this->q_0 * _unitConverter.getConversionFactorParticleDensity() << std::endl;
    clout << "lattice Equilibrium surface loading:    q_0=    " << this->q_0 << std::endl;
    clout << "Surface Mass transfer coefficient(m/s): k_s*=   " << getPhysSurfaceTransferConstant() << std::endl;
    clout << "-------------------------------------------------------------" << std::endl;
    this->isotherm->print(clout);
  }

  void write(std::string const& fileName) {
    std::string dataFile = singleton::directories().getLogOutDir() + fileName + ".dat";

    if (singleton::mpi().isMainProcessor()) {
      std::ofstream fout(dataFile.c_str(), std::ios::app);
      if (!fout) {
        std::cout << "error write() function: can not open std::ofstream" << std::endl;
      }
      else {
        this->print( fout );
        fout.close();
      }
    }
  }

  /// equations for surface concentration
  T filmDiffBoundary(T c_s_, T args[]) {
    T q = args[0];
    T c = args[1];
    T particleConcentration = args[2];
    auto c_s = c_s_;
    T q_s = this->isotherm->getLoading(c_s);
    T equation = this->k_f * (c - c_s) - this->k_s * (particleConcentration * q_s - q);
    return equation;
  }
};
}
#endif //OLB_APPS_FLORIAN_ADSORPTION3D_ADSORPTIONREACTION_H_
