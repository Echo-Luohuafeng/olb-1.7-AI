/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Marc Haußmann, Mathias J. Krause, Jonathan Jeppener-Haltenhoff
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

/* poiseuille3d.cpp:
 * This example examines a 3D Poseuille flow
 * It illustrates the usage of different flow setups, boundary conditions
 * and computation of error norms.
 * The following simulation options can be combined freely:
 * - if the compiler flag ENABLE_MRT is set, mrt collision operators are used
 * - forced/ nonForced flow
 * - different boundary conditions
 * - simulation only or eoc convergence analysis
 */

//#define ENABLE_MRT

#include "poiseuille3d.h"

/// User dialogue: read optional arguments
// Resolution, flow and boundary type are modified in place
int readParameters(int argc, char** argv,
  int& N, FlowType& flowType, BoundaryType& boundaryType, bool& eoc)
{
  if (argc > 1) {
    if (argv[1][0]=='-'&&argv[1][1]=='h') {
      OstreamManager clout( std::cout,"help" );
      clout<<"Usage: program [Resolution] [FlowType] [BoundaryType] [Eoc]"<<std::endl;
      clout<<"FlowType: 0=forced, 1=nonForced"<<std::endl;
      clout<<"BoundaryType: 0=bounceBack, 1=local, 2=interpolated, 3=bouzidi, 4=freeSlip, 5=partialSlip"<<std::endl;
      clout<<"Eoc: 0=false, 1=true"<<std::endl;
      clout<<"Default: Resolution=21, FlowType=nonForced, BoundaryType=interpolated, Eoc=false"<<std::endl;
      return 0;
    }
  }

  if (argc > 1) {
    N = atoi(argv[1]);
    if (N < 1) {
      std::cerr << "Fluid domain is too small" << std::endl;
      return 1;
    }
  }

  if (argc > 2) {
    int flowTypeNumber = atoi(argv[2]);
    if (flowTypeNumber < 0 || flowTypeNumber > (int)nonForced) {
      std::cerr << "Unknown fluid flow type" << std::endl;
      return 2;
    }
    flowType = (FlowType) flowTypeNumber;
  }

  if (argc > 3) {
    int boundaryTypeNumber = atoi(argv[3]);
    if (boundaryTypeNumber < 0 || boundaryTypeNumber > (int) partialSlip) {
      std::cerr << "Unknown boundary type" << std::endl;
      return 3;
    }
    boundaryType = (BoundaryType) boundaryTypeNumber;
  }

  if (argc > 4) {
    int eocNumber = atoi(argv[4]);
    eoc = (bool) eocNumber;
  }
  return 0;
}


int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  int status = readParameters(argc, argv, N, flowType, boundaryType, eoc);
  if (status != 0) {
    return status;
  }

  if (! eoc) {
    static Gnuplot<T> gplot("centerVelocity");
    simulatePoiseuille(N, gplot, eoc);
  }
  else {
    if ((boundaryType == freeSlip) || (boundaryType == partialSlip)) {
      throw std::invalid_argument(
        "eoc computation is currently not supported for slip boundary conditions");
    }
    static Gnuplot<T> gplot(
      "Velocity_and_StrainRate_eoc",
      false,
      "set terminal png size 720, 720 font 'Arial,10'",
      Gnuplot<T>::LOGLOG,
      Gnuplot<T>::LINREG);
    gplot.setLabel("Resolution", "Error");

    int maxN = N + 41;

    // loop over the different simulations
    for(int simuN = N; simuN < maxN; simuN += 10){

      /// Run the simulations
      clout << "Starting next simulation with N = " << simuN << std::endl;
      simulatePoiseuille(simuN, gplot, eoc);
    }

    gplot.writePNG();
  }
}
