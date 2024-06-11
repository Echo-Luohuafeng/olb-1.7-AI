/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Claudius Holeksa
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

#ifndef FREE_SURFACE_POST_PROCESSOR_3D_HH
#define FREE_SURFACE_POST_PROCESSOR_3D_HH

#include "freeSurfacePostProcessor3D.h"
#include "core/blockLattice.h"

#include <cmath>
#include <algorithm>

namespace olb {

// LocalProcessor 1

// Read

// range 0
// Mass
// CellType

// range 1
// DFs
// Epsilon
// CellType

// Write (always range 0)
// Mass
// CellFlags
// DFs (only replacing from incoming gas stream)

template <typename T, typename DESCRIPTOR>
template <typename CELL, typename PARAMETERS>
void FreeSurfaceMassFlowPostProcessor3D<T, DESCRIPTOR>::apply(CELL& cell, PARAMETERS& vars) {

  using namespace olb::FreeSurface;

  const bool drop_isolated_cells = vars.template get<FreeSurface::DROP_ISOLATED_CELLS>();
  const T transition = vars.template get<FreeSurface::TRANSITION>();
  const T lonely_threshold = vars.template get<FreeSurface::LONELY_THRESHOLD>();
  const bool has_surface_tension = vars.template get<FreeSurface::HAS_SURFACE_TENSION>();
  const T surface_tension_parameter = vars.template get<FreeSurface::SURFACE_TENSION_PARAMETER>();
  //const T force_conversion_factor = vars.template get<FreeSurface::FORCE_CONVERSION_FACTOR>();
  //const T lattice_size = vars.template get<FreeSurface::LATTICE_SIZE>();

  /*
  * Minor "hack". Remove all cell flags here, because it is needed in the last processor due to pulling steps in processor 6 and 7
  */
  setCellFlags(cell, static_cast<FreeSurface::Flags>(0));

  /*
  * This processor only works on interface types
  */
  if (isCellType(cell, FreeSurface::Type::Interface )) {
    T mass_tmp = cell.template getField<FreeSurface::MASS>();

    FreeSurface::NeighbourInfo neighbour_info = getNeighbourInfo(cell);

    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
      auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop, 0),
                                  descriptors::c<DESCRIPTOR>(iPop, 1),
                                  descriptors::c<DESCRIPTOR>(iPop, 2)});
      int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);

      /*
      * Iterate over neighbours and perform a mass exchange. Interface to fluid are simple cases.
      * Interface to interface has to be symmetric and multiple cases are for artifact reduction
      * from Thuerey
      * Added a distinction between the amount of interface nodes. Weight consideration seems to cause artifacts.
      */
      if( isCellType(cellC, FreeSurface::Type::Fluid)) {
        mass_tmp += cell[iPop_op] - cellC[iPop];
      } else if ( isCellType(cellC, FreeSurface::Type::Interface)) {
        FreeSurface::NeighbourInfo neighbour_neighbour_info = getNeighbourInfo(cellC);

        T mass_flow = 0.;

        if( !neighbour_info.has_fluid_neighbours){
          if(!neighbour_neighbour_info.has_fluid_neighbours){
            if(neighbour_info.interface_neighbours < neighbour_neighbour_info.interface_neighbours){
              mass_flow = -cellC[iPop];// - descriptors::t<T,DESCRIPTOR>(iPop);
            }else if(neighbour_info.interface_neighbours > neighbour_neighbour_info.interface_neighbours){
              mass_flow = cell[iPop_op];// + descriptors::t<T,DESCRIPTOR>(iPop_op);
            }else{
              mass_flow = cell[iPop_op] - cellC[iPop];
            }
          }else {
            mass_flow = -cellC[iPop];// - descriptors::t<T,DESCRIPTOR>(iPop);
          }
        }else if(!neighbour_info.has_gas_neighbours){
          if(!neighbour_neighbour_info.has_gas_neighbours){
            if(neighbour_info.interface_neighbours < neighbour_neighbour_info.interface_neighbours){
              mass_flow = cell[iPop_op];// + descriptors::t<T,DESCRIPTOR>(iPop_op);
            }else if(neighbour_info.interface_neighbours > neighbour_neighbour_info.interface_neighbours){
              mass_flow = -cellC[iPop];// - descriptors::t<T,DESCRIPTOR>(iPop);
            }else{
              mass_flow = cell[iPop_op] - cellC[iPop];
            }
          }else {
            mass_flow = cell[iPop_op];// + descriptors::t<T,DESCRIPTOR>(iPop_op);
          }
        }else {
          if(!neighbour_neighbour_info.has_fluid_neighbours){
            mass_flow = cell[iPop_op];// + descriptors::t<T,DESCRIPTOR>(iPop_op);
          }else if(!neighbour_neighbour_info.has_gas_neighbours){
            mass_flow = -cellC[iPop];// - descriptors::t<T,DESCRIPTOR>(iPop);
          }else {
            mass_flow = cell[iPop_op] - cellC[iPop];
          }
        }

        /*
        * Exchange depends on how filled the interfaces are
        */
        mass_tmp += mass_flow * 0.5 * (getClampedEpsilon(cell) + getClampedEpsilon(cellC));
      }
    }

    cell.template setField<FreeSurface::MASS>(mass_tmp);

    // Former 2 Step

    // Because I need the distribution functions of the last step I will write results in a temporary
    // array, before copying it back into the DFs

    Vector<T,DESCRIPTOR::q> dfs;

    T curvature = 0.;

    if(has_surface_tension){
      FreeSurface::NeighbourInfo info = getNeighbourInfo(cell);
      if(info.has_gas_neighbours){
        curvature = calculateSurfaceTensionCurvature(cell);
      }

    }

    // Gas pressure adjusting
    T gas_pressure = 1. - 6. * surface_tension_parameter * curvature;

    for(int iPop=1; iPop < DESCRIPTOR::q; iPop++) {
      auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop, 0),
                                  descriptors::c<DESCRIPTOR>(iPop, 1),
                                  descriptors::c<DESCRIPTOR>(iPop, 2)});
      int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);

      /*
      * Gas replacement
      */
      if ( isCellType(cellC, FreeSurface::Type::Gas )) {
        Vector<T, DESCRIPTOR::d> u_vel = cell.template getField<FreeSurface::PREVIOUS_VELOCITY>();
        T u[DESCRIPTOR::d];
        for(size_t u_i = 0; u_i < DESCRIPTOR::d; ++u_i){
          u[u_i] = u_vel[u_i];
        }
        //T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
        dfs[iPop_op] = equilibrium<DESCRIPTOR>::secondOrder(iPop, gas_pressure, u)
          + equilibrium<DESCRIPTOR>::secondOrder(iPop_op, gas_pressure, u)
          - cellC[iPop];
      }else {
        dfs[iPop_op] = cell[iPop_op];
      }
    }

    for(int iPop=1; iPop<DESCRIPTOR::q; iPop++) {
      cell[iPop] = dfs[iPop];
    }

    // Former 3 Step
    /*
    * Based on the mass calculation, flag this interface cell as toFluid or toGas if set boundaries are met
    */
    T rho = cell.computeRho();

    // Check if transition needs to happen.
    if ( mass_tmp < -transition * rho || (mass_tmp < lonely_threshold * rho && !neighbour_info.has_fluid_neighbours) ){
      setCellFlags(cell, FreeSurface::Flags::ToGas);
    }
    else if ( mass_tmp > (1. + transition)*rho  || ( mass_tmp > (1-lonely_threshold) * rho && !neighbour_info.has_gas_neighbours) ){
      setCellFlags(cell, FreeSurface::Flags::ToFluid);
    }else if(drop_isolated_cells && (neighbour_info.interface_neighbours == 0)){
      if(!neighbour_info.has_gas_neighbours){
        setCellFlags(cell, FreeSurface::Flags::ToFluid);
      }else if(!neighbour_info.has_fluid_neighbours){
        //setCellFlags(cell, FreeSurface::Flags::ToGas);
      }
    }
  }
}

// Processor 4

// Read
// range 0
// CellType
// CellFlags

// range 1
// CellType
// CellFlags
// DFs

// Write (always range 0)
// CellFlags
// DFs
template <typename T, typename DESCRIPTOR>
template <typename CELL>
void FreeSurfaceToFluidCellConversionPostProcessor3D<T, DESCRIPTOR>::apply(CELL& cell) {
  using namespace olb::FreeSurface;
  /*
  * Initializing new interface cells with DFs from surrounding fluid and interface cells
  * Just takes the arithmetic average.
  */
  if(isCellType(cell, FreeSurface::Type::Gas)){
    if(hasNeighbourFlags(cell, FreeSurface::Flags::ToFluid)){
      setCellFlags(cell, FreeSurface::Flags::NewInterface);
      T rho_avg = 0.;
      T u_avg[DESCRIPTOR::d] = {0., 0.};

      std::size_t ctr = 0;

      for(int iPop=1; iPop<DESCRIPTOR::q; iPop++) {
        auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop,0),
                                    descriptors::c<DESCRIPTOR>(iPop,1),
                                    descriptors::c<DESCRIPTOR>(iPop,2)});

        if (isCellType(cellC, FreeSurface::Type::Fluid) || isCellType(cellC, FreeSurface::Type::Interface)){
          T rho_tmp = 0.;
          T u_tmp[DESCRIPTOR::d] = {0., 0.};
          ++ctr;
          cellC.computeRhoU(rho_tmp, u_tmp);
          rho_avg += rho_tmp;
          for(size_t i = 0; i < DESCRIPTOR::d; ++i){
            u_avg[i] += u_tmp[i];
          }
        }
      }

      if(ctr > 0){
        rho_avg /= static_cast<T>(ctr);
        for(size_t i = 0; i < DESCRIPTOR::d; ++i){
          u_avg[i] /= static_cast<T>(ctr);
        }
      }

      cell.iniEquilibrium(rho_avg, u_avg);
    }
  }else if(hasCellFlags(cell, FreeSurface::Flags::ToGas)){
    /*
    * If a toGas cell has a neighbouring toFluid cell, unset the toGas flag
    */
    if(hasNeighbourFlags(cell, FreeSurface::Flags::ToFluid)){
      setCellFlags(cell, static_cast<FreeSurface::Flags>(0));
    }
  }
}

// LocalProcessor 5

// Read
// range 0
// CellType
// DFs

// range 1
// CellFlags

// Write (always range 0)
// CellFlags

template <typename T, typename DESCRIPTOR>
template <typename CELL>
void FreeSurfaceToGasCellConversionPostProcessor3D<T, DESCRIPTOR>::apply(CELL& cell) {
  using namespace olb::FreeSurface;
  /*
  * For the to be converted toGas cells, set the neighbours to interface cells
  */
  if (   isCellType(cell, FreeSurface::Type::Fluid)
      && hasNeighbourFlags(cell, FreeSurface::Flags::ToGas)){
    setCellFlags(cell, FreeSurface::Flags::NewInterface);
    T rho = cell.computeRho();
    cell.template setField<FreeSurface::MASS>(rho);
  }
}

// LocalProcessor 6

// Read
// range 0
// CellType
// CellFlags
// DFs
// Mass

// range 1
// Epsilon
// CellType
// CellFlags

// Write (always range 0)
// Mass
// TempMassExchange
template <typename T, typename DESCRIPTOR>
template <typename CELL>
void FreeSurfaceMassExcessPostProcessor3D<T, DESCRIPTOR>::apply(CELL& cell) {
  using namespace olb::FreeSurface;
  if( !isCellType(cell, FreeSurface::Type::Interface) ){
    return;
  }

  T rho = cell.computeRho();
  T mass = cell.template getField<FreeSurface::MASS>( );
  T mass_excess = 0.;

  auto normal = computeParkerYoungInterfaceNormal(cell);
  // redistribute excess mass

  /// @hint EPSILON of neighbours used here
  /// @hint Mass can be set in this processor, but not epsilon since it is needed for the normal computation. epsilon is set in the next processor
  /// Became untrue due to code section removal, but epsilon is still set in the next part because of pushed mass excess
  if(hasCellFlags(cell, FreeSurface::Flags::ToGas)){
    mass_excess = mass;
    cell.template setField<FreeSurface::MASS>( 0. );
    normal *= -1;
  } else if (hasCellFlags(cell, FreeSurface::Flags::ToFluid)) {
    mass_excess = mass - rho;
    cell.template setField<FreeSurface::MASS>( rho );
  } else {
    return;
  }

  std::array<T,DESCRIPTOR::q> products;
  products[0] = 0.;
  T product_sum = 0.;
  std::size_t product_total = 0;

  for(int iPop=1; iPop<DESCRIPTOR::q; iPop++) {
    auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop,0),
                                descriptors::c<DESCRIPTOR>(iPop,1),
                                descriptors::c<DESCRIPTOR>(iPop,2)});
    products[iPop] = 0.;

    // Thuerey Paper says we can't use new interface cells
    // or flagged cells
    // But surface tension showed us that it has anisotropic effects
    if( (isCellType(cellC, FreeSurface::Type::Interface) && (!hasCellFlags(cellC,
          static_cast<FreeSurface::Flags>(255)) /*|| hasCellFlags(cellC, FreeSurface::Flags::ToFluid)*/ ))
        /*|| isCellType(cellC, FreeSurface::Type::Fluid)*/
        ){
      products[iPop] = (normal[0] * descriptors::c<DESCRIPTOR>(iPop, 0) + normal[1] * descriptors::c<DESCRIPTOR>(iPop,1) + normal[2] * descriptors::c<DESCRIPTOR>(iPop,2));
      if(products[iPop] <= 0){
        products[iPop] = 0.;
      }
      ++product_total;
      product_sum += products[iPop];
    }
  }

  /* Prepare Mass excess push */
  Vector<T,DESCRIPTOR::q> mass_excess_vector{};
  mass_excess_vector[0] = 0.;
  /*
  if(product_sum > 0){
    T fraction = 1./ product_sum;

    for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
      mass_excess_vector[iPop] = mass_excess * products[iPop] * fraction;
    }
    cell.template setField<FreeSurface::TEMP_MASS_EXCHANGE>( mass_excess_vector );
  }
  else*/
  if (product_total > 0) {
    T product_fraction = 1. / product_total;
    for(int iPop=1; iPop < DESCRIPTOR::q; iPop++) {
      auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop,0),
                                  descriptors::c<DESCRIPTOR>(iPop,1),
                                  descriptors::c<DESCRIPTOR>(iPop,2)});
      if( (isCellType(cellC, FreeSurface::Type::Interface) && (!hasCellFlags(cellC,
          static_cast<FreeSurface::Flags>(255)) /*|| hasCellFlags(cellC, FreeSurface::Flags::ToFluid)*/ ))
          /*|| isCellType(cellC, FreeSurface::Type::Fluid)*/
          ){
        mass_excess_vector[iPop] = mass_excess * product_fraction;
      } else {
        mass_excess_vector[iPop] = 0.;
      }
    }
    cell.template setField<FreeSurface::TEMP_MASS_EXCHANGE>( mass_excess_vector );
  } else {
    mass_excess_vector[0] = mass_excess;
    for(int iPop=1; iPop < DESCRIPTOR::q; iPop++) {
      mass_excess_vector[iPop] = 0.;
    }
    cell.template setField<FreeSurface::TEMP_MASS_EXCHANGE>( mass_excess_vector );
  }
}


// LocalProcessor 7

// Read
// range 0
// CellFlags
// CellType
// DFs
// Mass
// Epsilon

// range 1
// CellFlags
// TempMassExchange

// Write (always range 0)
// Epsilon
// CellType
// Mass

template <typename T, typename DESCRIPTOR>
template <typename CELL>
void FreeSurfaceFinalizeConversionPostProcessor3D<T, DESCRIPTOR>::apply(CELL& cell) {
  using namespace olb::FreeSurface;
  /* Convert flagged cells to appropriate cell types */
  FreeSurface::Flags flags = static_cast<FreeSurface::Flags>(cell.template getField<FreeSurface::CELL_FLAGS>());

  switch(flags){
    case FreeSurface::Flags::ToFluid:
    {
      /// @hint moved flag removal to processor 1 without any negative effects
      setCellType(cell, FreeSurface::Type::Fluid);
      cell.template setField<FreeSurface::EPSILON>( 1. );
      T mass_tmp = cell.template getField<FreeSurface::MASS>();
      mass_tmp += cell.template getField<FreeSurface::TEMP_MASS_EXCHANGE>()[0];
      cell.template setField<FreeSurface::MASS>(mass_tmp);

    }
    break;
    case FreeSurface::Flags::ToGas:
    {
      /// @hint moved flag removal to processor 1 without any negative effects
      setCellType(cell, FreeSurface::Type::Gas);
      cell.template setField<FreeSurface::EPSILON>( 0. );
      T mass_tmp = cell.template getField<FreeSurface::MASS>();
      mass_tmp += cell.template getField<FreeSurface::TEMP_MASS_EXCHANGE>()[0];
      cell.template setField<FreeSurface::MASS>(mass_tmp);

    }
    break;
    case FreeSurface::Flags::NewInterface:
    {
      setCellType(cell, FreeSurface::Type::Interface);
    }
    break;
    default:
    break;
  }

  FreeSurface::Type type = static_cast<FreeSurface::Type>(cell.template getField<FreeSurface::CELL_TYPE>());

  /* Collection of mass excess in a pulling step */
  switch(type){
    case FreeSurface::Type::Interface:
    {
      T collected_excess = 0.;
      for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
        auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop,0),
                                    descriptors::c<DESCRIPTOR>(iPop,1),
                                    descriptors::c<DESCRIPTOR>(iPop,2)});
        auto tempMassExchange = cellC.template getFieldPointer<FreeSurface::TEMP_MASS_EXCHANGE>();

        if (hasCellFlags(cellC,
                FreeSurface::Flags::ToFluid
              | FreeSurface::Flags::ToGas)){
            int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);
            collected_excess += tempMassExchange[iPop_op];
          }
      }

      T mass_tmp = cell.template getField<FreeSurface::MASS>();

      mass_tmp += collected_excess;
      T rho;
      T u_tmp[DESCRIPTOR::d];
      cell.computeRhoU(rho, u_tmp);

      Vector<T,DESCRIPTOR::d> u_vel{u_tmp[0], u_tmp[1], u_tmp[2]};

      cell.template setField<FreeSurface::EPSILON>( mass_tmp / rho );
      cell.template setField<FreeSurface::MASS>(mass_tmp);
      cell.template setField<FreeSurface::PREVIOUS_VELOCITY>(u_vel);
    }
    break;
    case FreeSurface::Type::Fluid:
    {
      T collected_excess = 0.;

      for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
        auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop,0),
                                    descriptors::c<DESCRIPTOR>(iPop,1),
                                    descriptors::c<DESCRIPTOR>(iPop,2)});
        auto tempMassExchange = cellC.template getFieldPointer<FreeSurface::TEMP_MASS_EXCHANGE>();
        if (hasCellFlags(cellC,
              FreeSurface::Flags::ToFluid
            | FreeSurface::Flags::ToGas)) {
          int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);
          collected_excess += tempMassExchange[iPop_op];
        }
      }

      T mass_tmp = cell.template getField<FreeSurface::MASS>();
      mass_tmp += collected_excess;
      cell.template setField<FreeSurface::MASS>(mass_tmp);
    }
    break;
    default:
    break;
  }
}

// Setup
template<typename T, typename DESCRIPTOR>
FreeSurface3DSetup<T,DESCRIPTOR>::FreeSurface3DSetup(SuperLattice<T, DESCRIPTOR>& sLattice)
:
  sLattice{sLattice}
{}

template<typename T, typename DESCRIPTOR>
void FreeSurface3DSetup<T,DESCRIPTOR>::addPostProcessor(){
  sLattice.template addPostProcessor<FreeSurface::Stage0>(
    meta::id<FreeSurfaceMassFlowPostProcessor3D<T,DESCRIPTOR>>{});
  sLattice.template addPostProcessor<FreeSurface::Stage1>(
    meta::id<FreeSurfaceToFluidCellConversionPostProcessor3D<T,DESCRIPTOR>>{});
  sLattice.template addPostProcessor<FreeSurface::Stage2>(
    meta::id<FreeSurfaceToGasCellConversionPostProcessor3D<T,DESCRIPTOR>>{});
  sLattice.template addPostProcessor<FreeSurface::Stage3>(
    meta::id<FreeSurfaceMassExcessPostProcessor3D<T,DESCRIPTOR>>{});
  sLattice.template addPostProcessor<FreeSurface::Stage4>(
    meta::id<FreeSurfaceFinalizeConversionPostProcessor3D<T,DESCRIPTOR>>{});

  {
    // Communicate DFs, Epsilon and Cell Types
    auto& communicator = sLattice.getCommunicator(FreeSurface::Stage0());
    communicator.requestOverlap(2);
    communicator.template requestField<FreeSurface::EPSILON>();
    communicator.template requestField<FreeSurface::CELL_TYPE>();
    communicator.template requestField<descriptors::POPULATION>();
    communicator.exchangeRequests();
  }

  {
    // Communicate DFs, Cell Flags
    auto& communicator = sLattice.getCommunicator(FreeSurface::Stage1());
    communicator.requestOverlap(2);
    communicator.template requestField<FreeSurface::CELL_FLAGS>();
    communicator.template requestField<descriptors::POPULATION>();
    communicator.exchangeRequests();
  }

  {
    // Communicate Cell Flags
    auto& communicator = sLattice.getCommunicator(FreeSurface::Stage2());
    communicator.requestOverlap(2);
    communicator.template requestField<FreeSurface::CELL_FLAGS>();
    communicator.exchangeRequests();
  }

  {
    // Communicate Cell Flags
    auto& communicator = sLattice.getCommunicator(FreeSurface::Stage3());
    communicator.requestOverlap(2);
    communicator.template requestField<FreeSurface::CELL_FLAGS>();
    communicator.exchangeRequests();
  }

  {
    // Communicate TempMassExchange
    auto& communicator = sLattice.getCommunicator(FreeSurface::Stage4());
    communicator.requestOverlap(2);
    communicator.template requestField<FreeSurface::TEMP_MASS_EXCHANGE>();
    communicator.exchangeRequests();
  }

  sLattice.template addCustomTask<stage::PostStream>([&]() {
    sLattice.executePostProcessors(FreeSurface::Stage0());
    sLattice.executePostProcessors(FreeSurface::Stage1());
    sLattice.executePostProcessors(FreeSurface::Stage2());
    sLattice.executePostProcessors(FreeSurface::Stage3());
    sLattice.executePostProcessors(FreeSurface::Stage4());
  });
}
}
#endif
