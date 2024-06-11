/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef CORE_CONCEPTS_H
#define CORE_CONCEPTS_H

#if (defined __cpp_concepts && (__cplusplus >= __cpp_concepts))

#include <concepts>
#include <type_traits>

#include "dynamics/descriptorBase.h"

namespace olb {

namespace concepts {

// DESCRIPTOR describes DdQq lattice with associated fields
template <typename DESCRIPTOR>
concept LatticeDescriptor = std::derived_from<DESCRIPTOR, descriptors::LATTICE_DESCRIPTOR_BASE>;

/// Basic cell exposing value-typed population references
template <typename CELL>
concept MinimalCell = requires(CELL cell) {
  // CELL::value_t exists
  typename CELL::value_t;
  // CELL::descriptor_t is a lattice descriptor
  requires LatticeDescriptor<typename CELL::descriptor_t>;
  // Population references are exposed via operator[]
  { cell[unsigned{0}] } -> std::convertible_to<const typename CELL::value_t&>;
};

/// Full cell exposing populations, moments, dynamics and associated fields
template <typename CELL>
concept Cell = requires(CELL cell, typename CELL::value_t value) {
  // Population access
  requires MinimalCell<CELL>;

  // Computation of moments
  { cell.computeRho() } -> std::convertible_to<typename CELL::value_t>;
  cell.computeU(&value);
  cell.computeJ(&value);
  cell.computeRhoU(value, &value);
  cell.computeStress(&value);
  cell.computeAllMomenta(value,&value,&value);

  // Definition of moments
  cell.defineRho(value);
  cell.defineU(&value);
  cell.defineRhoU(value,&value);
  cell.defineAllMomenta(value,&value,&value);
  cell.inverseShiftRhoU(value,&value);

  // Definition of populations
  cell.definePopulations(&value);

  // Initialize to equilibrium
  cell.iniEquilibrium(value,&value);
  cell.iniRegularized(value,&value,&value);

  // Access to dynamics
  cell.getDynamics();

  // Access to associated fields by-value
  // (using guaranteed POPULATION field as placeholder)
  { cell.template getField<descriptors::POPULATION>() } -> std::same_as<
    FieldD<typename CELL::value_t,typename CELL::descriptor_t,descriptors::POPULATION>
  >;
  cell.template setField<descriptors::POPULATION>(&value);

  // Access to associated field components by-reference
  // (using guaranteed POPULATION field as placeholder)
  { cell.template getFieldComponent<descriptors::POPULATION>(0) } -> std::convertible_to<const typename CELL::value_t&>;

  // Access to associated fields by-pointer
  // (using guaranteed POPULATION field as placeholder, return type operator[] accessible)
  { cell.template getFieldPointer<descriptors::POPULATION>()[0] } -> std::convertible_to<typename CELL::value_t>;
};

/// Cell-wise operator with scope, priority and apply template
template <typename OPERATOR>
concept CellOperator = requires(OPERATOR op) {
  // Operator scope is defined
  { OPERATOR::scope } -> std::convertible_to<OperatorScope>;
  // Operator scope is cell-wise
  requires OPERATOR::scope == OperatorScope::PerCell
        || OPERATOR::scope == OperatorScope::PerCellWithParameters;
  // Operator priority is defined
  { op.getPriority() } -> std::same_as<int>;
  // Todo: Operator provides apply method template
};

/// Block-wise operator with scope, priority and apply template
template <typename OPERATOR>
concept BlockOperator = requires(OPERATOR op) {
  // Operator scope is defined
  { OPERATOR::scope } -> std::convertible_to<OperatorScope>;
  // Operator scope is block-wise
  requires OPERATOR::scope == OperatorScope::PerBlock;
  // Operator priority is defined
  { op.getPriority() } -> std::same_as<int>;
  // Todo: Operator provides apply method template
};

}

}

#define CONCEPT(C) olb::concepts::C

#else // Fallback if concepts are not available (standard < C++20)

#define CONCEPT(C) typename

#endif

#endif
