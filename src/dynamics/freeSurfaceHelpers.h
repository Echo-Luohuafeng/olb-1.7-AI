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

#pragma once

#include "dynamics/descriptorField.h"


#include <cstdint>
#include <array>

namespace olb {

namespace FreeSurface {

struct Stage0 { };
struct Stage1 { };
struct Stage2 { };
struct Stage3 { };
struct Stage4 { };

enum class Type : std::uint8_t {
  Gas = 0,
  Interface = 1,
  Fluid = 2,
  Solid = 4
};

enum class Flags : std::uint8_t {
  ToGas = 1,
  ToFluid = 2,
  NewInterface = 4
};

struct NeighbourInfo {
  bool has_fluid_neighbours = false;
  bool has_gas_neighbours = false;
  size_t interface_neighbours = 0;
};

struct CELL_TYPE : public descriptors::TYPED_FIELD_BASE<Type,1> { };
struct CELL_FLAGS : public descriptors::TYPED_FIELD_BASE<Flags,1> { };
struct MASS : public descriptors::FIELD_BASE<1> { };
struct EPSILON : public descriptors::FIELD_BASE<1> { };
struct PREVIOUS_VELOCITY : public descriptors::FIELD_BASE<0,1,0>{};
struct TEMP_MASS_EXCHANGE : public descriptors::FIELD_BASE<0,0,1>{};

// Variables for PostProcessor
struct DROP_ISOLATED_CELLS        : public descriptors::FIELD_BASE<1> { };
struct TRANSITION                 : public descriptors::FIELD_BASE<1> { };
struct LONELY_THRESHOLD           : public descriptors::FIELD_BASE<1> { };
struct HAS_SURFACE_TENSION        : public descriptors::FIELD_BASE<1> { };
struct SURFACE_TENSION_PARAMETER  : public descriptors::FIELD_BASE<1> { };
struct FORCE_CONVERSION_FACTOR    : public descriptors::FIELD_BASE<1> { };
struct LATTICE_SIZE               : public descriptors::FIELD_BASE<1> { };

template<typename T, size_t S>
static std::array<T,S> solvePivotedLU(std::array<std::array<T,S>,S>& matrix, const std::array<T,S>& b, size_t N = S) any_platform;

template<typename T, typename DESCRIPTOR>
void initialize(SuperLattice<T,DESCRIPTOR>& lattice);

template <typename CELL>
static bool isCellType(CELL& cell, const FreeSurface::Type& type) any_platform;

template <typename CELL>
static bool hasCellFlags(CELL& cell, const FreeSurface::Flags& type) any_platform;

template <typename CELL>
static bool hasNeighbour(CELL& cell, const FreeSurface::Type& type) any_platform;

template <typename CELL>
static bool hasNeighbourFlags(CELL& cell, const FreeSurface::Flags& flags) any_platform;

template <typename CELL, typename V=typename CELL::value_t>
static Vector<V,CELL::descriptor_t::d> computeInterfaceNormal(CELL& cell) any_platform;

template <typename CELL, typename V=typename CELL::value_t>
static Vector<V,CELL::descriptor_t::d> computeParkerYoungInterfaceNormal(CELL& cell) any_platform;

template <typename CELL, typename V=typename CELL::value_t>
static V getClampedEpsilon(CELL& cell) any_platform;

/*
template <typename CELL, typename V=typename CELL::value_t>
static V getClampedEpsilonCorrected(CELL& cell) any_platform;
*/

template<typename T, typename DESCRIPTOR>
static T calculateCubeOffset(T volume, const Vector<T,DESCRIPTOR::d>& normal) any_platform;

template<typename T, typename DESCRIPTOR>
static T calculateCubeOffsetOpt(T volume, const Vector<T,DESCRIPTOR::d>& normal) any_platform;

template <typename CELL, typename V=typename CELL::value_t>
static V calculateSurfaceTensionCurvature(CELL& cell) any_platform;

template <typename CELL, typename V=typename CELL::value_t>
static V calculateSurfaceTensionCurvature2D(CELL& cell) any_platform;

template <typename CELL, typename V=typename CELL::value_t>
static V calculateSurfaceTensionCurvature3D(CELL& cell) any_platform;

template<typename T, typename DESCRIPTOR>
static T plicInverse(T d, const Vector<T,DESCRIPTOR::d>& normal) any_platform;

template <typename CELL>
static NeighbourInfo getNeighbourInfo(CELL& cell) any_platform;

template <typename CELL, typename V=typename CELL::value_t>
static bool isHealthyInterface(CELL& cell) any_platform;

template <typename CELL, typename V=typename CELL::value_t>
static void setCellType(CELL& cell, const FreeSurface::Type& type) any_platform;

template <typename CELL, typename V=typename CELL::value_t>
static void setCellFlags(CELL& cell, const FreeSurface::Flags& flags) any_platform;

} // namespace FreeSurface

inline any_platform FreeSurface::Flags operator&(FreeSurface::Flags lhs, FreeSurface::Flags rhs) {
  return static_cast<FreeSurface::Flags>(static_cast<std::uint8_t>(lhs) & static_cast<std::uint8_t>(rhs));
}

inline any_platform FreeSurface::Flags operator|(FreeSurface::Flags lhs, FreeSurface::Flags rhs) {
  return static_cast<FreeSurface::Flags>(static_cast<std::uint8_t>(lhs) | static_cast<std::uint8_t>(rhs));
}

} // namespace olb

#include "freeSurfaceHelpers.hh"