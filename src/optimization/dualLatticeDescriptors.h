/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Mathias J. Krause
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

#ifndef DUAL_LATTICE_DESCRIPTORS_H
#define DUAL_LATTICE_DESCRIPTORS_H

#include "dynamics/latticeDescriptors.h"


namespace olb {

namespace descriptors {

////////////////////////////////////////////////////////////////////////////////
// descriptor for dual lattice - 3D

using DualForcedD3Q19Descriptor = D3Q19<FORCE,F,DJDF,DJDALPHA>;

using DualForcedMRTD3Q19Descriptor = D3Q19<tag::MRT,FORCE,COORDINATE,F,DJDF,DJDALPHA>;

////////////////////////////////////////////////////////////////////////////////
// descriptor for dual lattice - 3D porous

using DualPorousD3Q19Descriptor = D3Q19<POROSITY,F,DJDF,DJDALPHA>;

} // namespace descriptors

} // namespace olb

#endif
