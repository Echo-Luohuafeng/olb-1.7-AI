/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Fedor Bukreev, Adrian Kummerländer
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

//This file contains the AdvectionDiffusionZeroGradientBoundary
//This is a new version of the Boundary, which only contains free floating functions

#ifndef SET_NEW_SLIP_BOUNDARY_3D_H
#define SET_NEW_SLIP_BOUNDARY_3D_H

namespace olb {

template<typename T, typename DESCRIPTOR>
void setNewSlipBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,
                        SuperGeometry<T,3>& superGeometry,
                        int material);

template<typename T, typename DESCRIPTOR>
void setNewSlipBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,
                        FunctorPtr<SuperIndicatorF3D<T>>&& indicator);

template<typename T, typename DESCRIPTOR>
void setNewSlipBoundary(BlockLattice<T,DESCRIPTOR>& _block,
                        BlockIndicatorF3D<T>& indicator,
                        bool includeOuterCells);

}

#endif
