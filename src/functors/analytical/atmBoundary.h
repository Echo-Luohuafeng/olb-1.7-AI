/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2024 Dennis Teutscher
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
#ifndef ATMBOUNDARY_H
#define ATMBOUNDARY_H


namespace olb {
template <typename T>
class AnalyticalWindProfileF3D : public AnalyticalF3D<T,T> {

protected:
  T _u0;
  T _z0;
  T _z_ref;
  T _kappa;
  T _d =0;
  int _windDirection;
  int _normalToEarth;


public:
  AnalyticalWindProfileF3D(T u0, T z0, T z_ref,T kappa, T d, int windDirection = 0, int normalToEarth =2) : AnalyticalF3D<T,T>( 3 )
  ,_u0(u0), _z0(z0), _z_ref(z_ref), _kappa(kappa), _d(d),_windDirection(windDirection), _normalToEarth(normalToEarth) {
    this->getName() = "AnalyticalWindProfileF3D";
  };

  bool operator()( T output[3], const T input[3] ) override {
    output[_windDirection] = (_u0*_kappa)/std::log10((_z_ref+_z0)/_z0) / _kappa *std::log10(((input[_normalToEarth]-_d)+_z0)/_z0);
    return true;
  };

};
}
#endif
