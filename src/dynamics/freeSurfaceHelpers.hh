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

namespace olb {

namespace FreeSurface {

template<typename T, typename DESCRIPTOR>
std::enable_if_t<DESCRIPTOR::d == 2,T>
offsetHelper(T volume, const Vector<T,DESCRIPTOR::d>& sorted_normal) any_platform {
  T d2 = volume * sorted_normal[1] + 0.5 * sorted_normal[0];
  if(d2 >= sorted_normal[0]){
    return d2;
  }

  T d1 = util::sqrt(2. * sorted_normal[0] * sorted_normal[1] * volume);

  return d1;
}


// A lot of magic numbers are happening here. Optimized algorithm taken from Moritz Lehmann
template<typename T, typename DESCRIPTOR>
std::enable_if_t<DESCRIPTOR::d == 3,T>
offsetHelper(T volume, const Vector<T,DESCRIPTOR::d>& sorted_normal) any_platform {
  T sn0_plus_sn1 = sorted_normal[0] + sorted_normal[1];
  T sn0_times_sn1 = sorted_normal[0] * sorted_normal[1];
  T sn2_volume = sorted_normal[2] * volume;

  T min_sn0_plus_sn1_and_sn2 = util::min(sn0_plus_sn1, sorted_normal[2]);

  T d5 = sn2_volume + 0.5 * sn0_plus_sn1;
  if(d5 > min_sn0_plus_sn1_and_sn2 && d5 <= sorted_normal[2]){
    return d5;
  }

  T d2 = 0.5 * sorted_normal[0] + 0.28867513 * util::sqrt( util::max(0., 24. * sorted_normal[1] * sn2_volume - sorted_normal[0]*sorted_normal[0]) );

  if(d2 > sorted_normal[0] && d2 <= sorted_normal[1]){
    return d2;
  }

  T d1  = std::cbrt(6.0 * sn0_times_sn1 * sn2_volume);
  if(d1 <= sorted_normal[0]){
    return d1;
  }

  T x3 = 81.0  * sn0_times_sn1 * (sn0_plus_sn1 - 2. * sn2_volume);
  T y3 = util::sqrt(util::max(0., 23328. * sn0_times_sn1*sn0_times_sn1*sn0_times_sn1 - x3*x3 ));
  T u3 = std::cbrt(x3*x3 + y3*y3);
  T d3 = sn0_plus_sn1 - (7.5595264 * sn0_times_sn1 + 0.26456684 * u3) * (1./util::sqrt(u3)) * util::sin(0.5235988 - 0.3333334 * util::atan(y3 / x3));
  if(d3 > sorted_normal[1] && d3 <= min_sn0_plus_sn1_and_sn2){
    return d3;
  }

  T t4 = 9. * util::pow(sn0_plus_sn1 + sorted_normal[2], 2) - 18.;
  T x4 = util::max(sn0_times_sn1 * sorted_normal[2] * (324. - 648. * volume), 1.1754944e-38);
  T y4 = util::sqrt(util::max(4. * t4*t4*t4 - x4*x4, 0.));
  T u4 = std::cbrt(x4*x4 + y4*y4);
  T d4  = 0.5 * (sn0_plus_sn1 + sorted_normal[2]) - (0.20998684 * t4 + 0.13228342 * u4) * (1./util::sqrt(u4)) * util::sin(0.5235988- 0.3333334 * util::atan(y4/x4));

  return d4;
}

// A lot of magic numbers are happening here. Optimized algorithm taken from Moritz Lehmann
template<typename T, typename DESCRIPTOR>
T offsetHelperOpt(T vol, const Vector<T,DESCRIPTOR::d>& sn) any_platform {
  const T sn0_p_sn1 = sn[0] + sn[1];
  const T sn2_t_V = sn[2] * vol;

  if(sn0_p_sn1 <= 2. * sn2_t_V){
    return sn2_t_V + 0.5 * sn0_p_sn1;
  }

  const T sq_sn0 = util::pow(sn[0],2), sn1_6 = 6. * sn[1], v1 = sq_sn0 / sn1_6;

  if(v1 <= sn2_t_V && sn2_t_V < v1 + 0.5 * (sn[1]-sn[0])){
    return 0.5 *(sn[0] + util::sqrt(sq_sn0 + 8.0 * sn[1] * (sn2_t_V - v1)));
  }

  const T v6 = sn[0] * sn1_6 * sn2_t_V;
  if(sn2_t_V < v1){
    return std::cbrt(v6);
  }

  const T v3 = sn[2] < sn0_p_sn1 ? (util::pow(sn[2],2) * (3. * sn0_p_sn1 - sn[2]) + sq_sn0 *(sn[0] - 3.0 * sn[2]) + util::pow(sn[1],2)*(sn[1]-3.0 * sn[2])) / (sn[0] * sn1_6) : 0.5 * sn0_p_sn1;

  const T sq_sn0_sq_sn1 = sq_sn0 + util::pow(sn[1],2), v6_cb_sn0_sn1 = v6 - util::pow(sn[0],3) - util::pow(sn[1],3);

  const bool case34 = sn2_t_V < v3;
  const T a = case34 ? v6_cb_sn0_sn1 : 0.5 * (v6_cb_sn0_sn1 - util::pow(sn[2], 3));
  const T b = case34 ? sq_sn0_sq_sn1 : 0.5 * (sq_sn0_sq_sn1 + util::pow(sn[2], 2));
  const T c = case34 ? sn0_p_sn1 : 0.5;
  const T t = util::sqrt(util::pow(c,2) - b);
  return c - 2.0 * t * util::sin(0.33333334 * util::asin((util::pow(c,3) - 0.5 * a - 1.5 * b * c) / util::pow(t,3)));
}

template<typename T, size_t S>
std::array<T,S> solvePivotedLU(std::array<std::array<T,S>,S>& matrix, const std::array<T,S>& b, size_t N) {
  std::array<T,S> x;
  std::array<T,S> pivots;
  for(size_t i = 0; i < S; ++i){
    pivots[i] = i;
    x[i] = 0.;
  }

  N = std::min(N,S);

  for(size_t i = 0; i < N; ++i){

    T max = 0.;
    size_t max_index = i;

    for(size_t j = i; j < N; ++j){
      T abs = std::abs(matrix[pivots[j]][i]);
      if(abs > max){
        max_index = j;
        max = abs;
      }
    }

    if(max_index != i){
      size_t tmp_index = pivots[i];
      pivots[i] = pivots[max_index];
      pivots[max_index] = tmp_index;
    }

    for(size_t j = i + 1; j < N; ++j){
      matrix[pivots[j]][i] /= matrix[pivots[i]][i];

      for(size_t k = i + 1; k < N; ++k){

        matrix[pivots[j]][k] -= matrix[pivots[j]][i] * matrix[pivots[i]][k];
      }
    }
  }

  for(size_t i = 0; i  < N; ++i){
    x[i] = b[pivots[i]];

    for(size_t j = 0; j < i; ++j){
      x[i] -= matrix[pivots[i]][j] * x[j];
    }
  }

  for(size_t i = N; i > 0; --i){
    for(size_t j = i; j < N; ++j){
      x[i-1] -= matrix[pivots[i-1]][j] * x[j];
    }

    x[i-1] /= matrix[pivots[i-1]][i-1];
  }

  return x;
}

template<typename T, typename DESCRIPTOR>
void initialize(SuperLattice<T,DESCRIPTOR>& lattice) {
  lattice.executePostProcessors(FreeSurface::Stage0());
  lattice.executePostProcessors(FreeSurface::Stage1());
  lattice.executePostProcessors(FreeSurface::Stage2());
  lattice.executePostProcessors(FreeSurface::Stage3());
  lattice.executePostProcessors(FreeSurface::Stage4());
}

template <typename CELL>
bool isCellType(CELL& cell, const FreeSurface::Type& type) {
  return cell.template getField<FreeSurface::CELL_TYPE>() == type;
}

template <typename CELL>
bool hasCellFlags(CELL& cell, const FreeSurface::Flags& flags) {
  return static_cast<bool>(cell.template getField<FreeSurface::CELL_FLAGS>() & flags);
}

template <typename CELL>
bool hasNeighbour(CELL& cell, const FreeSurface::Type& type) {
  using DESCRIPTOR = typename CELL::descriptor_t;
  for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto cellC = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));
    if(isCellType(cellC, type)) {
      return true;
    }
  }

  return false;
}


template <typename CELL>
bool hasNeighbourFlags(CELL& cell, const FreeSurface::Flags& flags) {
  using DESCRIPTOR = typename CELL::descriptor_t;
  for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto cellC = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));
    if(hasCellFlags(cellC, flags)) {
      return true;
    }
  }

  return false;
}

template <typename CELL, typename V>
Vector<V,CELL::descriptor_t::d> computeInterfaceNormal(CELL& cell) {
  Vector<V,CELL::descriptor_t::d> normal{};
  normal[0] = 0.5 * (  getClampedEpsilon(cell.neighbor({-1, 0, 0}))
                     - getClampedEpsilon(cell.neighbor({ 1, 0, 0})));
  normal[1] = 0.5 * (  getClampedEpsilon(cell.neighbor({ 0,-1, 0}))
                     - getClampedEpsilon(cell.neighbor({ 0, 1, 0})));
  normal[2] = 0.5 * (  getClampedEpsilon(cell.neighbor({ 0, 0,-1}))
                     - getClampedEpsilon(cell.neighbor({ 0, 0, 1})));
  return normal;
}

template <typename CELL, typename V>
Vector<V,CELL::descriptor_t::d> computeParkerYoungInterfaceNormal(CELL& cell) {
  using DESCRIPTOR = typename CELL::descriptor_t;

  Vector<V,CELL::descriptor_t::d> normal{};

  for(size_t dim = 0; dim < CELL::descriptor_t::d; ++dim){
    normal[dim] = 0;
  }

  for(int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
      int omega_weight = 1;
      if(descriptors::c<DESCRIPTOR>(iPop, 0) != 0){
        omega_weight *= 2;
      }
      if(descriptors::c<DESCRIPTOR>(iPop, 1) != 0){
        omega_weight *= 2;
      }

      // For the 3D case
      if(CELL::descriptor_t::d == 3 && descriptors::c<DESCRIPTOR>(iPop)[2] != 0){
        omega_weight *= 2;
      }

      omega_weight /= 2;

      auto cellC = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      V epsilon = getClampedEpsilon(cellC);

      normal[0] -= omega_weight * (descriptors::c<DESCRIPTOR>(iPop, 0) * epsilon);
      normal[1] -= omega_weight * (descriptors::c<DESCRIPTOR>(iPop, 1) * epsilon);
      if(CELL::descriptor_t::d == 3){
         normal[2] -= omega_weight * (descriptors::c<DESCRIPTOR>(iPop, 2) * epsilon);
      }

  }

  return normal;
}

template <typename CELL, typename V>
V getClampedEpsilon(CELL& cell) {

  V epsilon = cell.template getField<FreeSurface::EPSILON>();
  return util::max(0., util::min(1., epsilon));
}

/*
template <typename CELL, typename V>
V getClampedEpsilonCorrected(CELL& cell) {

  if(isCellType(cell, FreeSurface::Type::Interface) && !hasNeighbour(cell, FreeSurface::Type::Gas)){
    return 1.0;
  }

  V epsilon = cell.template getField<FreeSurface::EPSILON>();

  return util::max(0., util::min(1., epsilon));
}
*/

template<typename T, typename DESCRIPTOR>
T calculateCubeOffset(T volume, const Vector<T,DESCRIPTOR::d>& normal) {
  std::vector<T> abs_normal(DESCRIPTOR::d, T{0});
  for(int i = 0; i < DESCRIPTOR::d; i++){
    abs_normal[i] = util::abs(normal[i]);
  }

  T volume_symmetry = 0.5 - util::abs(volume - 0.5);

  std::sort(abs_normal.begin(), abs_normal.end());

  if constexpr (DESCRIPTOR::d == 2) {
    abs_normal[0] = util::max(normal[0], 1e-5);
  } else if (DESCRIPTOR::d == 3){
    abs_normal[0] = util::max(normal[0], 1e-12);
    abs_normal[1] = util::max(normal[1], 1e-12);
  }

  T d = offsetHelper<T,DESCRIPTOR>(volume_symmetry, abs_normal);

  T sorted_normal_acc = 0;
  for(int i = 0; i < DESCRIPTOR::d; i++){
    sorted_normal_acc += abs_normal[i];
  }

  return std::copysign(d - 0.5 * sorted_normal_acc, volume - 0.5);
}

// Optimized version of function calculateCubeOffset for 3D
template<typename T, typename DESCRIPTOR>
T calculateCubeOffsetOpt(T volume, const Vector<T,DESCRIPTOR::d>& normal) {
  olb::Vector<T,DESCRIPTOR::d> abs_normal;

  abs_normal[0] = util::abs(normal[0]);
  abs_normal[1] = util::abs(normal[1]);
  abs_normal[2] = util::abs(normal[2]);

  T a_l1 = abs_normal[0] + abs_normal[1] + abs_normal[2];

  T volume_symmetry = 0.5 - util::abs(volume - 0.5);

  olb::Vector<T,DESCRIPTOR::d> sorted_normal;
  sorted_normal[0] = util::min(util::min(abs_normal[0], abs_normal[1]), abs_normal[2]) / a_l1;
  sorted_normal[1] = 0.;
  sorted_normal[2] = util::max(util::max(abs_normal[0], abs_normal[1]), abs_normal[2]) / a_l1;

  sorted_normal[1] = util::max(1. - sorted_normal[0] - sorted_normal[2], 0.);

  T d = offsetHelperOpt<T,DESCRIPTOR>(volume_symmetry, sorted_normal);

  return a_l1 * std::copysign(0.5 - d, volume - 0.5);
}

template <typename CELL, typename V>
V calculateSurfaceTensionCurvature(CELL& cell) {
  using DESCRIPTOR = typename CELL::descriptor_t;

  if constexpr (DESCRIPTOR::d == 2) {
    return calculateSurfaceTensionCurvature2D(cell);
  } else if (DESCRIPTOR::d == 3){
    return calculateSurfaceTensionCurvature3D(cell);
  }

  return 0;
}

template <typename CELL, typename V>
V calculateSurfaceTensionCurvature2D(CELL& cell) {
  auto normal = computeParkerYoungInterfaceNormal(cell);

  using DESCRIPTOR = typename CELL::descriptor_t;
  {
    V norm = 0.;
    for(size_t i = 0; i < DESCRIPTOR::d; ++i){
      norm += normal[i] * normal[i];
    }

    norm = util::sqrt(norm);

    if(norm < 1e-6){
      return 0.;
    }

    for(size_t i = 0; i <DESCRIPTOR::d; ++i){
      normal[i] /= norm;
    }
  }

  // Rotation matrix is
  // ( n1 | -n0 )
  // ( n0 |  n1 )

  // It is 2 because of the amount of fitting parameters. Not because of the dimension
  constexpr size_t S = 2;
  std::array<std::array<V,S>, S> lq_matrix;
  std::array<V,S> b_rhs;
  for(size_t i = 0; i < S; ++i){
    for(size_t j = 0; j < S; ++j){
      lq_matrix[i][j] = 0.;
    }
    b_rhs[i] = 0.;
  }

  // Offset for the plic correction
  V origin_offset = 0.;
  {
    V fill_level = getClampedEpsilon(cell);
    origin_offset = calculateCubeOffset<V,DESCRIPTOR>(fill_level, normal);
  }

  // The amount of neighbouring interfaces. if less are available we will solve a reduced curve by setting the less important parameters to zero
  std::size_t healthy_interfaces = 0;

  for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
    auto cellC = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

    if(   !isCellType(cellC, FreeSurface::Type::Interface)
       || !hasNeighbour(cellC, FreeSurface::Type::Gas)) {
      continue;
    }

    ++healthy_interfaces;

    V fill_level = getClampedEpsilon(cellC);

    V cube_offset = calculateCubeOffset<V,DESCRIPTOR>(fill_level, normal);

    V x_pos = descriptors::c<DESCRIPTOR>(iPop,0);
    V y_pos = descriptors::c<DESCRIPTOR>(iPop,1);

    // Rotation
    V rot_x_pos = x_pos * normal[1] - y_pos * normal[0];
    V rot_y_pos = x_pos * normal[0] + y_pos * normal[1] + (cube_offset - origin_offset);

    V rot_x_pos_2 = rot_x_pos * rot_x_pos;
    V rot_x_pos_3 = rot_x_pos_2 * rot_x_pos;
    V rot_x_pos_4 = rot_x_pos_3 * rot_x_pos;

    lq_matrix[1][1] += rot_x_pos_2;
    lq_matrix[1][0] += rot_x_pos_3;
    lq_matrix[0][0] += rot_x_pos_4;

    b_rhs[0] += rot_x_pos_2*(rot_y_pos);
    b_rhs[1] += rot_x_pos*(rot_y_pos);
  }

  lq_matrix[0][1] = lq_matrix[1][0];

  // Thikonov regularization parameter
  V alpha = 0.0;
  for(size_t i = 0; i < DESCRIPTOR::d; ++i){
    lq_matrix[i][i] += alpha;
  }

  // It is 2 because of the fitting parameters. Not dependent on the dimension
  std::array<V,S> solved_fit = FreeSurface::solvePivotedLU<V,S>(lq_matrix, b_rhs, healthy_interfaces);

  // signed curvature -> kappa = y'' / ( (1 + y'²)^(3/2) )
  V denom = std::sqrt(1. + solved_fit[1]*solved_fit[1]);
  denom = denom * denom * denom;
  V curvature = 2.*solved_fit[0] / denom;
  return util::max(-1., util::min(1., curvature));
}

template <typename CELL, typename V>
V calculateSurfaceTensionCurvature3D(CELL& cell){
  // This is b_z
  auto normal = computeParkerYoungInterfaceNormal(cell);

  using DESCRIPTOR = typename CELL::descriptor_t;
  {
    V norm = 0.;
    for(size_t i = 0; i < DESCRIPTOR::d; ++i){
      norm += normal[i] * normal[i];
    }

    norm = util::sqrt(norm);

    if(norm < 1e-12){
      return 0.;
    }

    for(size_t i = 0; i <DESCRIPTOR::d; ++i){
      normal[i] /= norm;
    }
  }

  std::array<V,3> r_vec{
    0.56270900, 0.32704452, 0.75921047
  };
  /*
  std::array<T,DESCRIPTOR::d> r_vec{
    0.,0.,1.
  };
  */
  std::array<std::array<V,3>,3> rotation{{
    {{0., 0., 0.}},
    //{{normal[1], -normal[0], 0.}},
    {{normal[1] * r_vec[2] - normal[2] * r_vec[1], normal[2] * r_vec[0] - normal[0] * r_vec[2], normal[0] * r_vec[1] - normal[1] * r_vec[0]}},
    {{normal[0], normal[1], normal[2]}}
  }};

  // Cross product with (0,0,1) x normal
  // This is b_y

  // (normal[0], normal[1], normal[2])

  V cross_norm = 0.;
  for(size_t i = 0; i < DESCRIPTOR::d; ++i){
    cross_norm += rotation[1][i] * rotation[1][i];
  }

  // If too close too each other use the identity matrix
  if(cross_norm > 1e-6){

    cross_norm = util::sqrt(cross_norm);

    for(size_t i = 0; i <DESCRIPTOR::d; ++i){
      rotation[1][i] /= cross_norm;
    }
  }else {

    rotation[1] = {{
      -normal[2],
      0.,
      normal[0]
    }};

    cross_norm = 0.;
    for(size_t i = 0; i < DESCRIPTOR::d; ++i){
      cross_norm += rotation[1][i] * rotation[1][i];
    }

    cross_norm = util::sqrt(cross_norm);

    for(size_t i = 0; i <DESCRIPTOR::d; ++i){
      rotation[1][i] /= cross_norm;
    }
  }

  // Cross product of ((0,0,1) x normal / | (0,0,1) x normal |) x normal
  // This is b_x
  rotation[0] = {{
    rotation[1][1] * normal[2] - rotation[1][2] * normal[1],
    rotation[1][2] * normal[0] - rotation[1][0] * normal[2],
    rotation[1][0] * normal[1] - rotation[1][1] * normal[0]
  }};

  // These three form a matrix and are entered into each row
  // ( b_x )
  // ( b_y )
  // ( b_z )

  constexpr size_t S = 5;
  std::array<std::array<V,S>, S> lq_matrix;
  std::array<V,S> b_rhs;
  for(size_t i = 0; i < S; ++i){
    for(size_t j = 0; j < S; ++j){
      lq_matrix[i][j] = 0.;
    }
    b_rhs[i] = 0.;
  }
  V origin_offset = 0.;
  {
    V fill_level = getClampedEpsilon(cell);
    origin_offset = calculateCubeOffsetOpt<V,DESCRIPTOR>(fill_level, normal);
  }

  size_t healthy_interfaces = 0;
  for(int iPop = 1; iPop < DESCRIPTOR::q; iPop++){
    auto cellC = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

    if(!isCellType(cellC, FreeSurface::Type::Interface) || !hasNeighbour(cellC, FreeSurface::Type::Gas)){
      continue;
    }

    ++healthy_interfaces;

    V fill_level = getClampedEpsilon(cellC);

    V cube_offset = calculateCubeOffsetOpt<V, DESCRIPTOR>(fill_level, normal);

    int i = descriptors::c<DESCRIPTOR>(iPop)[0];
    int j = descriptors::c<DESCRIPTOR>(iPop)[1];
    int k = descriptors::c<DESCRIPTOR>(iPop)[2];

    std::array<V,3> pos{static_cast<V>(i),static_cast<V>(j),static_cast<V>(k)};
    std::array<V,3> r_pos{0.,0.,cube_offset - origin_offset};

    for(size_t a = 0; a < DESCRIPTOR::d; ++a){
      for(size_t b = 0; b < DESCRIPTOR::d; ++b){
        r_pos[a] += rotation[a][b] * pos[b];
      }
    }

    V r_x_2 = r_pos[0] * r_pos[0];
    V r_x_3 = r_x_2 * r_pos[0];
    V r_x_4 = r_x_3 * r_pos[0];

    V r_y_2 = r_pos[1] * r_pos[1];
    V r_y_3 = r_y_2 * r_pos[1];
    V r_y_4 = r_y_3 * r_pos[1];

    V r_x_2_y_2 = r_x_2 * r_y_2;
    V r_x_3_y = r_x_3 * r_pos[1];
    V r_x_2_y = r_x_2 * r_pos[1];

    V r_x_y_3 = r_pos[0] * r_y_3;
    V r_x_y_2 = r_pos[0] * r_y_2;

    V r_x_y = r_pos[0] * r_pos[1];

    lq_matrix[0][0] += r_x_4;
    lq_matrix[1][1] += r_y_4;
    lq_matrix[2][2] += r_x_2_y_2;
    lq_matrix[3][3] += r_x_2;
    lq_matrix[4][4] += r_y_2;

    // skip [1][0] copy later from [2][2]
    lq_matrix[2][0] += r_x_3_y;
    lq_matrix[3][0] += r_x_3;
    lq_matrix[4][0] += r_x_2_y;

    lq_matrix[2][1] += r_x_y_3;
    lq_matrix[3][1] += r_x_y_2;
    lq_matrix[4][1] += r_y_3;

    // skip [3][2] copy from [4][0]
    // skip [4][2] copy from [3][1]

    lq_matrix[4][3] += r_x_y;

    b_rhs[0] +=  r_x_2 * r_pos[2];
    b_rhs[1] +=  r_y_2 * r_pos[2];
    b_rhs[2] +=  r_x_y * r_pos[2];
    b_rhs[3] +=  r_pos[0] * r_pos[2];
    b_rhs[4] +=  r_pos[1] * r_pos[2];
  }

  lq_matrix[1][0] = lq_matrix[2][2];
  lq_matrix[3][2] = lq_matrix[4][0];
  lq_matrix[4][2] = lq_matrix[3][1];

  for(size_t i = 0; i < S; ++i){
    for(size_t j = i + 1; j < S; ++j){
      lq_matrix[i][j] = lq_matrix[j][i];
    }
  }

  // Consider using Thikonov regularization?
  //T alpha = 1e-8;
  V alpha = 0.0;
  for(size_t i = 0; i < S; ++i){
    lq_matrix[i][i] += alpha;
  }

  std::array<V,S> solved_fit = FreeSurface::solvePivotedLU<V,S>(lq_matrix, b_rhs, healthy_interfaces);

  V denom = std::sqrt(1. + solved_fit[3]*solved_fit[3] + solved_fit[4]*solved_fit[4]);
  denom = denom * denom * denom;
  V curvature = ( (1.+solved_fit[4]*solved_fit[4]) * solved_fit[0] + (1. + solved_fit[3]*solved_fit[3] ) * solved_fit[1] - solved_fit[3] * solved_fit[4] * solved_fit[2] ) / denom;

  return util::max(-1., util::min(1., curvature));
}

template<typename T, typename DESCRIPTOR>
T plicInverse(T d_o, const Vector<T,DESCRIPTOR::d>& normal){

  Vector<T,DESCRIPTOR::d> abs_normal;
  for(int i = 0; i < DESCRIPTOR::d; i++){
    abs_normal[i] = util::abs(normal[i]);
  }

  //const T n1 = std::min_element(abs_normal.begin(), abs_normal.end());
  //const T n2 = std::max_element(abs_normal.begin(), abs_normal.end());

  const T n1 = abs_normal[0];
  const T n2 = abs_normal[1];

  for(int i = 0; i < DESCRIPTOR::d; i++){
    if(abs_normal[i] < n1){
      n1 = abs_normal[i];
    }

    if(abs_normal[i] > n2){
      n2 = abs_normal[i];
    }
  }

  T abs_normal_acc = 0;
  for(int i = 0; i < DESCRIPTOR::d; i++){
    abs_normal_acc += abs_normal[i];
  }
  const T n3 = abs_normal_acc - n1 - n2;
  const T d = 0.5 * (n1+n2+n3) - util::abs(d_o);
  T vol;

  if(DESCRIPTOR::d == 2){
    if(d < n1){
      vol = d * d / (2. * n1 * n2);
    } else if(d >= n1){
      vol = d / n2 - n1 / (2. * n2);
    }
  }

  if(DESCRIPTOR::d == 3){
    if(util::min(n1+n3,n2) <= d && d <= n2){
      vol = (d-0.5 *(n1+n3))/n2;
    } else if(d < n1){
      vol = util::pow(d,3) / (6. * n1 * n2 * n3);
    } else if(d <= n3){
      vol = (3.0 * d * (d-n1) + std::pow(n1,2))/(6. * n2 * n3);
    } else {
      vol = (util::pow(d,3) - util::pow(d-n1,3) - util::pow(d-n3,3) - util::pow(util::max(0., d-n2),3)) / (6. * n1* n2 * n3);
    }
  }

  return std::copysign(0.5 - vol, d_o) + 0.5;
}

template <typename CELL>
NeighbourInfo getNeighbourInfo(CELL& cell) {
  NeighbourInfo info{};
  using DESCRIPTOR = typename CELL::descriptor_t;

  for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
    auto cellC = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

    if(isCellType(cellC, FreeSurface::Type::Gas)){
      info.has_gas_neighbours = true;
    }
    else if(isCellType(cellC, FreeSurface::Type::Fluid)){
      info.has_fluid_neighbours = true;
    }
    else if(isCellType(cellC, FreeSurface::Type::Interface)){
      ++info.interface_neighbours;
    }
  }
  return info;
}

template <typename CELL, typename V>
bool isHealthyInterface(CELL& cell) {
  bool has_fluid_neighbours = false;
  bool has_gas_neighbours = false;

  if(!isCellType(cell, FreeSurface::Type::Interface)){
    return false;
  }

  using DESCRIPTOR = typename CELL::descriptor_t;
  for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
    auto cellC = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));
    if(isCellType(cellC, FreeSurface::Type::Gas)){
      has_gas_neighbours = true;
      if(has_fluid_neighbours){
        return true;
      }
    }
    else if(isCellType(cellC, FreeSurface::Type::Fluid)){
      has_fluid_neighbours = true;
      if(has_gas_neighbours){
        return true;
      }
    }
  }
  return false;
}

template <typename CELL, typename V>
void setCellType(CELL& cell, const FreeSurface::Type& type) {
  cell.template setField<FreeSurface::CELL_TYPE>(type);
}

template <typename CELL, typename V>
void setCellFlags(CELL& cell, const FreeSurface::Flags& flags){
  cell.template setField<FreeSurface::CELL_FLAGS>(flags);
}

} // namespace FreeSurface

} // namespace olb
