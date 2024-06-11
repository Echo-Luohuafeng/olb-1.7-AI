/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef DYNAMICS_GUO_ZHAO_DYNAMICS_CSE_H
#define DYNAMICS_GUO_ZHAO_DYNAMICS_CSE_H


#ifndef DISABLE_CSE

#include "equilibrium.h"
#include "collision.h"
#include "latticeDescriptors.h"

namespace olb {

namespace collision {

namespace detail {

template <typename... FIELDS>
struct SmagorinskyEffectiveOmega<BGK,descriptors::D3Q19<FIELDS...>,momenta::Forced<momenta::BulkTuple>,guoZhao::GuoZhaoSecondOrder> {

template <CONCEPT(MinimalCell) CELL, typename PARAMETERS, typename V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform
{
auto x21 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x20 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x28 = parameters.template get<olb::descriptors::OMEGA>();
auto x29 = parameters.template get<olb::collision::LES::Smagorinsky>();
auto x22 = cell.template getFieldComponent<olb::descriptors::EPSILON>(0);
auto x19 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x23 = cell[14] + cell[5];
auto x24 = cell[18] + cell[9];
auto x25 = x23 + x24 + cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[15] + cell[16] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[7] + cell[8] + V{1};
auto x26 = V{1} / (x25);
auto x27 = V{0.333333333333333}*cell[4];
auto x30 = V{0.333333333333333}*cell[5];
auto x31 = V{0.333333333333333}*cell[13];
auto x32 = V{0.333333333333333}*cell[14];
auto x33 = V{0.5}*x21;
auto x34 = V{1}*cell[12];
auto x35 = V{1}*cell[15];
auto x36 = V{1}*cell[17];
auto x37 = V{1}*cell[3];
auto x38 = V{1}*cell[16];
auto x39 = V{1}*cell[18];
auto x40 = V{1}*cell[9];
auto x41 = V{1}*cell[8];
auto x42 = -x41;
auto x43 = x40 + x42;
auto x44 = V{1}*cell[7];
auto x45 = V{1}*cell[6];
auto x46 = -x45;
auto x47 = x44 + x46;
auto x48 = x26*(x34 + x35 + x36 - x37 - x38 - x39 + x43 + x47);
auto x49 = x33 + x48;
auto x50 = x49*x49;
auto x51 = V{0.333333333333333}*cell[0];
auto x52 = V{0.333333333333333}*cell[1];
auto x53 = V{0.333333333333333}*cell[10];
auto x54 = x51 + x52 + x53 - V{0.666666666666667}*cell[17] - V{0.666666666666667}*cell[18] - V{0.666666666666667}*cell[8] - V{0.666666666666667}*cell[9];
auto x55 = V{0.333333333333333}*cell[2];
auto x56 = V{0.333333333333333}*cell[11];
auto x57 = x55 + x56 - V{0.666666666666667}*cell[15] - V{0.666666666666667}*cell[16] - V{0.666666666666667}*cell[6] - V{0.666666666666667}*cell[7];
auto x58 = x25*x50 + x27 + x30 + x31 + x32 + x54 + x57 - V{0.666666666666667}*cell[12] - V{0.666666666666667}*cell[3];
auto x59 = V{0.333333333333333}*cell[6];
auto x60 = V{0.333333333333333}*cell[7];
auto x61 = V{0.333333333333333}*cell[15];
auto x62 = V{0.333333333333333}*cell[16];
auto x63 = V{0.5}*x20;
auto x64 = V{1}*cell[11];
auto x65 = V{1}*cell[13];
auto x66 = V{1}*cell[2];
auto x67 = V{1}*cell[14];
auto x68 = V{1}*cell[5];
auto x69 = V{1}*cell[4];
auto x70 = -x69;
auto x71 = x68 + x70;
auto x72 = x26*(x36 + x39 - x40 + x42 + x64 + x65 - x66 - x67 + x71);
auto x73 = x63 + x72;
auto x74 = x73*x73;
auto x75 = V{0.333333333333333}*cell[3];
auto x76 = V{0.333333333333333}*cell[12];
auto x77 = x75 + x76 - V{0.666666666666667}*cell[13] - V{0.666666666666667}*cell[14] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[5];
auto x78 = x25*x74 + x54 + x59 + x60 + x61 + x62 + x77 - V{0.666666666666667}*cell[11] - V{0.666666666666667}*cell[2];
auto x79 = V{0.333333333333333}*cell[8];
auto x80 = V{0.333333333333333}*cell[9];
auto x81 = V{0.333333333333333}*cell[17];
auto x82 = V{0.333333333333333}*cell[18];
auto x83 = V{1}*cell[10];
auto x84 = V{1}*cell[1];
auto x85 = V{0.5}*x19 + x26*(x35 + x38 - x44 + x46 + x65 + x67 - x68 + x70 + x83 - x84);
auto x86 = x85*x85;
auto x87 = x25*x86 + x51 + x57 + x77 + x79 + x80 + x81 + x82 - V{0.666666666666667}*cell[10] - V{0.666666666666667}*cell[1];
auto x88 = x25*x73;
auto x89 = x85*x88;
auto x90 = x49*x88;
auto x91 = x25*x49*x85 - x35 + x38 + x47;
auto x92 = V{1} / (V{3.00000046417339}*util::sqrt(x26*(x29*x29)*util::sqrt((x23 + x89 - cell[13] - cell[4])*(-x65 + x67 + x71 + x89) + (x24 + x90 - cell[17] - cell[8])*(-x36 + x39 + x43 + x90) + V{0.5}*(x58*x58) + V{0.5}*(x78*x78) + V{0.5}*(x87*x87) + x91*x91) + V{0.0277777691819762}/((x28)*(x28))) + V{0.5}/x28);
auto x93 = V{1} / (x22);
auto x94 = x93*(V{1.5}*x50 + V{1.5}*x74 + V{1.5}*x86);
auto x95 = V{1} - x94;
auto x96 = V{1} - x92;
auto x97 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x98 = V{4.5}*cell[14];
auto x99 = V{4.5}*cell[16];
auto x100 = V{4.5}*cell[5];
auto x101 = V{4.5}*cell[7];
auto x102 = V{4.5}*cell[13] - V{4.5}*cell[4];
auto x103 = V{4.5}*cell[15] - V{4.5}*cell[6];
auto x104 = V{2.25}*x19 + x26*(-x100 - x101 + x102 + x103 + x98 + x99 + V{4.5}*cell[10] - V{4.5}*cell[1]);
auto x105 = x104*x85*x93;
auto x106 = x94 + V{-1};
auto x107 = V{1.5}*x19;
auto x108 = V{3}*cell[14];
auto x109 = V{3}*cell[16];
auto x110 = V{3}*cell[5];
auto x111 = V{3}*cell[7];
auto x112 = V{3}*cell[13] - V{3}*cell[4];
auto x113 = V{3}*cell[15] - V{3}*cell[6];
auto x114 = x26*(x108 + x109 - x110 - x111 + x112 + x113 + V{3}*cell[10] - V{3}*cell[1]);
auto x115 = x107 + x114;
auto x116 = x106 + x115;
auto x117 = V{2.25}*x20;
auto x118 = V{4.5}*cell[18];
auto x119 = V{4.5}*cell[9];
auto x120 = V{4.5}*cell[17] - V{4.5}*cell[8];
auto x121 = x26*(x100 + x102 + x118 - x119 + x120 - x98 + V{4.5}*cell[11] - V{4.5}*cell[2]);
auto x122 = x117 + x121;
auto x123 = x122*x73*x93;
auto x124 = V{1.5}*x20;
auto x125 = V{3}*cell[18];
auto x126 = V{3}*cell[9];
auto x127 = V{3}*cell[17] - V{3}*cell[8];
auto x128 = x26*(-x108 + x110 + x112 + x125 - x126 + x127 + V{3}*cell[11] - V{3}*cell[2]);
auto x129 = x124 + x128;
auto x130 = x106 + x129;
auto x131 = V{2.25}*x21;
auto x132 = x26*(x101 + x103 - x118 + x119 + x120 - x99 + V{4.5}*cell[12] - V{4.5}*cell[3]);
auto x133 = x131 + x132;
auto x134 = x133*x49*x93;
auto x135 = V{1.5}*x21;
auto x136 = x26*(-x109 + x111 + x113 - x125 + x126 + x127 + V{3}*cell[12] - V{3}*cell[3]);
auto x137 = x135 + x136;
auto x138 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x139 = x93*(x104 + x122)*(x73 + x85);
auto x140 = x93*(x104 - x117 - x121)*(-x63 - x72 + x85);
auto x141 = x129 + x95;
auto x142 = -x107 - x114;
auto x143 = x93*(x104 + x133)*(x49 + x85);
auto x144 = -x33 - x48;
auto x145 = -x131 - x132;
auto x146 = x93*(x104 + x145)*(x144 + x85);
auto x147 = x137 + x95;
auto x148 = x93*(x122 + x133)*(x49 + x73);
auto x149 = x93*(x122 + x145)*(x144 + x73);
auto x150 = -x124 - x128;
auto x151 = x115 + x95;
auto x152 = -x135 - x136;
cell[0] = x92*(x95*(x27 + x30 + x31 + x32 + x51 + x52 + x53 + x55 + x56 + x59 + x60 + x61 + x62 + x75 + x76 + x79 + x80 + x81 + x82 + V{0.333333333333333}) + V{-0.333333333333333}) + V{1}*x96*cell[0];
cell[1] = x84*x96 - x92*(x97*(-x105 + x116) + V{0.0555555555555556});
cell[2] = x66*x96 - x92*(x97*(-x123 + x130) + V{0.0555555555555556});
cell[3] = x37*x96 - x92*(x97*(x106 - x134 + x137) + V{0.0555555555555556});
cell[4] = x69*x96 - x92*(x138*(x116 + x129 - x139) + V{0.0277777777777778});
cell[5] = x68*x96 + x92*(x138*(x140 + x141 + x142) + V{-0.0277777777777778});
cell[6] = x45*x96 - x92*(x138*(x116 + x137 - x143) + V{0.0277777777777778});
cell[7] = x44*x96 + x92*(x138*(x142 + x146 + x147) + V{-0.0277777777777778});
cell[8] = x41*x96 - x92*(x138*(x130 + x137 - x148) + V{0.0277777777777778});
cell[9] = x40*x96 + x92*(x138*(x147 + x149 + x150) + V{-0.0277777777777778});
cell[10] = x83*x96 + x92*(x97*(x105 + x151) + V{-0.0555555555555556});
cell[11] = x64*x96 + x92*(x97*(x123 + x141) + V{-0.0555555555555556});
cell[12] = x34*x96 + x92*(x97*(x134 + x147) + V{-0.0555555555555556});
cell[13] = x65*x96 + x92*(x138*(x129 + x139 + x151) + V{-0.0277777777777778});
cell[14] = x67*x96 + x92*(x138*(x140 + x150 + x151) + V{-0.0277777777777778});
cell[15] = x35*x96 + x92*(x138*(x137 + x143 + x151) + V{-0.0277777777777778});
cell[16] = x38*x96 + x92*(x138*(x146 + x151 + x152) + V{-0.0277777777777778});
cell[17] = x36*x96 + x92*(x138*(x137 + x141 + x148) + V{-0.0277777777777778});
cell[18] = x39*x96 + x92*(x138*(x141 + x149 + x152) + V{-0.0277777777777778});
return { x25, x50 + x74 + x86 };
}

template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
auto computeEffectiveOmega(CELL& cell, PARAMETERS& parameters) any_platform
{
auto x19 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x21 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x28 = parameters.template get<olb::descriptors::OMEGA>();
auto x20 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x29 = parameters.template get<olb::collision::LES::Smagorinsky>();
auto x0 = cell[14] + cell[5];
auto x1 = cell[18] + cell[9];
auto x2 = x0 + x1 + cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[15] + cell[16] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[7] + cell[8] + V{1};
auto x3 = V{1} / (x2);
auto x4 = V{1}*cell[15];
auto x5 = V{1}*cell[17];
auto x6 = V{1}*cell[16];
auto x7 = V{1}*cell[18];
auto x8 = V{1}*cell[9];
auto x9 = -V{1}*cell[8];
auto x10 = x8 + x9;
auto x11 = V{1}*cell[7];
auto x12 = -V{1}*cell[6];
auto x13 = x11 + x12;
auto x14 = V{0.5}*x21 + x3*(x10 + x13 + x4 + x5 - x6 - x7 + V{1}*cell[12] - V{1}*cell[3]);
auto x15 = V{0.333333333333333}*cell[0];
auto x16 = x15 + V{0.333333333333333}*cell[10] - V{0.666666666666667}*cell[17] - V{0.666666666666667}*cell[18] + V{0.333333333333333}*cell[1] - V{0.666666666666667}*cell[8] - V{0.666666666666667}*cell[9];
auto x17 = V{0.333333333333333}*cell[11] - V{0.666666666666667}*cell[15] - V{0.666666666666667}*cell[16] + V{0.333333333333333}*cell[2] - V{0.666666666666667}*cell[6] - V{0.666666666666667}*cell[7];
auto x18 = x16 + x17 + x2*(x14*x14) - V{0.666666666666667}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5];
auto x22 = V{1}*cell[13];
auto x23 = V{1}*cell[14];
auto x24 = V{1}*cell[5];
auto x25 = -V{1}*cell[4];
auto x26 = x24 + x25;
auto x27 = V{0.5}*x20 + x3*(x22 - x23 + x26 + x5 + x7 - x8 + x9 + V{1}*cell[11] - V{1}*cell[2]);
auto x30 = V{0.333333333333333}*cell[12] - V{0.666666666666667}*cell[13] - V{0.666666666666667}*cell[14] + V{0.333333333333333}*cell[3] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[5];
auto x31 = x16 + x2*(x27*x27) + x30 - V{0.666666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7];
auto x32 = V{0.5}*x19 + x3*(-x11 + x12 + x22 + x23 - x24 + x25 + x4 + x6 + V{1}*cell[10] - V{1}*cell[1]);
auto x33 = x15 + x17 + x2*(x32*x32) + x30 - V{0.666666666666667}*cell[10] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9];
auto x34 = x2*x27;
auto x35 = x32*x34;
auto x36 = x14*x34;
auto x37 = x13 + x14*x2*x32 - x4 + x6;
return V{1} / (V{3.00000046417339}*util::sqrt(x3*(x29*x29)*util::sqrt((x0 + x35 - cell[13] - cell[4])*(-x22 + x23 + x26 + x35) + (x1 + x36 - cell[17] - cell[8])*(x10 + x36 - x5 + x7) + V{0.5}*(x18*x18) + V{0.5}*(x31*x31) + V{0.5}*(x33*x33) + x37*x37) + V{0.0277777691819762}/((x28)*(x28))) + V{0.5}/x28);
}

};


}

}

}

#endif

#endif