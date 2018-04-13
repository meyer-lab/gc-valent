#include <algorithm>
#include <cstdio>
#include <numeric>
#include <array>
#include <vector>
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode.h>            /* prototypes for CVODE fcts., consts. */
#include <string>
#include <sundials/sundials_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <cvode/cvode_direct.h>
#include <iostream>
#include "model.hpp"

using std::array;
using std::copy;
using std::vector;
using std::fill;
using std::string;


ratesS param(const double * const rxn, const double * const tfR) {
	ratesS r;

	r.IL2 = rxn[0];
	r.IL15 = rxn[1];
	r.IL7 = rxn[2];
	r.IL9 = rxn[3];
	r.kfwd = rxn[4];
	r.k5rev = rxn[5];
	r.k6rev = rxn[6];
	r.k15rev = rxn[7];
	r.k17rev = rxn[8];
	r.k18rev = rxn[9];
	r.k22rev = rxn[10];
	r.k23rev = rxn[11];
	r.k27rev = rxn[12];
	r.k29rev = rxn[13];
	r.k31rev = rxn[14];
	
	// These are probably measured in the literature
	r.k10rev = 12.0 * r.k5rev / 1.5; // doi:10.1016/j.jmb.2004.04.038
	r.k11rev = 63.0 * r.k5rev / 1.5; // doi:10.1016/j.jmb.2004.04.038
	// To satisfy detailed balance these relationships should hold
	// _Based on initial assembly steps
	r.k4rev = kfbnd * r.k6rev * k3rev / k1rev / k3fwd;
	r.k7rev = k3fwd * k2rev * r.k5rev / kfbnd / k3rev;
	r.k12rev = k1rev * r.k11rev / k2rev;
	// _Based on formation of full complex
	r.k9rev = k2rev * r.k10rev * r.k12rev / kfbnd / k3rev / r.k6rev * k3fwd;
	r.k8rev = k2rev * r.k10rev * r.k12rev / kfbnd / r.k7rev / k3rev * k3fwd;

	// IL15
	// To satisfy detailed balance these relationships should hold
	// _Based on initial assembly steps
	r.k16rev = r.kfwd * r.k18rev * r.k15rev / k13rev / kfbnd;
	r.k19rev = r.kfwd * k14rev * r.k17rev / kfbnd / r.k15rev;
	r.k24rev = k13rev * r.k23rev / k14rev;

	// _Based on formation of full complex
	r.k21rev = k14rev * r.k22rev * r.k24rev / r.kfwd / r.k15rev / r.k18rev * kfbnd;
	r.k20rev = k14rev * r.k22rev * r.k24rev / r.k19rev / r.k15rev;

	// _One detailed balance IL7/9 loop
	r.k32rev = r.k29rev * r.k31rev / k30rev;
	r.k28rev = k25rev * r.k27rev / k26rev;

	// Set the rates
	r.endo = tfR[0];
	r.activeEndo = tfR[1];
	r.sortF = tfR[2];
	r.kRec = tfR[3];
	r.kDeg = tfR[4];

	std::copy_n(tfR + 5, 6, r.Rexpr.begin());

	return r;
}

void dy_dt(const double * const y, const ratesS * const r, double * const dydt, double IL2, double IL15, double IL7, double IL9) {
	// IL2 in nM
	double IL2Ra = y[0];
	double IL2Rb = y[1];
	double gc = y[2];
	double IL2_IL2Ra = y[3];
	double IL2_IL2Rb = y[4];
	double IL2_gc = y[5];
	double IL2_IL2Ra_IL2Rb = y[6];
	double IL2_IL2Ra_gc = y[7];
	double IL2_IL2Rb_gc = y[8];
	double IL2_IL2Ra_IL2Rb_gc = y[9];
	
	// IL15 in nM
	double IL15Ra = y[10];
	double IL15_IL15Ra = y[11];
	double IL15_IL2Rb = y[12];
	double IL15_gc = y[13];
	double IL15_IL15Ra_IL2Rb = y[14];
	double IL15_IL15Ra_gc = y[15];
	double IL15_IL2Rb_gc = y[16];
	double IL15_IL15Ra_IL2Rb_gc = y[17];
	
	// IL7, IL9 in nM
	double IL7Ra = y[18];
	double IL7Ra_IL7 = y[19];
	double gc_IL7 = y[20];
	double IL7Ra_gc_IL7 = y[21];
	double IL9R = y[22];
	double IL9R_IL9 = y[23];
	double gc_IL9 = y[24];
	double IL9R_gc_IL9 = y[25];
	
	// IL2
	dydt[0] = -kfbnd * IL2Ra * IL2 + k1rev * IL2_IL2Ra - r->kfwd * IL2Ra * IL2_gc + r->k6rev * IL2_IL2Ra_gc - r->kfwd * IL2Ra * IL2_IL2Rb_gc + r->k8rev * IL2_IL2Ra_IL2Rb_gc - r->kfwd * IL2Ra * IL2_IL2Rb + r->k12rev * IL2_IL2Ra_IL2Rb;
	dydt[1] = -kfbnd * IL2Rb * IL2 + k2rev * IL2_IL2Rb - r->kfwd * IL2Rb * IL2_gc + r->k7rev * IL2_IL2Rb_gc - r->kfwd * IL2Rb * IL2_IL2Ra_gc + r->k9rev * IL2_IL2Ra_IL2Rb_gc - r->kfwd * IL2Rb * IL2_IL2Ra + r->k11rev * IL2_IL2Ra_IL2Rb;
	dydt[2] = -k3fwd * IL2 * gc + k3rev * IL2_gc - r->kfwd * IL2_IL2Rb * gc + r->k5rev * IL2_IL2Rb_gc - r->kfwd * IL2_IL2Ra * gc + r->k4rev * IL2_IL2Ra_gc - r->kfwd * IL2_IL2Ra_IL2Rb * gc + r->k10rev * IL2_IL2Ra_IL2Rb_gc;
	dydt[3] = -r->kfwd * IL2_IL2Ra * IL2Rb + r->k11rev * IL2_IL2Ra_IL2Rb - r->kfwd * IL2_IL2Ra * gc + r->k4rev * IL2_IL2Ra_gc + kfbnd * IL2 * IL2Ra - k1rev * IL2_IL2Ra;
	dydt[4] = -r->kfwd * IL2_IL2Rb * IL2Ra + r->k12rev * IL2_IL2Ra_IL2Rb - r->kfwd * IL2_IL2Rb * gc + r->k5rev * IL2_IL2Rb_gc + kfbnd * IL2 * IL2Rb - k2rev * IL2_IL2Rb;
	dydt[5] = -r->kfwd * IL2_gc * IL2Ra + r->k6rev * IL2_IL2Ra_gc - r->kfwd * IL2_gc * IL2Rb + r->k7rev * IL2_IL2Rb_gc + k3fwd * IL2 * gc - k3rev * IL2_gc;
	dydt[6] = -r->kfwd * IL2_IL2Ra_IL2Rb * gc + r->k10rev * IL2_IL2Ra_IL2Rb_gc + r->kfwd * IL2_IL2Ra * IL2Rb - r->k11rev * IL2_IL2Ra_IL2Rb + r->kfwd * IL2_IL2Rb * IL2Ra - r->k12rev * IL2_IL2Ra_IL2Rb;
	dydt[7] = -r->kfwd * IL2_IL2Ra_gc * IL2Rb + r->k9rev * IL2_IL2Ra_IL2Rb_gc + r->kfwd * IL2_IL2Ra * gc - r->k4rev * IL2_IL2Ra_gc + r->kfwd * IL2_gc * IL2Ra - r->k6rev * IL2_IL2Ra_gc;
	dydt[8] = -r->kfwd * IL2_IL2Rb_gc * IL2Ra + r->k8rev * IL2_IL2Ra_IL2Rb_gc + r->kfwd * gc * IL2_IL2Rb - r->k5rev * IL2_IL2Rb_gc + r->kfwd * IL2_gc * IL2Rb - r->k7rev * IL2_IL2Rb_gc;
	dydt[9] = r->kfwd * IL2_IL2Rb_gc * IL2Ra - r->k8rev * IL2_IL2Ra_IL2Rb_gc + r->kfwd * IL2_IL2Ra_gc * IL2Rb - r->k9rev * IL2_IL2Ra_IL2Rb_gc + r->kfwd * IL2_IL2Ra_IL2Rb * gc - r->k10rev * IL2_IL2Ra_IL2Rb_gc;

	// IL15
	dydt[10] = -kfbnd * IL15Ra * IL15 + k13rev * IL15_IL15Ra - kfbnd * IL15Ra * IL15_gc + r->k18rev * IL15_IL15Ra_gc - r->kfwd * IL15Ra * IL15_IL2Rb_gc + r->k20rev * IL15_IL15Ra_IL2Rb_gc - r->kfwd * IL15Ra * IL15_IL2Rb + r->k24rev * IL15_IL15Ra_IL2Rb;
	dydt[11] = -r->kfwd * IL15_IL15Ra * IL2Rb + r->k23rev * IL15_IL15Ra_IL2Rb - r->kfwd * IL15_IL15Ra * gc + r->k16rev * IL15_IL15Ra_gc + kfbnd * IL15 * IL15Ra - k13rev * IL15_IL15Ra;
	dydt[12] = -r->kfwd * IL15_IL2Rb * IL15Ra + r->k24rev * IL15_IL15Ra_IL2Rb - kfbnd * IL15_IL2Rb * gc + r->k17rev * IL15_IL2Rb_gc + kfbnd * IL15 * IL2Rb - k14rev * IL15_IL2Rb;
	dydt[13] = -kfbnd * IL15_gc * IL15Ra + r->k18rev * IL15_IL15Ra_gc - r->kfwd * IL15_gc * IL2Rb + r->k19rev * IL15_IL2Rb_gc + kfbnd * IL15 * gc - r->k15rev * IL15_gc;
	dydt[14] = -r->kfwd * IL15_IL15Ra_IL2Rb * gc + r->k22rev * IL15_IL15Ra_IL2Rb_gc + r->kfwd * IL15_IL15Ra * IL2Rb - r->k23rev * IL15_IL15Ra_IL2Rb + r->kfwd * IL15_IL2Rb * IL15Ra - r->k24rev * IL15_IL15Ra_IL2Rb;
	dydt[15] = -r->kfwd * IL15_IL15Ra_gc * IL2Rb + r->k21rev * IL15_IL15Ra_IL2Rb_gc + r->kfwd * IL15_IL15Ra * gc - r->k16rev * IL15_IL15Ra_gc + kfbnd * IL15_gc * IL15Ra - r->k18rev * IL15_IL15Ra_gc;
	dydt[16] = -r->kfwd * IL15_IL2Rb_gc * IL15Ra + r->k20rev * IL15_IL15Ra_IL2Rb_gc + kfbnd * gc * IL15_IL2Rb - r->k17rev * IL15_IL2Rb_gc + r->kfwd * IL15_gc * IL2Rb - r->k19rev * IL15_IL2Rb_gc;
	dydt[17] =  r->kfwd * IL15_IL2Rb_gc * IL15Ra - r->k20rev * IL15_IL15Ra_IL2Rb_gc + r->kfwd * IL15_IL15Ra_gc * IL2Rb - r->k21rev * IL15_IL15Ra_IL2Rb_gc + r->kfwd * IL15_IL15Ra_IL2Rb * gc - r->k22rev * IL15_IL15Ra_IL2Rb_gc;
	
	dydt[1] = dydt[1] - kfbnd * IL2Rb * IL15 + k14rev * IL15_IL2Rb - r->kfwd * IL2Rb * IL15_gc + r->k19rev * IL15_IL2Rb_gc - r->kfwd * IL2Rb * IL15_IL15Ra_gc + r->k21rev * IL15_IL15Ra_IL2Rb_gc - r->kfwd * IL2Rb * IL15_IL15Ra + r->k23rev * IL15_IL15Ra_IL2Rb;
	dydt[2] = dydt[2] - kfbnd * IL15 * gc + r->k15rev * IL15_gc - kfbnd * IL15_IL2Rb * gc + r->k17rev * IL15_IL2Rb_gc - r->kfwd * IL15_IL15Ra * gc + r->k16rev * IL15_IL15Ra_gc - r->kfwd * IL15_IL15Ra_IL2Rb * gc + r->k22rev * IL15_IL15Ra_IL2Rb_gc;
	
	// IL7
	dydt[2] = dydt[2] - kfbnd * IL7 * gc + k26rev * gc_IL7 - r->kfwd * gc * IL7Ra_IL7 + r->k27rev * IL7Ra_gc_IL7;
	dydt[18] = -kfbnd * IL7Ra * IL7 + k25rev * IL7Ra_IL7 - r->kfwd * IL7Ra * gc_IL7 + r->k28rev * IL7Ra_gc_IL7;
	dydt[19] = kfbnd * IL7Ra * IL7 - k25rev * IL7Ra_IL7 - r->kfwd * gc * IL7Ra_IL7 + r->k27rev * IL7Ra_gc_IL7;
	dydt[20] = -r->kfwd * IL7Ra * gc_IL7 + r->k28rev * IL7Ra_gc_IL7 + kfbnd * IL7 * gc - k26rev * gc_IL7;
	dydt[21] = r->kfwd * IL7Ra * gc_IL7 - r->k28rev * IL7Ra_gc_IL7 + r->kfwd * gc * IL7Ra_IL7 - r->k27rev * IL7Ra_gc_IL7;

	// IL9
	dydt[2] = dydt[2] - kfbnd * IL9 * gc + k30rev * gc_IL9 - r->kfwd * gc * IL9R_IL9 + r->k31rev * IL9R_gc_IL9;
	dydt[22] = -kfbnd * IL9R * IL9 + r->k29rev * IL9R_IL9 - r->kfwd * IL9R * gc_IL9 + r->k32rev * IL9R_gc_IL9;
	dydt[23] = kfbnd * IL9R * IL9 - r->k29rev * IL9R_IL9 - r->kfwd * gc * IL9R_IL9 + r->k31rev * IL9R_gc_IL9;
	dydt[24] = -r->kfwd * IL9R * gc_IL9 + r->k32rev * IL9R_gc_IL9 + kfbnd * IL9 * gc - k30rev * gc_IL9;
	dydt[25] = r->kfwd * IL9R * gc_IL9 - r->k32rev * IL9R_gc_IL9 + r->kfwd * gc * IL9R_IL9 - r->k31rev * IL9R_gc_IL9;
}


void jacobian(const double * const y, const ratesS * const r, double * const dydt, double IL2, double IL15, double IL7, double IL9) {
    // IL2 in nM
	double IL2Ra = y[0];
	double IL2Rb = y[1];
	double gc = y[2];
	double IL2_IL2Ra = y[3];
	double IL2_IL2Rb = y[4];
	double IL2_gc = y[5];
	double IL2_IL2Ra_IL2Rb = y[6];
	double IL2_IL2Ra_gc = y[7];
	double IL2_IL2Rb_gc = y[8];
	
	// IL15 in nM
	double IL15Ra = y[10];
	double IL15_IL15Ra = y[11];
	double IL15_IL2Rb = y[12];
	double IL15_gc = y[13];
	double IL15_IL15Ra_IL2Rb = y[14];
	double IL15_IL15Ra_gc = y[15];
	double IL15_IL2Rb_gc = y[16];
	
	// IL7, IL9 in nM
	double IL7Ra = y[18];
	double IL7Ra_IL7 = y[19];
	double gc_IL7 = y[20];
	double IL9R = y[22];
	double IL9R_IL9 = y[23];
	double gc_IL9 = y[24];
	
    array<array<double, 26>, 26> out;
    
    // unless otherwise specified, assume all partial derivatives are 0
    for (array<double, 26> &aa : out)
    	fill(aa.begin(), aa.end(), 0.0);
        
    // IL2Ra
    
    // partial derivative of IL2Ra with respect to IL2Ra
    out[0][0] = -kfbnd * IL2 - r->kfwd * IL2_gc - r->kfwd * IL2_IL2Rb_gc - r->kfwd * IL2_IL2Rb;
    
    // partial derivative of IL2Ra with respect to IL2_IL2Ra
    out[0][3] = k1rev;
    
    // partial derivative of IL2Ra with respect to IL2_IL2Rb
    out[0][4] = - r->kfwd * IL2Ra;
    
    // partial derivative of IL2Ra with respect to IL2_gc
    out[0][5] = - r->kfwd * IL2Ra;
    
    // partial derivative of IL2Ra with respect to IL2_IL2Ra_IL2Rb
    out[0][6] = r->k12rev;
    
    // partial derivative of IL2Ra with respect to IL2_IL2Ra_gc
    out[0][7] = r->k6rev;
    
    // partial derivative of IL2Ra with respect to IL2_IL2Rb_gc
    out[0][8] = - r->kfwd * IL2Ra;
    
    // partial derivative of IL2Ra with respect to IL2_IL2Ra_IL2Rb_gc
    out[0][9] = r->k8rev;
    
    // IL2Rb
    
    // partial derivative of IL2Rb with respect to IL2Rb
    out[1][1] = -kfbnd * IL2 - r->kfwd * IL2_gc - r->kfwd * IL2_IL2Ra_gc - r->kfwd * IL2_IL2Ra - kfbnd * IL15 - r->kfwd * IL15_gc - r->kfwd * IL15_IL15Ra_gc - r->kfwd * IL15_IL15Ra;
    
    // partial derivative of IL2Rb with respect to IL2_IL2Ra
    out[1][3] = - r->kfwd * IL2Rb;
    
    // partial derivative of IL2Rb with respect to IL2_IL2Rb
    out[1][4] = k2rev;
    
    // partial derivative of IL2Rb with respect to IL2_gc
    out[1][5] = - r->kfwd * IL2Rb;
    
    // partial derivative of IL2Rb with respect to IL2_IL2Ra_IL2Rb
    out[1][6] = r->k11rev;
    
    // partial derivative of IL2Rb with respect to IL2_IL2Ra_gc
    out[1][7] = - r->kfwd * IL2Rb;
    
    // partial derivative of IL2Rb with respect to IL2_IL2Rb_gc
    out[1][8] = r->k7rev;
    
    // partial derivative of IL2Rb with respect to IL2_IL2Ra_IL2Rb_gc
    out[1][9] = r->k9rev;
    
    // partial derivative of IL2Rb with respect to IL15_IL15Ra
    out[1][11] = - r->kfwd * IL2Rb;
    
    // partial derivative of IL2Rb with respect to IL15_IL2Rb
    out[1][12] = k14rev;
    
    // partial derivative of IL2Rb with respect to IL15_gc
    out[1][13] = - r->kfwd * IL2Rb;
    
    // partial derivative of IL2Rb with respect to IL15_IL15Ra_IL2Rb
    out[1][14] = r->k23rev;
    
    // partial derivative of IL2Rb with respect to IL15_IL15Ra_gc
    out[1][15] = - r->kfwd * IL2Rb;
    
    // partial derivative of IL2Rb with respect to IL15_IL2Rb_gc
    out[1][16] = r->k19rev;
    
    // partial derivative of IL2Rb with respect to IL15_IL15Ra_IL2Rb_gc
    out[1][17] = r->k21rev;
    
    // gc
    
    // partial derivative of gc with respect to gc
    out[2][2] = -k3fwd * IL2 - r->kfwd * IL2_IL2Rb - r->kfwd * IL2_IL2Ra - r->kfwd * IL2_IL2Ra_IL2Rb - kfbnd * IL15 - kfbnd * IL15_IL2Rb - r->kfwd * IL15_IL15Ra - r->kfwd * IL15_IL15Ra_IL2Rb - kfbnd * IL7 - r->kfwd * IL7Ra_IL7 - kfbnd * IL9 - r->kfwd * IL9R_IL9;
    
    // partial derivative of gc with respect to IL2_IL2Ra
    out[2][3] = - r->kfwd * gc;
    
    // partial derivative of gc with respect to IL2_IL2Rb
    out[2][4] = - r->kfwd * gc;
    
    // partial derivative of gc with respect to IL2_gc
    out[2][5] = k3rev;
    
    // partial derivative of gc with respect to IL2_IL2Ra_IL2Rb
    out[2][6] = - r->kfwd * gc;
    
    // partial derivative of gc with respect to IL2_IL2Ra_gc
    out[2][7] = r->k4rev;    
        
    // partial derivative of gc with respect to IL2_IL2Rb_gc
    out[2][8] = r->k5rev;
    
    // partial derivative of gc with respect to IL2_IL2Ra_IL2Rb_gc
    out[2][9] = r->k10rev;
    
    // partial derivative of gc with respect to IL15_IL15Ra
    out[2][11] = - r->kfwd;
    
    // partial derivative of gc with respect to IL15_IL2Rb
    out[2][12] = - kfbnd * gc;
    
    // partial derivative of gc with respect to IL15_gc
    out[2][13] = r->k15rev;
    
    // partial derivative of gc with respect to IL15_IL15Ra_IL2Rb
    out[2][14] = - r->kfwd;
    
    // partial derivative of gc with respect to IL15_IL15Ra_gc
    out[2][15] = r->k16rev;
    
    // partial derivative of gc with respect to IL15_IL2Rb_gc
    out[2][16] = r->k17rev;
    
    // partial derivative of gc with respect to IL15_IL15Ra_IL2Rb_gc
    out[2][17] = r->k22rev;
    
    // partial derivative of gc with respect to IL7Ra_IL7
    out[2][19] = - r->kfwd * gc;
    
    // partial derivative of gc with respect to gc_IL7
    out[2][20] = k26rev;
    
    // partial derivative of gc with respect to IL7Ra_gc_IL7
    out[2][21] = r->k27rev;
    
    // partial derivative of gc with respect to IL9R_IL9
    out[2][23] = - r->kfwd * gc;
    
    // partial derivative of gc with respect to gc_IL9
    out[2][24] = k30rev;
    
    // partial derivative of gc with respect to IL9R_gc_IL9
    out[2][25] = r->k31rev;
    
    // IL2_IL2Ra
    
    // partial derivative of IL2_IL2Ra with respect to IL2Ra
    out[3][0] = kfbnd * IL2;
    
    // partial derivative of IL2_IL2Ra with respect to IL2Rb
    out[3][1] = -r->kfwd * IL2_IL2Ra;
    
    // partial derivative of IL2_IL2Ra with respect to gc
    out[3][2] = - r->kfwd * IL2_IL2Ra;
    
    // partial derivative of IL2_IL2Ra with respect to IL2_IL2Ra
    out[3][3] = -r->kfwd * IL2Rb - r->kfwd * gc - k1rev;
    
    // partial derivative of IL2_IL2Ra with respect to IL2_IL2Ra_IL2Rb
    out[3][6] = r->k11rev;
    
    // partial derivative of IL2_IL2Ra with respect to IL2_IL2Ra_gc
    out[3][7] = r->k4rev;
    
    // IL2_IL2Rb
    
    // partial derivative of IL2_IL2Rb with respect to IL2Ra
    out[4][0] = -r->kfwd * IL2_IL2Rb;
    
    // partial derivative of IL2_IL2Rb with respect to IL2Rb
    out[4][1] = kfbnd * IL2;
    
    // partial derivative of IL2_IL2Rb with respect to gc
    out[4][2] = - r->kfwd * IL2_IL2Rb;
    
    // partial derivative of IL2_IL2Rb with respect to IL2_IL2Rb
    out[4][4] = -r->kfwd * IL2Ra - r->kfwd * gc - k2rev;
    
    // partial derivative of IL2_IL2Rb with respect to IL2_IL2Ra_IL2Rb
    out[4][6] = r->k12rev;
    
    // partial derivative of IL2_IL2Rb with respect to IL2_IL2Rb_gc
    out[4][8] = r->k5rev;
        
    // IL2_gc
    
    // partial derivative of IL2_gc with respect to IL2Ra
    out[5][0] = -r->kfwd * IL2_gc;
    
    // partial derivative of IL2_gc with respect to IL2Rb
    out[5][1] = - r->kfwd * IL2_gc;
    
    // partial derivative of IL2_gc with respect to gc
    out[5][2] = k3fwd * IL2;
    
    // partial derivative of IL2_gc with respect to IL2_gc
    out[5][5] = -r->kfwd * IL2Ra - r->kfwd * IL2Rb - k3rev;
    
    // partial derivative of IL2_gc with respect to IL2_IL2Ra_gc
    out[5][7] = r->k6rev;
    
    // partial derivative of IL2_gc with respect to IL2_IL2Rb_gc
    out[5][8] = r->k7rev;
    
    // IL2_IL2Ra_IL2Rb
    
    // partial derivative of IL2_IL2Ra_IL2Rb with respect to IL2Ra
    out[6][0] = r->kfwd * IL2_IL2Rb;
    
    // partial derivative of IL2_IL2Ra_IL2Rb with respect to IL2Rb
    out[6][1] = r->kfwd * IL2_IL2Ra;
    
    // partial derivative of IL2_IL2Ra_IL2Rb with respect to gc
    out[6][2] = -r->kfwd * IL2_IL2Ra_IL2Rb;
    
    // partial derivative of IL2_IL2Ra_IL2Rb with respect to IL2_IL2Ra
    out[6][3] = r->kfwd * IL2Rb; 
    
    // partial derivative of IL2_IL2Ra_IL2Rb with respect to IL2_IL2Rb
    out[6][4] = r->kfwd * IL2Ra;
    
    // partial derivative of IL2_IL2Ra_IL2Rb with respect to IL2_IL2Ra_IL2Rb
    out[6][6] = -r->kfwd * gc - r->k11rev - r->k12rev;
    
    // partial derivative of IL2_IL2Ra_IL2Rb with respect to IL2_IL2Ra_IL2Rb_gc
    out[6][9] = r->k10rev;
    
    // IL2_IL2Ra_gc
    
    // partial derivative of IL2_IL2Ra_gc with respect to IL2Ra
    out[7][0] = r->kfwd * IL2_gc;
    
    // partial derivative of IL2_IL2Ra_gc with respect to IL2Rb
    out[7][1] = -r->kfwd * IL2_IL2Ra_gc;
    
    // partial derivative of IL2_IL2Ra_gc with respect to gc
    out[7][2] = r->kfwd * IL2_IL2Ra;
    
    // partial derivative of IL2_IL2Ra_gc with respect to IL2_IL2Ra
    out[7][3] = r->kfwd * gc;
    
    // partial derivative of IL2_IL2Ra_gc with respect to IL2_gc
    out[7][5] = r->kfwd * IL2Ra;
    
    // partial derivative of IL2_IL2Ra_gc with respect to IL2_IL2Ra_gc
    out[7][7] = -r->kfwd * IL2Rb - r->k4rev - r->k6rev;
    
    // partial derivative of IL2_IL2Ra_gc with respect to IL2_IL2Ra_IL2Rb_gc
    out[7][9] = r->k9rev;
    
    // IL2_IL2Rb_gc
    
    // partial derivative for IL2_IL2Rb_gc with respect to IL2Ra
    out[8][0] = -r->kfwd * IL2_IL2Rb_gc;
    
    // partial derivative for IL2_IL2Rb_gc with respect to IL2Rb
    out[8][1] = r->kfwd * IL2_gc;
    
    // partial derivative for IL2_IL2Rb_gc with respect to gc
    out[8][2] = r->kfwd * IL2_IL2Rb;
    
    // partial derivative for IL2_IL2Rb_gc with respect to IL2_IL2Rb
    out[8][4] = r->kfwd * gc;
    
    // partial derivative for IL2_IL2Rb_gc with respect to IL2_gc
    out[8][5] = r->kfwd * IL2Rb;
    
    // partial derivative for IL2_IL2Rb_gc with respect to IL2_IL2Rb_gc
    out[8][8] = -r->kfwd * IL2Ra - r->k5rev - r->k7rev;
    
    // partial derivative for IL2_IL2Rb_gc with respect to IL2_IL2Ra_IL2Rb_gc
    out[8][9] = r->k8rev;
    
    // IL2_IL2Ra_IL2Rb_gc

    // partial derivative of IL2_IL2Ra_IL2Rb_gc with respect to IL2Ra
    out[9][0] = r->kfwd * IL2_IL2Rb_gc;
    
    // partial derivative of IL2_IL2Ra_IL2Rb_gc with respect to IL2Rb
    out[9][1] = r->kfwd * IL2_IL2Ra_gc;
    
    // partial derivative of IL2_IL2Ra_IL2Rb_gc with respect to gc
    out[9][2] = r->kfwd * IL2_IL2Ra_IL2Rb;
    
    // partial derivative of IL2_IL2Ra_IL2Rb_gc with respect to IL2_IL2Ra_IL2Rb
    out[9][6] = r->kfwd * gc;
    
    // partial derivative of IL2_IL2Ra_IL2Rb_gc with respect to IL2_IL2Ra_gc
    out[9][7] = r->kfwd * IL2Rb;
    
    // partial derivative of IL2_IL2Ra_IL2Rb_gc with respect to IL2_IL2Rb_gc
    out[9][8] = r->kfwd * IL2Ra;
    
    // partial derivative of IL2_IL2Ra_IL2Rb_gc with respect to IL2_IL2Ra_IL2Rb_gc
    out[9][9] = - r->k8rev - r->k9rev - r->k10rev;
    
    // IL15Ra
    
    // partial derivative of IL15Ra with respect to IL15Ra
    out[10][10] = -kfbnd * IL15 - kfbnd * IL15_gc - r->kfwd * IL15_IL2Rb_gc - r->kfwd * IL15_IL2Rb;
    
    // partial derivative of IL15Ra with respect to IL15_IL15Ra
    out[10][11] = k13rev;
    
    // partial derivative of IL15Ra with respect to IL15_IL2Rb
    out[10][12] = - r->kfwd * IL15Ra;
    
    // partial derivative of IL15Ra with respect to IL15_gc
    out[10][13] = - kfbnd * IL15Ra;
    
    // partial derivative of IL15Ra with respect to IL15_IL15Ra_IL2Rb
    out[10][14] = r->k24rev;
    
    // partial derivative of IL15Ra with respect to IL15_IL15Ra_gc
    out[10][15] = r->k18rev;
    
    // partial derivative of IL15Ra with respect to IL15_IL2Rb_gc
    out[10][16] = - r->kfwd * IL15Ra;
    
    // partial derivative of IL15Ra with respect to IL15_IL15Ra_IL2Rb_gc
    out[10][17] = r->k20rev;
    
    // IL15_IL15Ra
    
    // partial derivative of IL15_IL15Ra with respect to IL2Rb
    out[11][1] = -r->kfwd * IL15_IL15Ra;
    
    // partial derivative of IL15_IL15Ra with respect to gc
    out[11][2] = - r->kfwd * IL15_IL15Ra;
    
    // partial derivative of IL15_IL15Ra with respect to IL15Ra
    out[11][10] = kfbnd * IL15;
    
    // partial derivative of IL15_IL15Ra with respect to IL15_IL15Ra
    out[11][11] = -r->kfwd * IL2Rb- r->kfwd * gc - k13rev;
    
    // partial derivative of IL15_IL15Ra with respect to IL15_IL15Ra_IL2Rb
    out[11][14] = r->k23rev;
    
    // partial derivative of IL15_IL15Ra with respect to IL15_IL15Ra_gc
    out[11][15] = r->k16rev;
    
    // IL15_IL2Rb
    
    // partial derivative of IL15_IL2Rb with respect to IL2Rb
    out[12][1] = kfbnd * IL15;
    
    // partial derivative of IL15_IL2Rb with respect to gc
    out[12][2] = - kfbnd * IL15_IL2Rb;
    
    // partial derivative of IL15_IL2Rb with respect to IL15Ra
    out[12][10] = -r->kfwd * IL15_IL2Rb;
    
    // partial derivative of IL15_IL2Rb with respect to IL15_IL2Rb
    out[12][12] = -r->kfwd * IL15Ra - kfbnd * gc - k14rev;
    
    // partial derivative of IL15_IL2Rb with respect to IL15_IL15Ra_IL2Rb
    out[12][14] = r->k24rev;
    
    // partial derivative of IL15_IL2Rb with respect to IL15_IL2Rb_gc
    out[12][16] = r->k17rev;
    
    // IL15_gc
    
    // partial derivative of IL15_gc with respect to IL2Rb
    out[13][1] = - r->kfwd * IL15_gc;
    
    // partial derivative of IL15_gc with respect to gc
    out[13][2] = kfbnd * IL15;
    
    // partial derivative of IL15_gc with respect to IL15Ra
    out[13][10] = -kfbnd * IL15_gc;
    
    // partial derivative of IL15_gc with respect to IL15_gc
    out[13][13] = -kfbnd * IL15Ra - r->kfwd * IL2Rb - r->k15rev * IL15_gc;
    
    // partial derivative of IL15_gc with respect to IL15_IL15Ra_gc
    out[13][15] = r->k18rev;
    
    // partial derivative of IL15_gc with respect to IL15_IL2Rb_gc
    out[13][16] = r->k19rev; 
    
    // IL15_IL15Ra_IL2Rb
    
    // partial derivative of IL15_IL15Ra_IL2Rb with respect to IL2Rb
    out[14][1] = r->kfwd * IL15_IL15Ra;
    
    // partial derivative of IL15_IL15Ra_IL2Rb with respect to gc
    out[14][2] = -r->kfwd * IL15_IL15Ra_IL2Rb;
    
    // partial derivative of IL15_IL15Ra_IL2Rb with respect to IL15Ra
    out[14][10] = r->kfwd * IL15_IL2Rb;
    
    // partial derivative of IL15_IL15Ra_IL2Rb with respect to IL15_IL15Ra
    out[14][11] = r->kfwd * IL2Rb;
    
    // partial derivative of IL15_IL15Ra_IL2Rb with respect to IL15_IL2Rb
    out[14][12] = r->kfwd * IL15Ra;
    
    // partial derivative of IL15_IL15Ra_IL2Rb with respect to IL15_IL15Ra_IL2Rb
    out[14][14] = -r->kfwd * gc - r->k23rev - r->k24rev;
    
    // partial derivative of IL15_IL15Ra_IL2Rb with respect to IL15_IL15Ra_IL2Rb_gc
    out[14][17] = r->k22rev;
    
    // IL15_IL15Ra_gc
    
    // partial derivative of IL15_IL15Ra_gc with respect to IL2Rb
    out[15][1] = -r->kfwd * IL15_IL15Ra_gc;
    
    // partial derivative of IL15_IL15Ra_gc with respect to gc
    out[15][2] = r->kfwd * IL15_IL15Ra;
    
    // partial derivative of IL15_IL15Ra_gc with respect to IL15Ra
    out[15][10] = kfbnd * IL15_gc;
    
    // partial derivative of IL15_IL15Ra_gc with respect to IL15_IL15Ra
    out[15][11] = r->kfwd * gc;
    
    // partial derivative of IL15_IL15Ra_gc with respect to IL15_gc
    out[15][13] = kfbnd * IL15Ra;
    
    // partial derivative of IL15_IL15Ra_gc with respect to IL15_IL15Ra_gc
    out[15][15] = -r->kfwd * IL2Rb - r->k16rev - r->k18rev;
    
    // partial derivative of IL15_IL15Ra_gc with respect to IL15_IL15Ra_IL2Rb_gc
    out[15][17] = r->k21rev;
    
    // IL15_IL2Rb_gc
    
    // partial derivative of IL15_IL2Rb_gc with respect to IL2Rb
    out[16][1] = r->kfwd * IL15_gc;
    
    // partial derivative of IL15_IL2Rb_gc with respect to gc
    out[16][2] = kfbnd * IL15_IL2Rb;
    
    // partial derivative of IL15_IL2Rb_gc with respect to IL15Ra
    out[16][10] = -r->kfwd * IL15_IL2Rb_gc;
    
    // partial derivative of IL15_IL2Rb_gc with respect to IL15_IL2Rb
    out[16][12] = kfbnd * gc;
    
    // partial derivative of IL15_IL2Rb_gc with respect to IL15_gc
    out[16][13] = r->kfwd * IL2Rb;
    
    // partial derivative of IL15_IL2Rb_gc with respect to IL15_IL2Rb_gc
    out[16][16] = -r->kfwd * IL15Ra - r->k17rev - r->k19rev;
    
    // partial derivative of IL15_IL2Rb_gc with respect to IL15_IL15Ra_IL2Rb_gc
    out[16][17] = r->k20rev;
    
    // IL15_IL15Ra_IL2Rb_gc
    
    // partial derivative of IL15_IL15Ra_IL2Rb_gc with respect to IL2Rb
    out[17][1] = r->kfwd * IL15_IL15Ra_gc;
    
    // partial derivative of IL15_IL15Ra_IL2Rb_gc with respect to gc
    out[17][2] = r->kfwd * IL15_IL15Ra_IL2Rb;
    
    // partial derivative of IL15_IL15Ra_IL2Rb_gc with respect to IL15Ra
    out[17][10] = r->kfwd * IL15_IL2Rb_gc;
    
    // partial derivative of IL15_IL15Ra_IL2Rb_gc with respect to IL15_IL15Ra_IL2Rb
    out[17][14] = r->kfwd * gc;
    
    // partial derivative of IL15_IL15Ra_IL2Rb_gc with respect to IL15_IL15Ra_gc
    out[17][15] = r->kfwd * IL2Rb;
    
    // partial derivative of IL15_IL15Ra_IL2Rb_gc with respect to IL15_IL2Rb_gc
    out[17][16] = r->kfwd * IL15Ra;
    
    // partial derivative of IL15_IL15Ra_IL2Rb_gc with respect to IL15_IL15Ra_IL2Rb_gc
    out[17][17] = - r->k20rev - r->k21rev - r->k22rev;
    
    // IL7Ra 
    
    // partial derivative of IL7Ra with respect to IL7Ra 
    out[18][18] = -kfbnd * IL7 - r->kfwd * gc_IL7;
    
    // partial derivative of IL7Ra with respect to IL7Ra_IL7
    out[18][19] = k25rev;
    
    // partial derivative of IL7Ra with respect to gc_IL7
    out[18][20] = - r->kfwd * IL7Ra;
    
    // partial derivative of IL7Ra with respect to IL7Ra_gc_IL7
    out[18][21] = r->k28rev;
    
    // IL7Ra_IL7
    
    // partial derivative of IL7Ra_IL7 with respect to gc
    out[19][2] = - r->kfwd * IL7Ra_IL7;
    
    // partial derivative of IL7Ra_IL7 with respect to IL7Ra
    out[19][18] = kfbnd * IL7;
    
    // partial derivative of IL7Ra_IL7 with respect to IL7Ra_IL7
    out[19][19] = - k25rev - r->kfwd * gc;
    
    // partial derivative of IL7Ra_IL7 with respect to IL7Ra_gc_IL7
    out[19][21] = r->k27rev;
    
    // gc_IL7
    
    // partial derivative of gc_IL7 with respect to gc
    out[20][2] = kfbnd * IL7;
    
    // partial derivative of gc_IL7 with respect to IL7Ra
    out[20][18] = -r->kfwd * gc_IL7;
    
    // partial derivative of gc_IL7 with respect to gc_IL7
    out[20][20] = -r->kfwd * IL7Ra - k26rev;
    
    // partial derivative of gc_IL7 with respect to IL7Ra_gc_IL7
    out[20][21] = r->k28rev;
    
    // IL7Ra_gc_IL7
    
    // partial derivative of IL7Ra_gc_IL7 with respect to gc
    out[21][2] = r->kfwd * IL7Ra_IL7;
    
    // partial derivative of IL7Ra_gc_IL7 with respect to IL7Ra
    out[21][18] = r->kfwd * gc_IL7;
    
    // partial derivative of IL7Ra_gc_IL7 with respect to IL7Ra_IL7
    out[21][19] = r->kfwd * gc;
    
    // partial derivative of IL7Ra_gc_IL7 with respect to gc_IL7
    out[21][20] = r->kfwd * IL7Ra;
    
    // partial derivative of IL7Ra_gc_IL7 with respect to IL7Ra_gc_IL7
    out[21][21] = - r->k28rev - r->k27rev;
    
    // IL9R
    
    // partial derivative of IL9R with respect to IL9R
    out[22][22] = -kfbnd * IL9 - r->kfwd * gc_IL9;
    
    // partial derivative of IL9R with respect to IL9R_IL9
    out[22][23] = r->k29rev;
    
    // partial derivative of IL9R with respect to gc_IL9
    out[22][24] = - r->kfwd * IL9R;
    
    // partial derivative of IL9R with respect to IL9R_gc_IL9
    out[22][25] = r->k32rev;
    
    // IL9R_IL9 
   
    // partial derivative of IL9R_IL9 with respect to gc
    out[23][2] = - r->kfwd * IL9R_IL9;
    
    // partial derivative of IL9R_IL9 with respect to IL9R
    out[23][22] = kfbnd * IL9; 
    
    // partial derivative of IL9R_IL9 with respect to IL9R_IL9
    out[23][23] = - r->k29rev - r->kfwd * gc;
    
    // partial derivative of IL9R_IL9 with respect to IL9R_gc_IL9
    out[23][25] = r->k31rev;
    
    // gc_IL9

    // partial derivative of gc_IL9 with respect to gc
    out[24][2] = kfbnd * IL9;
    
    // partial derivative of gc_IL9 with respect to IL9R
    out[24][22] = -r->kfwd * gc_IL9;
    
    // partial derivative of gc_IL9 with respect to gc_IL9
    out[24][24] = -r->kfwd * IL9R - k30rev;
    
    // partial derivative of gc_IL9 with respect to IL9R_gc_IL9
    out[24][25] = r->k32rev;
    
    // IL9R_gc_IL9
   
    // partial derivative of IL9R_gc_IL9 with respect to gc
    out[25][2] = r->kfwd * IL9R_IL9;
    
    // partial derivative of IL9R_gc_IL9 with respect to IL9R
    out[25][22] = r->kfwd * gc_IL9;
    
    // partial derivative of IL9R_gc_IL9 with respect to IL9R_IL9
    out[25][23] = r->kfwd * gc;
    
    // partial derivative of IL9R_gc_IL9 with respect to gc_IL9
    out[25][24] = r->kfwd * IL9R;
    
    // partial derivative of IL9R_gc_IL9 with respect to IL9R_gc_IL9
    out[25][25] = - r->k32rev - r->k31rev;

    // Copy all the data from out to dydt
    for (size_t ii = 0; ii < out.size(); ii++)
    	copy(out[ii].begin(), out[ii].end(), dydt + ii*out.size());
}


extern "C" void jacobian_C(double *y_in, double, double *dydt_out, double *rxn_in) {
	array<double, 12> tfr;

	ratesS r = param(rxn_in, tfr.data());

	dy_dt(y_in, &r, dydt_out, r.IL2, r.IL15, r.IL7, r.IL9);
}