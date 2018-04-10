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
#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <cvode/cvode_direct.h>
#include <iostream>
#include "model.hpp"

using std::array;
using std::copy;
using std::vector;
using std::fill;
using std::string;

constexpr size_t Nparams = 26;

std::array<bool, 26> __active_species_IDX() {
	std::array<bool, 26> __active_species_IDX;
	std::fill(__active_species_IDX.begin(), __active_species_IDX.end(), false);

	__active_species_IDX[8] = true;
	__active_species_IDX[9] = true;
	__active_species_IDX[16] = true;
	__active_species_IDX[17] = true;
	__active_species_IDX[21] = true;
	__active_species_IDX[25] = true;

	return __active_species_IDX;
}

const std::array<bool, 26> activeV = __active_species_IDX();

ratesS param(const double * const rxntfR) {
	ratesS r;

	r.IL2 = rxntfR[0];
	r.IL15 = rxntfR[1];
	r.IL7 = rxntfR[2];
	r.IL9 = rxntfR[3];
	r.kfwd = rxntfR[4];
	r.k5rev = rxntfR[5];
	r.k6rev = rxntfR[6];
	r.k15rev = rxntfR[7];
	r.k17rev = rxntfR[8];
	r.k18rev = rxntfR[9];
	r.k22rev = rxntfR[10];
	r.k23rev = rxntfR[11];
	r.k27rev = rxntfR[12];
	r.k29rev = rxntfR[13];
	r.k31rev = rxntfR[14];

	// Set the rates
	r.endo = rxntfR[15];
	r.activeEndo = rxntfR[16];
	r.sortF = rxntfR[17];
	r.kRec = rxntfR[18];
	r.kDeg = rxntfR[19];

	std::copy_n(rxntfR + 20, 6, r.Rexpr.begin());

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

	// These are probably measured in the literature
	const double k10rev = 12.0 * r->k5rev / 1.5; // doi:10.1016/j.jmb.2004.04.038
	const double k11rev = 63.0 * r->k5rev / 1.5; // doi:10.1016/j.jmb.2004.04.038
	// To satisfy detailed balance these relationships should hold
	// _Based on initial assembly steps
	const double k4rev = kfbnd * r->k6rev * k3rev / k1rev / k3fwd;
	const double k7rev = k3fwd * k2rev * r->k5rev / kfbnd / k3rev;
	const double k12rev = k1rev * k11rev / k2rev;
	// _Based on formation of full complex
	const double k9rev = k2rev * k10rev * k12rev / kfbnd / k3rev / r->k6rev * k3fwd;
	const double k8rev = k2rev * k10rev * k12rev / kfbnd / k7rev / k3rev * k3fwd;

	// IL15
	// To satisfy detailed balance these relationships should hold
	// _Based on initial assembly steps
	const double k16rev = r->kfwd * r->k18rev * r->k15rev / k13rev / kfbnd;
	const double k19rev = r->kfwd * k14rev * r->k17rev / kfbnd / r->k15rev;
	const double k24rev = k13rev * r->k23rev / k14rev;

	// _Based on formation of full complex
	const double k21rev = k14rev * r->k22rev * k24rev / r->kfwd / r->k15rev / r->k18rev * kfbnd;
	const double k20rev = k14rev * r->k22rev * k24rev / k19rev / r->k15rev;

	// _One detailed balance IL7/9 loop
	const double k32rev = r->k29rev * r->k31rev / k30rev;
	const double k28rev = k25rev * r->k27rev / k26rev;
	
	// IL2
	dydt[0] = -kfbnd * IL2Ra * IL2 + k1rev * IL2_IL2Ra - r->kfwd * IL2Ra * IL2_gc + r->k6rev * IL2_IL2Ra_gc - r->kfwd * IL2Ra * IL2_IL2Rb_gc + k8rev * IL2_IL2Ra_IL2Rb_gc - r->kfwd * IL2Ra * IL2_IL2Rb + k12rev * IL2_IL2Ra_IL2Rb;
	dydt[1] = -kfbnd * IL2Rb * IL2 + k2rev * IL2_IL2Rb - r->kfwd * IL2Rb * IL2_gc + k7rev * IL2_IL2Rb_gc - r->kfwd * IL2Rb * IL2_IL2Ra_gc + k9rev * IL2_IL2Ra_IL2Rb_gc - r->kfwd * IL2Rb * IL2_IL2Ra + k11rev * IL2_IL2Ra_IL2Rb;
	dydt[2] = -k3fwd * IL2 * gc + k3rev * IL2_gc - r->kfwd * IL2_IL2Rb * gc + r->k5rev * IL2_IL2Rb_gc - r->kfwd * IL2_IL2Ra * gc + k4rev * IL2_IL2Ra_gc - r->kfwd * IL2_IL2Ra_IL2Rb * gc + k10rev * IL2_IL2Ra_IL2Rb_gc;
	dydt[3] = -r->kfwd * IL2_IL2Ra * IL2Rb + k11rev * IL2_IL2Ra_IL2Rb - r->kfwd * IL2_IL2Ra * gc + k4rev * IL2_IL2Ra_gc + kfbnd * IL2 * IL2Ra - k1rev * IL2_IL2Ra;
	dydt[4] = -r->kfwd * IL2_IL2Rb * IL2Ra + k12rev * IL2_IL2Ra_IL2Rb - r->kfwd * IL2_IL2Rb * gc + r->k5rev * IL2_IL2Rb_gc + kfbnd * IL2 * IL2Rb - k2rev * IL2_IL2Rb;
	dydt[5] = -r->kfwd * IL2_gc * IL2Ra + r->k6rev * IL2_IL2Ra_gc - r->kfwd * IL2_gc * IL2Rb + k7rev * IL2_IL2Rb_gc + k3fwd * IL2 * gc - k3rev * IL2_gc;
	dydt[6] = -r->kfwd * IL2_IL2Ra_IL2Rb * gc + k10rev * IL2_IL2Ra_IL2Rb_gc + r->kfwd * IL2_IL2Ra * IL2Rb - k11rev * IL2_IL2Ra_IL2Rb + r->kfwd * IL2_IL2Rb * IL2Ra - k12rev * IL2_IL2Ra_IL2Rb;
	dydt[7] = -r->kfwd * IL2_IL2Ra_gc * IL2Rb + k9rev * IL2_IL2Ra_IL2Rb_gc + r->kfwd * IL2_IL2Ra * gc - k4rev * IL2_IL2Ra_gc + r->kfwd * IL2_gc * IL2Ra - r->k6rev * IL2_IL2Ra_gc;
	dydt[8] = -r->kfwd * IL2_IL2Rb_gc * IL2Ra + k8rev * IL2_IL2Ra_IL2Rb_gc + r->kfwd * gc * IL2_IL2Rb - r->k5rev * IL2_IL2Rb_gc + r->kfwd * IL2_gc * IL2Rb - k7rev * IL2_IL2Rb_gc;
	dydt[9] = r->kfwd * IL2_IL2Rb_gc * IL2Ra - k8rev * IL2_IL2Ra_IL2Rb_gc + r->kfwd * IL2_IL2Ra_gc * IL2Rb - k9rev * IL2_IL2Ra_IL2Rb_gc + r->kfwd * IL2_IL2Ra_IL2Rb * gc - k10rev * IL2_IL2Ra_IL2Rb_gc;

	// IL15
	dydt[10] = -kfbnd * IL15Ra * IL15 + k13rev * IL15_IL15Ra - kfbnd * IL15Ra * IL15_gc + r->k18rev * IL15_IL15Ra_gc - r->kfwd * IL15Ra * IL15_IL2Rb_gc + k20rev * IL15_IL15Ra_IL2Rb_gc - r->kfwd * IL15Ra * IL15_IL2Rb + k24rev * IL15_IL15Ra_IL2Rb;
	dydt[11] = -r->kfwd * IL15_IL15Ra * IL2Rb + r->k23rev * IL15_IL15Ra_IL2Rb - r->kfwd * IL15_IL15Ra * gc + k16rev * IL15_IL15Ra_gc + kfbnd * IL15 * IL15Ra - k13rev * IL15_IL15Ra;
	dydt[12] = -r->kfwd * IL15_IL2Rb * IL15Ra + k24rev * IL15_IL15Ra_IL2Rb - kfbnd * IL15_IL2Rb * gc + r->k17rev * IL15_IL2Rb_gc + kfbnd * IL15 * IL2Rb - k14rev * IL15_IL2Rb;
	dydt[13] = -kfbnd * IL15_gc * IL15Ra + r->k18rev * IL15_IL15Ra_gc - r->kfwd * IL15_gc * IL2Rb + k19rev * IL15_IL2Rb_gc + kfbnd * IL15 * gc - r->k15rev * IL15_gc;
	dydt[14] = -r->kfwd * IL15_IL15Ra_IL2Rb * gc + r->k22rev * IL15_IL15Ra_IL2Rb_gc + r->kfwd * IL15_IL15Ra * IL2Rb - r->k23rev * IL15_IL15Ra_IL2Rb + r->kfwd * IL15_IL2Rb * IL15Ra - k24rev * IL15_IL15Ra_IL2Rb;
	dydt[15] = -r->kfwd * IL15_IL15Ra_gc * IL2Rb + k21rev * IL15_IL15Ra_IL2Rb_gc + r->kfwd * IL15_IL15Ra * gc - k16rev * IL15_IL15Ra_gc + kfbnd * IL15_gc * IL15Ra - r->k18rev * IL15_IL15Ra_gc;
	dydt[16] = -r->kfwd * IL15_IL2Rb_gc * IL15Ra + k20rev * IL15_IL15Ra_IL2Rb_gc + kfbnd * gc * IL15_IL2Rb - r->k17rev * IL15_IL2Rb_gc + r->kfwd * IL15_gc * IL2Rb - k19rev * IL15_IL2Rb_gc;
	dydt[17] =  r->kfwd * IL15_IL2Rb_gc * IL15Ra - k20rev * IL15_IL15Ra_IL2Rb_gc + r->kfwd * IL15_IL15Ra_gc * IL2Rb - k21rev * IL15_IL15Ra_IL2Rb_gc + r->kfwd * IL15_IL15Ra_IL2Rb * gc - r->k22rev * IL15_IL15Ra_IL2Rb_gc;
	
	dydt[1] = dydt[1] - kfbnd * IL2Rb * IL15 + k14rev * IL15_IL2Rb - r->kfwd * IL2Rb * IL15_gc + k19rev * IL15_IL2Rb_gc - r->kfwd * IL2Rb * IL15_IL15Ra_gc + k21rev * IL15_IL15Ra_IL2Rb_gc - r->kfwd * IL2Rb * IL15_IL15Ra + r->k23rev * IL15_IL15Ra_IL2Rb;
	dydt[2] = dydt[2] - kfbnd * IL15 * gc + r->k15rev * IL15_gc - kfbnd * IL15_IL2Rb * gc + r->k17rev * IL15_IL2Rb_gc - r->kfwd * IL15_IL15Ra * gc + k16rev * IL15_IL15Ra_gc - r->kfwd * IL15_IL15Ra_IL2Rb * gc + r->k22rev * IL15_IL15Ra_IL2Rb_gc;
	
	// IL7
	dydt[2] = dydt[2] - kfbnd * IL7 * gc + k26rev * gc_IL7 - r->kfwd * gc * IL7Ra_IL7 + r->k27rev * IL7Ra_gc_IL7;
	dydt[18] = -kfbnd * IL7Ra * IL7 + k25rev * IL7Ra_IL7 - r->kfwd * IL7Ra * gc_IL7 + k28rev * IL7Ra_gc_IL7;
	dydt[19] = kfbnd * IL7Ra * IL7 - k25rev * IL7Ra_IL7 - r->kfwd * gc * IL7Ra_IL7 + r->k27rev * IL7Ra_gc_IL7;
	dydt[20] = -r->kfwd * IL7Ra * gc_IL7 + k28rev * IL7Ra_gc_IL7 + kfbnd * IL7 * gc - k26rev * gc_IL7;
	dydt[21] = r->kfwd * IL7Ra * gc_IL7 - k28rev * IL7Ra_gc_IL7 + r->kfwd * gc * IL7Ra_IL7 - r->k27rev * IL7Ra_gc_IL7;

	// IL9
	dydt[2] = dydt[2] - kfbnd * IL9 * gc + k30rev * gc_IL9 - r->kfwd * gc * IL9R_IL9 + r->k31rev * IL9R_gc_IL9;
	dydt[22] = -kfbnd * IL9R * IL9 + r->k29rev * IL9R_IL9 - r->kfwd * IL9R * gc_IL9 + k32rev * IL9R_gc_IL9;
	dydt[23] = kfbnd * IL9R * IL9 - r->k29rev * IL9R_IL9 - r->kfwd * gc * IL9R_IL9 + r->k31rev * IL9R_gc_IL9;
	dydt[24] = -r->kfwd * IL9R * gc_IL9 + k32rev * IL9R_gc_IL9 + kfbnd * IL9 * gc - k30rev * gc_IL9;
	dydt[25] = r->kfwd * IL9R * gc_IL9 - k32rev * IL9R_gc_IL9 + r->kfwd * gc * IL9R_IL9 - r->k31rev * IL9R_gc_IL9;
}


extern "C" void dydt_C(double *y_in, double, double *dydt_out, double *rxn_in) {
	ratesS r = param(rxn_in);

	dy_dt(y_in, &r, dydt_out, r.IL2, r.IL15, r.IL7, r.IL9);
}


/**
 * @brief      Solve for the ligand consumption rate in the endosome.
 *
 * @param      dydt  The rate of change vector solved for the receptor species.
 */
void findLigConsume(double *dydt) {
	// Calculate the ligand consumption.
	dydt[52] -= std::accumulate(dydt+3, dydt+10, 0) / internalV;
	dydt[53] -= std::accumulate(dydt+11, dydt+18, 0) / internalV;
	dydt[54] -= std::accumulate(dydt+19, dydt+22, 0) / internalV;
	dydt[55] -= std::accumulate(dydt+23, dydt+26, 0) / internalV;
}


void trafficking(const double * const y, const ratesS * const r, double * const dydt) {
	// Implement trafficking.
	size_t halfL = activeV.size();

	// Actually calculate the trafficking
	for (size_t ii = 0; ii < halfL; ii++) {
		if (activeV[ii]) {
			dydt[ii] += -y[ii]*(r->endo + r->activeEndo); // Endocytosis
			dydt[ii+halfL] += y[ii]*(r->endo + r->activeEndo)/internalFrac - r->kDeg*y[ii+halfL]; // Endocytosis, degradation
		} else {
			dydt[ii] += -y[ii]*r->endo + r->kRec*(1.0-r->sortF)*y[ii+halfL]*internalFrac; // Endocytosis, recycling
			dydt[ii+halfL] += y[ii]*r->endo/internalFrac - r->kRec*(1.0-r->sortF)*y[ii+halfL] - (r->kDeg*r->sortF)*y[ii+halfL]; // Endocytosis, recycling, degradation
		}
	}

	// Expression: IL2Ra, IL2Rb, gc, IL15Ra, IL7Ra, IL9R
	dydt[0] += r->Rexpr[0];
	dydt[1] += r->Rexpr[1];
	dydt[2] += r->Rexpr[2];
	dydt[10] += r->Rexpr[3];
	dydt[18] += r->Rexpr[4];
	dydt[22] += r->Rexpr[5];

	// Degradation does lead to some clearance of ligand in the endosome
	for (size_t ii = 0; ii < 4; ii++) {
		dydt[52 + ii] -= y[52 + ii] * r->kDeg;
	}
}


void fullModel(const double * const y, const ratesS * const r, double *dydt) {
	// Implement full model.
	fill(dydt, dydt + 56, 0.0);

	// Calculate cell surface and endosomal reactions
	dy_dt(y,      r, dydt,     r->IL2, r->IL15, r->IL7, r->IL9);
	dy_dt(y + 26, r, dydt + 26, y[52],   y[53],  y[54],  y[55]);

	// Handle trafficking
	trafficking(y, r, dydt);

	// Handle endosomal ligand balance.
	findLigConsume(dydt);
}


int fullModelCVode (const double, const N_Vector xx, N_Vector dxxdt, void *user_data) {
	ratesS rattes = param(static_cast<double *>(user_data));

	// Get the data in the right form
	fullModel(NV_DATA_S(xx), &rattes, NV_DATA_S(dxxdt));

	return 0;
}


extern "C" void fullModel_C(const double * const y_in, double, double *dydt_out, double *rxn_in) {
	ratesS r = param(rxn_in);

	fullModel(y_in, &r, dydt_out);
}

// TODO: Provide sensitivity of autocrine state.
array<double, 56> solveAutocrine(const ratesS * const r) {
	array<double, 56> y0;
	fill(y0.begin(), y0.end(), 0.0);

	array<size_t, 6> recIDX = {{0, 1, 2, 10, 18, 22}};

	// Expand out trafficking terms
	double kRec = r->kRec*(1-r->sortF);
	double kDeg = r->kDeg*r->sortF;

	// Assuming no autocrine ligand, so can solve steady state
	// Add the species
	for (size_t ii = 0; ii < recIDX.size(); ii++) {
		y0[recIDX[ii] + 26] = r->Rexpr[ii] / kDeg / internalFrac;
		y0[recIDX[ii]] = (r->Rexpr[ii] + kRec*y0[recIDX[ii] + 26]*internalFrac)/r->endo;
	}

	return y0;
}


static void errorHandler(int error_code, const char *module, const char *function, char *msg, void *) {
	if (error_code == CV_WARNING) return;

	std::cout << "Internal CVode error in " << function << std::endl;
	std::cout << msg << std::endl;
	std::cout << "In module: " << module << std::endl;
	std::cout << "Error code: " << error_code << std::endl;
}

struct solver {
	void *cvode_mem;
	SUNLinearSolver LS;
	N_Vector state;
	N_Vector *yS;
	SUNMatrix A;
	bool sensi;
};


void solverFree(solver *sMem) {
	if (sMem->sensi) {
		CVodeSensFree(sMem->cvode_mem);
	}

	N_VDestroy_Serial(sMem->state);
	CVodeFree(&sMem->cvode_mem);
	SUNLinSolFree(sMem->LS);
	SUNMatDestroy(sMem->A);
}


void solver_setup(solver *sMem, double *params, bool sensi) {
	// Copy in whether we're doing a sensitivity analysis
	sMem->sensi = sensi;

	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration */
	sMem->cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (sMem->cvode_mem == nullptr) {
		solverFree(sMem);
		throw std::runtime_error(string("Error calling CVodeCreate in solver_setup."));
	}
	
	CVodeSetErrHandlerFn(sMem->cvode_mem, &errorHandler, static_cast<void *>(params));

	/* Call CVodeInit to initialize the integrator memory and specify the
	 * user's right hand side function in y'=f(t,y), the inital time T0, and
	 * the initial dependent variable vector y. */
	if (CVodeInit(sMem->cvode_mem, fullModelCVode, 0.0, sMem->state) < 0) {
		solverFree(sMem);
		throw std::runtime_error(string("Error calling CVodeInit in solver_setup."));
	}
	
	N_Vector abbstol = N_VNew_Serial(NV_LENGTH_S(sMem->state));
	N_VConst(abstolIn, abbstol);
	
	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	if (CVodeSVtolerances(sMem->cvode_mem, reltolIn, abbstol) < 0) {
		N_VDestroy_Serial(abbstol);
		solverFree(sMem);
		throw std::runtime_error(string("Error calling CVodeSVtolerances in solver_setup."));
	}
	N_VDestroy_Serial(abbstol);

	sMem->A = SUNDenseMatrix(NV_LENGTH_S(sMem->state), NV_LENGTH_S(sMem->state));
	sMem->LS = SUNDenseLinearSolver(sMem->state, sMem->A);
	
	// Call CVDense to specify the CVDENSE dense linear solver
	// Also SUNSPBCGS and SUNSPTFQMR options
	if (CVDlsSetLinearSolver(sMem->cvode_mem, sMem->LS, sMem->A) < 0) {
		solverFree(sMem);
		throw std::runtime_error(string("Error calling CVDlsSetLinearSolver in solver_setup."));
	}
	
	// Pass along the parameter structure to the differential equations
	if (CVodeSetUserData(sMem->cvode_mem, static_cast<void *>(params)) < 0) {
		solverFree(sMem);
		throw std::runtime_error(string("Error calling CVodeSetUserData in solver_setup."));
	}

	CVodeSetMaxNumSteps(sMem->cvode_mem, 2000000);

	if (sensi) {
		// Set sensitivity initial conditions
		sMem->yS = N_VCloneVectorArray(Nparams, sMem->state);
		for (size_t is = 0; is < Nparams; is++) N_VConst(0.0, sMem->yS[is]);

		// Call CVodeSensInit1 to activate forward sensitivity computations
		// and allocate internal memory for CVODES related to sensitivity
		// calculations. Computes the right-hand sides of the sensitivity
		// ODE, one at a time
		if (CVodeSensInit(sMem->cvode_mem, Nparams, CV_SIMULTANEOUS, nullptr, sMem->yS) < 0) {
			solverFree(sMem);
			throw std::runtime_error(string("Error calling CVodeSensInit in solver_setup."));
		}

		array<double, Nparams> abs;
		fill(abs.begin(), abs.end(), 1.0E-2);

		// Call CVodeSensSStolerances to estimate tolerances for sensitivity 
		// variables based on the rolerances supplied for states variables and 
		// the scaling factor pbar
		if (CVodeSensSStolerances(sMem->cvode_mem, 1.0E-2, abs.data()) < 0) {
			solverFree(sMem);
			throw std::runtime_error(string("Error calling CVodeSensSStolerances in solver_setup."));
		}

		// Specify problem parameter information for sensitivity calculations
		if (CVodeSetSensParams(sMem->cvode_mem, params, nullptr, nullptr) < 0) {
			solverFree(sMem);
			throw std::runtime_error(string("Error calling CVodeSetSensParams in solver_setup."));
		}
	}
}


extern "C" int runCkine (double *tps, size_t ntps, double *out, double *rxnRatesIn, bool sensi) {
	ratesS rattes = param(rxnRatesIn);
	size_t itps = 0;

	array<double, 56> y0 = solveAutocrine(&rattes);

	solver sMem;

	if (tps[0] < std::numeric_limits<double>::epsilon()) {
		std::copy_n(y0.begin(), y0.size(), out);
		itps = 1;
	}

	// Fill output values with 0's
	fill(out, out + ntps*y0.size(), 0.0);

	// Just the full model
	sMem.state = N_VMake_Serial(static_cast<long>(y0.size()), y0.data());

	solver_setup(&sMem, rxnRatesIn, sensi);

	double tret = 0;

	for (; itps < ntps; itps++) {
		if (tps[itps] < tret) {
			std::cout << "Can't go backwards." << std::endl;
			solverFree(&sMem);
			return -1;
		}

		int returnVal = CVode(sMem.cvode_mem, tps[itps], sMem.state, &tret, CV_NORMAL);
		
		if (returnVal < 0) {
			std::cout << "CVode error in CVode. Code: " << returnVal << std::endl;
			solverFree(&sMem);
			return returnVal;
		}

		// Copy out result
		std::copy_n(NV_DATA_S(sMem.state), y0.size(), out + y0.size()*itps);
	}

	solverFree(&sMem);
	return 0;
}
