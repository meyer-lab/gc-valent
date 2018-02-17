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
#include <cvodes/cvodes.h>

using std::array;
using std::copy;
using std::copy_n;
using std::vector;
using std::fill;
using std::string;
using std::runtime_error;

const size_t plen = 28;

struct ratesS {
	std::array<double, plen> p;
	double *rxn;
	double *trafRates;
};


const double abstolIn = 1E-5;
const double reltolIn = 1E-7;
const double internalV = 623.0; // Same as that used in TAM model
const double internalFrac = 0.5; // Same as that used in TAM model

// The indices carried over in the reduced IL2 model
const array<size_t, 21> IL2_assoc = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 52}};

array<bool, 26> __active_species_IDX() {
	array<bool, 26> __active_species_IDX;
	fill(__active_species_IDX.begin(), __active_species_IDX.end(), false);

	__active_species_IDX[8] = true;
	__active_species_IDX[9] = true;
	__active_species_IDX[16] = true;
	__active_species_IDX[17] = true;
	__active_species_IDX[21] = true;
	__active_species_IDX[25] = true;

	return __active_species_IDX;
}

const array<bool, 26> activeV = __active_species_IDX();


template<typename DD>
void dy_dt(const DD * const y, const DD * const rxn, DD *dydt) {
	// Set the constant inputs
	DD IL2 = rxn[0];
	DD IL15 = rxn[1];
	DD IL7 = rxn[2];
	DD IL9 = rxn[3];
	DD kfwd = rxn[4];
	DD k5rev = rxn[5];
	DD k6rev = rxn[6];
	DD k15rev = rxn[7];
	DD k17rev = rxn[8];
	DD k18rev = rxn[9];
	DD k22rev = rxn[10];
	DD k23rev = rxn[11];
	DD k26rev = rxn[12];
	DD k27rev = rxn[13];
	DD k29rev = rxn[14];
	DD k30rev = rxn[15];
	DD k31rev = rxn[16];

	// IL2 in nM
	DD IL2Ra = y[0];
	DD IL2Rb = y[1];
	DD gc = y[2];
	DD IL2_IL2Ra = y[3];
	DD IL2_IL2Rb = y[4];
	DD IL2_gc = y[5];
	DD IL2_IL2Ra_IL2Rb = y[6];
	DD IL2_IL2Ra_gc = y[7];
	DD IL2_IL2Rb_gc = y[8];
	DD IL2_IL2Ra_IL2Rb_gc = y[9];
	
	// IL15 in nM
	DD IL15Ra = y[10];
	DD IL15_IL15Ra = y[11];
	DD IL15_IL2Rb = y[12];
	DD IL15_gc = y[13];
	DD IL15_IL15Ra_IL2Rb = y[14];
	DD IL15_IL15Ra_gc = y[15];
	DD IL15_IL2Rb_gc = y[16];
	DD IL15_IL15Ra_IL2Rb_gc = y[17];
	
	// IL7, IL9 in nM
	DD IL7Ra = y[18];
	DD IL7Ra_IL7 = y[19];
	DD gc_IL7 = y[20];
	DD IL7Ra_gc_IL7 = y[21];
	DD IL9R = y[22];
	DD IL9R_IL9 = y[23];
	DD gc_IL9 = y[24];
	DD IL9R_gc_IL9 = y[25];

	// These are probably measured in the literature
	DD kfbnd = 0.01; // Assuming on rate of 10^7 M-1 sec-1
	DD k1rev = kfbnd * 10; // doi:10.1016/j.jmb.2004.04.038, 10 nM

	DD k2rev = kfbnd * 144; // doi:10.1016/j.jmb.2004.04.038, 144 nM
	DD k3fwd = kfbnd / 10.0; // Very weak, > 50 uM. Voss, et al (1993). PNAS. 90, 2428–2432.
	DD k3rev = 50000 * k3fwd;
	DD k10rev = 12.0 * k5rev / 1.5; // doi:10.1016/j.jmb.2004.04.038
	DD k11rev = 63.0 * k5rev / 1.5; // doi:10.1016/j.jmb.2004.04.038
	
	// Literature values for k values for IL-15
	DD k13rev = kfbnd * 0.065; // based on the multiple papers suggesting 30-100 pM
	DD k14rev = kfbnd * 438; // doi:10.1038/ni.2449, 438 nM
	
	// Literature values for IL-7
	DD k25rev = kfbnd * 59; // DOI:10.1111/j.1600-065X.2012.01160.x, 59 nM
	
	// To satisfy detailed balance these relationships should hold
	// _Based on initial assembly steps
	DD k4rev = kfbnd * k6rev * k3rev / k1rev / k3fwd;
	DD k7rev = k3fwd * k2rev * k5rev / kfbnd / k3rev;
	DD k12rev = k1rev * k11rev / k2rev;
	// _Based on formation of full complex
	DD k9rev = k2rev * k10rev * k12rev / kfbnd / k3rev / k6rev * k3fwd;
	DD k8rev = k2rev * k10rev * k12rev / kfbnd / k7rev / k3rev * k3fwd;

	// IL15
	// To satisfy detailed balance these relationships should hold
	// _Based on initial assembly steps
	DD k16rev = kfwd * k18rev * k15rev / k13rev / kfbnd;
	DD k19rev = kfwd * k14rev * k17rev / kfbnd / k15rev;
	DD k24rev = k13rev * k23rev / k14rev;
	// _Based on formation of full complex

	DD k21rev = k14rev * k22rev * k24rev / kfwd / k15rev / k18rev * kfbnd;
	DD k20rev = k14rev * k22rev * k24rev / k19rev / k15rev;

	// _One detailed balance IL7/9 loop
	DD k32rev = k29rev * k31rev / k30rev;
	DD k28rev = k25rev * k27rev / k26rev;
	
	// IL2
	dydt[0] = -kfbnd * IL2Ra * IL2 + k1rev * IL2_IL2Ra - kfwd * IL2Ra * IL2_gc + k6rev * IL2_IL2Ra_gc - kfwd * IL2Ra * IL2_IL2Rb_gc + k8rev * IL2_IL2Ra_IL2Rb_gc - kfwd * IL2Ra * IL2_IL2Rb + k12rev * IL2_IL2Ra_IL2Rb;
	dydt[1] = -kfbnd * IL2Rb * IL2 + k2rev * IL2_IL2Rb - kfwd * IL2Rb * IL2_gc + k7rev * IL2_IL2Rb_gc - kfwd * IL2Rb * IL2_IL2Ra_gc + k9rev * IL2_IL2Ra_IL2Rb_gc - kfwd * IL2Rb * IL2_IL2Ra + k11rev * IL2_IL2Ra_IL2Rb;
	dydt[2] = -k3fwd * IL2 * gc + k3rev * IL2_gc - kfwd * IL2_IL2Rb * gc + k5rev * IL2_IL2Rb_gc - kfwd * IL2_IL2Ra * gc + k4rev * IL2_IL2Ra_gc - kfwd * IL2_IL2Ra_IL2Rb * gc + k10rev * IL2_IL2Ra_IL2Rb_gc;
	dydt[3] = -kfwd * IL2_IL2Ra * IL2Rb + k11rev * IL2_IL2Ra_IL2Rb - kfwd * IL2_IL2Ra * gc + k4rev * IL2_IL2Ra_gc + kfbnd * IL2 * IL2Ra - k1rev * IL2_IL2Ra;
	dydt[4] = -kfwd * IL2_IL2Rb * IL2Ra + k12rev * IL2_IL2Ra_IL2Rb - kfwd * IL2_IL2Rb * gc + k5rev * IL2_IL2Rb_gc + kfbnd * IL2 * IL2Rb - k2rev * IL2_IL2Rb;
	dydt[5] = -kfwd *IL2_gc * IL2Ra + k6rev * IL2_IL2Ra_gc - kfwd * IL2_gc * IL2Rb + k7rev * IL2_IL2Rb_gc + k3fwd * IL2 * gc - k3rev * IL2_gc;
	dydt[6] = -kfwd * IL2_IL2Ra_IL2Rb * gc + k10rev * IL2_IL2Ra_IL2Rb_gc + kfwd * IL2_IL2Ra * IL2Rb - k11rev * IL2_IL2Ra_IL2Rb + kfwd * IL2_IL2Rb * IL2Ra - k12rev * IL2_IL2Ra_IL2Rb;
	dydt[7] = -kfwd * IL2_IL2Ra_gc * IL2Rb + k9rev * IL2_IL2Ra_IL2Rb_gc + kfwd * IL2_IL2Ra * gc - k4rev * IL2_IL2Ra_gc + kfwd * IL2_gc * IL2Ra - k6rev * IL2_IL2Ra_gc;
	dydt[8] = -kfwd * IL2_IL2Rb_gc * IL2Ra + k8rev * IL2_IL2Ra_IL2Rb_gc + kfwd * gc * IL2_IL2Rb - k5rev * IL2_IL2Rb_gc + kfwd * IL2_gc * IL2Rb - k7rev * IL2_IL2Rb_gc;
	dydt[9] = kfwd * IL2_IL2Rb_gc * IL2Ra - k8rev * IL2_IL2Ra_IL2Rb_gc + kfwd * IL2_IL2Ra_gc * IL2Rb - k9rev * IL2_IL2Ra_IL2Rb_gc + kfwd * IL2_IL2Ra_IL2Rb * gc - k10rev * IL2_IL2Ra_IL2Rb_gc;

	// IL15
	dydt[10] = -kfbnd * IL15Ra * IL15 + k13rev * IL15_IL15Ra - kfbnd * IL15Ra * IL15_gc + k18rev * IL15_IL15Ra_gc - kfwd * IL15Ra * IL15_IL2Rb_gc + k20rev * IL15_IL15Ra_IL2Rb_gc - kfwd * IL15Ra * IL15_IL2Rb + k24rev * IL15_IL15Ra_IL2Rb;
	dydt[11] = -kfwd * IL15_IL15Ra * IL2Rb + k23rev * IL15_IL15Ra_IL2Rb - kfwd * IL15_IL15Ra * gc + k16rev * IL15_IL15Ra_gc + kfbnd * IL15 * IL15Ra - k13rev * IL15_IL15Ra;
	dydt[12] = -kfwd * IL15_IL2Rb * IL15Ra + k24rev * IL15_IL15Ra_IL2Rb - kfbnd * IL15_IL2Rb * gc + k17rev * IL15_IL2Rb_gc + kfbnd * IL15 * IL2Rb - k14rev * IL15_IL2Rb;
	dydt[13] = -kfbnd * IL15_gc * IL15Ra + k18rev * IL15_IL15Ra_gc - kfwd * IL15_gc * IL2Rb + k19rev * IL15_IL2Rb_gc + kfbnd * IL15 * gc - k15rev * IL15_gc;
	dydt[14] = -kfwd * IL15_IL15Ra_IL2Rb * gc + k22rev * IL15_IL15Ra_IL2Rb_gc + kfwd * IL15_IL15Ra * IL2Rb - k23rev * IL15_IL15Ra_IL2Rb + kfwd * IL15_IL2Rb * IL15Ra - k24rev * IL15_IL15Ra_IL2Rb;
	dydt[15] = -kfwd * IL15_IL15Ra_gc * IL2Rb + k21rev * IL15_IL15Ra_IL2Rb_gc + kfwd * IL15_IL15Ra * gc - k16rev * IL15_IL15Ra_gc + kfbnd * IL15_gc * IL15Ra - k18rev * IL15_IL15Ra_gc;
	dydt[16] = -kfwd * IL15_IL2Rb_gc * IL15Ra + k20rev * IL15_IL15Ra_IL2Rb_gc + kfbnd * gc * IL15_IL2Rb - k17rev * IL15_IL2Rb_gc + kfwd * IL15_gc * IL2Rb - k19rev * IL15_IL2Rb_gc;
	dydt[17] =  kfwd * IL15_IL2Rb_gc * IL15Ra - k20rev * IL15_IL15Ra_IL2Rb_gc + kfwd * IL15_IL15Ra_gc * IL2Rb - k21rev * IL15_IL15Ra_IL2Rb_gc + kfwd * IL15_IL15Ra_IL2Rb * gc - k22rev * IL15_IL15Ra_IL2Rb_gc;
	
	dydt[1] = dydt[1] - kfbnd * IL2Rb * IL15 + k14rev * IL15_IL2Rb - kfwd * IL2Rb * IL15_gc + k19rev * IL15_IL2Rb_gc - kfwd * IL2Rb * IL15_IL15Ra_gc + k21rev * IL15_IL15Ra_IL2Rb_gc - kfwd * IL2Rb * IL15_IL15Ra + k23rev * IL15_IL15Ra_IL2Rb;
	dydt[2] = dydt[2] - kfbnd * IL15 * gc + k15rev * IL15_gc - kfbnd * IL15_IL2Rb * gc + k17rev * IL15_IL2Rb_gc - kfwd * IL15_IL15Ra * gc + k16rev * IL15_IL15Ra_gc - kfwd * IL15_IL15Ra_IL2Rb * gc + k22rev * IL15_IL15Ra_IL2Rb_gc;
	
	// IL7
	dydt[2] = dydt[2] - kfbnd * IL7 * gc + k26rev * gc_IL7 - kfwd * gc * IL7Ra_IL7 + k27rev * IL7Ra_gc_IL7;
	dydt[18] = -kfbnd * IL7Ra * IL7 + k25rev * IL7Ra_IL7 - kfwd * IL7Ra * gc_IL7 + k28rev * IL7Ra_gc_IL7;
	dydt[19] = kfbnd * IL7Ra * IL7 - k25rev * IL7Ra_IL7 - kfwd * gc * IL7Ra_IL7 + k27rev * IL7Ra_gc_IL7;
	dydt[20] = -kfwd * IL7Ra * gc_IL7 + k28rev * IL7Ra_gc_IL7 + kfbnd * IL7 * gc - k26rev * gc_IL7;
	dydt[21] = kfwd * IL7Ra * gc_IL7 - k28rev * IL7Ra_gc_IL7 + kfwd * gc * IL7Ra_IL7 - k27rev * IL7Ra_gc_IL7;

	// IL9
	dydt[2] = dydt[2] - kfbnd * IL9 * gc + k30rev * gc_IL9 - kfwd * gc * IL9R_IL9 + k31rev * IL9R_gc_IL9;
	dydt[22] = -kfbnd * IL9R * IL9 + k29rev * IL9R_IL9 - kfwd * IL9R * gc_IL9 + k32rev * IL9R_gc_IL9;
	dydt[23] = kfbnd * IL9R * IL9 - k29rev * IL9R_IL9 - kfwd * gc * IL9R_IL9 + k31rev * IL9R_gc_IL9;
	dydt[24] = -kfwd * IL9R * gc_IL9 + k32rev * IL9R_gc_IL9 + kfbnd * IL9 * gc - k30rev * gc_IL9;
	dydt[25] = kfwd * IL9R * gc_IL9 - k32rev * IL9R_gc_IL9 + kfwd * gc * IL9R_IL9 - k31rev * IL9R_gc_IL9;
}


extern "C" void dydt_C(double *y_in, double t, double *dydt_out, double *rxn_in) {
	dy_dt(y_in, rxn_in, dydt_out);
}


template<typename DD>
void findLigConsume(DD *dydt) {
	// Calculate the ligand consumption.
	dydt[52] -= std::accumulate(dydt+3, dydt+10, 0) / internalV;
	dydt[53] -= std::accumulate(dydt+11, dydt+18, 0) / internalV;
	dydt[54] -= std::accumulate(dydt+19, dydt+22, 0) / internalV;
	dydt[55] -= std::accumulate(dydt+23, dydt+26, 0) / internalV;
}


template<typename DD>
void trafficking(const DD * const y, const DD * const tfR, DD *dydt) {
	// Implement trafficking.

	// Set the rates
	DD endo = tfR[0];
	DD activeEndo = tfR[1];
	DD sortF = tfR[2];
	DD kRec = tfR[3];
	DD kDeg = tfR[4];

	size_t halfL = activeV.size();

	// Actually calculate the trafficking
	for (size_t ii = 0; ii < halfL; ii++) {
		if (activeV[ii]) {
			dydt[ii] += -y[ii]*(endo + activeEndo); // Endocytosis
			dydt[ii+halfL] += y[ii]*(endo + activeEndo)/internalFrac - kDeg*y[ii+halfL]; // Endocytosis, degradation
		} else {
			dydt[ii] += -y[ii]*endo + kRec*(1.0-sortF)*y[ii+halfL]*internalFrac; // Endocytosis, recycling
			dydt[ii+halfL] += y[ii]*endo/internalFrac - kRec*(1.0-sortF)*y[ii+halfL] - (kDeg*sortF)*y[ii+halfL]; // Endocytosis, recycling, degradation
		}
	}

	// Expression: IL2Ra, IL2Rb, gc, IL15Ra, IL7Ra, IL9R
	dydt[0] += tfR[5];
	dydt[1] += tfR[6];
	dydt[2] += tfR[7];
	dydt[10] += tfR[8];
	dydt[18] += tfR[9];
	dydt[22] += tfR[10];

	// Degradation does lead to some clearance of ligand in the endosome
	for (size_t ii = 0; ii < 4; ii++) {
		dydt[52 + ii] -= dydt[52 + ii] * kDeg;
	}
}


template<typename DD>
void fullModel(const DD * const y, const DD * const r, const DD * const tfR, DD *dydt) {
	// Implement full model.
	fill(dydt, dydt + 56, 0.0);

	// Calculate endosomal reactions
	array<DD, 17> rr;
	copy_n(r, rr.size(), rr.begin());
	copy_n(y + 52, 4, rr.begin());

	// Calculate cell surface and endosomal reactions
	dy_dt(y, r, dydt);
	dy_dt(y + 26, rr.data(), dydt + 26);

	// Handle trafficking
	trafficking(y, tfR, dydt);

	// Handle endosomal ligand balance.
	findLigConsume(dydt);
}


int fullModelCVode (const double, const N_Vector xx, N_Vector dxxdt, void *user_data) {
	ratesS *rIn = static_cast<ratesS *>(user_data);

	array<double, 56> xxArr;

	// Get the data in the right form
	if (NV_LENGTH_S(xx) == xxArr.size()) { // If we're using the full model
		fullModel(xxArr.data(), rIn->rxn, rIn->trafRates, NV_DATA_S(dxxdt));
	} else if (NV_LENGTH_S(xx) == IL2_assoc.size()) { // If it looks like we're using the IL2 model
		fill(xxArr.begin(), xxArr.end(), 0.0);

		for (size_t ii = 0; ii < IL2_assoc.size(); ii++) {
			xxArr[IL2_assoc[ii]] = NV_Ith_S(xx, ii);
		}

		array<double, 56> dydt;

		fullModel(xxArr.data(), rIn->rxn, rIn->trafRates, dydt.data());

		for (size_t ii = 0; ii < IL2_assoc.size(); ii++) {
			NV_Ith_S(dxxdt, ii) = dydt[IL2_assoc[ii]];
		}
	} else {
		throw runtime_error(string("Failed to find the right wrapper."));
	}

	return 0;
}


extern "C" void fullModel_C(const double * const y_in, double t, double *dydt_out, double *rxn_in, double *tfr_in) {
	fullModel(y_in, rxn_in, tfr_in, dydt_out);
}


array<double, 56> solveAutocrine(const double * const trafRates) {
	array<double, 56> y0;
	fill(y0.begin(), y0.end(), 0.0);

	array<size_t, 6> recIDX = {{0, 1, 2, 10, 18, 22}};

	double internalFrac = 0.5; // Same as that used in TAM model

	// Expand out trafficking terms
	double endo = trafRates[0];
	double sortF = trafRates[2];
	double kRec = trafRates[3]*(1-sortF);
	double kDeg = trafRates[4]*sortF;

	// Assuming no autocrine ligand, so can solve steady state
	// Add the species
	for (size_t ii = 0; ii < recIDX.size(); ii++) {
		y0[recIDX[ii] + 26] = trafRates[5 + ii] / kDeg / internalFrac;
		y0[recIDX[ii]] = (trafRates[5 + ii] + kRec*y0[recIDX[ii] + 26]*internalFrac)/endo;
	}

	return y0;
}


static void errorHandler(int error_code, const char *module, const char *function, char *msg, void *user_data) {
	if (error_code == CV_WARNING) return;

	std::cout << "Internal CVode error in " << function << std::endl;
	std::cout << msg << std::endl;
	std::cout << "In module: " << module << std::endl;
	std::cout << "Error code: " << error_code << std::endl;
}

void* solver_setup(N_Vector init, void * params) {
	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration */
	void *cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

	if (cvode_mem == NULL) {
		CVodeFree(&cvode_mem);
		throw runtime_error(string("Error calling CVodeCreate in solver_setup."));
	}

	auto checkReturn = [&cvode_mem] (int retVal, string name) {
		if (retVal < 0) {
			CVodeFree(&cvode_mem);
			throw runtime_error(string("Error calling ") + name + string(" in solver_setup."));
		}
	};
	
	CVodeSetErrHandlerFn(cvode_mem, &errorHandler, params);

	/* Call CVodeInit to initialize the integrator memory and specify the
	 * user's right hand side function in y'=f(t,y), the inital time T0, and
	 * the initial dependent variable vector y. */
	checkReturn(CVodeInit(cvode_mem, fullModelCVode, 0.0, init), "CVodeInit");
	
	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	checkReturn(CVodeSStolerances(cvode_mem, reltolIn, abstolIn), "CVodeSVtolerances");

	SUNMatrix A = SUNDenseMatrix((int) NV_LENGTH_S(init), (int) NV_LENGTH_S(init));

	SUNLinearSolver LS = SUNDenseLinearSolver(init, A);
	
	// Call CVDense to specify the CVDENSE dense linear solver
	checkReturn(CVDlsSetLinearSolver(cvode_mem, LS, A), "CVDlsSetLinearSolver");
	
	// Pass along the parameter structure to the differential equations
	checkReturn(CVodeSetUserData(cvode_mem, params), "CVodeSetUserData");

	CVodeSetMaxNumSteps(cvode_mem, 2000000);
	
	return cvode_mem;
}


void sensi_setup(void *cvode_mem, N_Vector init, ratesS *params) {
	// Setup all the initial sensitivities
	N_Vector *uS = N_VCloneVectorArray(plen, init);
	for(size_t is = 0; is < plen; is++) N_VConst(0.0, uS[is]);

	auto checkReturn = [&cvode_mem] (int retVal, string name) {
		if (retVal < 0) {
			CVodeFree(&cvode_mem);
			throw runtime_error(string("Error calling ") + name + string(" in solver_setup."));
		}
	};

	checkReturn(CVodeSensInit1(cvode_mem, plen, CV_SIMULTANEOUS, NULL, uS), "CVodeSensInit1");
	N_VDestroyVectorArray_Serial(uS, plen);

	array<double, plen> abstol;
	fill(abstol.begin(), abstol.end(), 1.0E-2);
	checkReturn(CVodeSensSStolerances(cvode_mem, 1.0E-2, abstol.data()), "CVodeSensSStolerances");

	checkReturn(CVodeSetSensErrCon(cvode_mem, true), "CVodeSetSensErrCon");

	checkReturn(CVodeSetSensDQMethod(cvode_mem, CV_CENTERED, 0.0), "CVodeSetSensDQMethod");

	checkReturn(CVodeSetSensParams(cvode_mem, params->p.data(), NULL, NULL), "CVodeSetSensParams");
}


extern "C" int runCkine (double *tps, size_t ntps, double *out, double *rxnRatesIn, double *trafRatesIn) {
	ratesS rattes;
	rattes.rxn = rattes.p.data();
	rattes.trafRates = rattes.p.data() + 17;

	copy_n(rxnRatesIn, 17, rattes.rxn);
	copy_n(trafRatesIn, 11, rattes.trafRates);

	array<double, 56> y0 = solveAutocrine(rattes.trafRates);
	N_Vector state;

	// Fill output values with 0's
	fill(out, out + ntps*y0.size(), 0.0);

	// Can we use the reduced IL2 only model
	if (rattes.trafRates[8] + rattes.trafRates[9] + rattes.trafRates[10] == 0.0) {
		state = N_VNew_Serial((long) IL2_assoc.size());

		for (size_t ii = 0; ii < IL2_assoc.size(); ii++) {
			NV_Ith_S(state, ii) = y0[IL2_assoc[ii]];
		}
	} else { // Just the full model
		state = N_VMake_Serial((long) y0.size(), y0.data());
	}

	void *cvode_mem = solver_setup(state, (void *) &rattes);

	sensi_setup(cvode_mem, state, &rattes);

	double tret = 0;

	for (size_t itps = 0; itps < ntps; itps++) {
		if (tps[itps] < tret) {
			std::cout << "Can't go backwards." << std::endl;
			N_VDestroy_Serial(state);
			CVodeFree(&cvode_mem);
			return -1;
		}

		int returnVal = CVode(cvode_mem, tps[itps], state, &tret, CV_NORMAL);
		
		if (returnVal < 0) {
			std::cout << "CVode error in CVode. Code: " << returnVal << std::endl;
			N_VDestroy_Serial(state);
			CVodeFree(&cvode_mem);
			return returnVal;
		}

		// Copy out result
		if (NV_LENGTH_S(state) == y0.size()) {
			copy_n(NV_DATA_S(state), y0.size(), out + y0.size()*itps);
		} else { // If we're dealing with a reduced model
			for (size_t ii = 0; ii < IL2_assoc.size(); ii++) {
				out[y0.size()*itps + IL2_assoc[ii]] = NV_Ith_S(state, ii);
			}
		}
	}

	N_VDestroy_Serial(state);
	CVodeFree(&cvode_mem);
	return 0;
}
