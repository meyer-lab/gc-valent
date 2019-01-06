#include <algorithm>
#include <cstdio>
#include <numeric>
#include <array>
#include <thread>
#include <vector>
#include <list>
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode.h>            /* prototypes for CVODE fcts., consts. */
#include <string>
#include <sundials/sundials_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <cvode/cvode_direct.h>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cppad/cppad.hpp>
#include "model.hpp"
#include "reaction.hpp"
#include "jacobian.hpp"
#include "thread_pool.hpp"

using std::array;
using std::copy;
using std::vector;
using std::fill;
using std::string;
using std::endl;
using std::cout;

static void errorHandler(int, const char *, const char *, char *, void *);
int ewt(N_Vector, N_Vector, void *);
int Jac(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector, N_Vector, N_Vector);
int fullModelCVode (const double, const N_Vector, N_Vector, void *);

std::mutex print_mutex; // mutex to prevent threads printing on top of each other


extern "C" void dydt_C(double *y_in, double, double *dydt_out, double *rxn_in) {
	std::vector<double> v(rxn_in, rxn_in + Nparams);
	ratesS<double> r(v);

	dy_dt(y_in, &r.surface, dydt_out, r.ILs.data());
}


extern "C" void fullModel_C(const double * const y_in, double, double *dydt_out, double *rxn_in) {
	std::vector<double> v(rxn_in, rxn_in + Nparams);
	ratesS<double> r(v);

	fullModel(y_in, &r, dydt_out);
}


class solver {
public:
	void *cvode_mem;
	SUNLinearSolver LS;
	N_Vector state, qB;
	SUNMatrix A;
	bool sensi;
	int ncheck;
	double tret;
	vector<double> params;
	array<double, Nspecies> activities;

	void commonSetup(vector<double> paramsIn) {
		tret = 0.0;
		params = paramsIn;

		// Setup state variable by solving for autocrine
		ratesS<double> rattes(params);
		array<double, Nspecies> y0 = solveAutocrine(&rattes);
		state = N_VNew_Serial(static_cast<long>(Nspecies));
		std::copy_n(y0.data(), Nspecies, NV_DATA_S(state));

		/* Call CVodeCreate to create the solver memory and specify the
		 * Backward Differentiation Formula and the use of a Newton iteration */
		cvode_mem = CVodeCreate(CV_BDF);
		if (cvode_mem == nullptr) {
			throw std::runtime_error(string("Error calling CVodeCreate in solver_setup."));
		}
		
		CVodeSetErrHandlerFn(cvode_mem, &errorHandler, static_cast<void *>(this));

		// Pass along the parameter structure to the differential equations
		if (CVodeSetUserData(cvode_mem, static_cast<void *>(this)) < 0) {
			throw std::runtime_error(string("Error calling CVodeSetUserData in solver_setup."));
		}

		/* Call CVodeInit to initialize the integrator memory and specify the
		 * user's right hand side function in y'=f(t,y), the inital time T0, and
		 * the initial dependent variable vector y. */
		if (CVodeInit(cvode_mem, fullModelCVode, 0.0, state) < 0) {
			throw std::runtime_error(string("Error calling CVodeInit in solver_setup."));
		}
		
		// Call CVodeWFtolerances to specify the tolerances
		if (CVodeWFtolerances(cvode_mem, ewt) < 0) {
			throw std::runtime_error(string("Error calling CVodeWFtolerances in solver_setup."));
		}

		A = SUNDenseMatrix(NV_LENGTH_S(state), NV_LENGTH_S(state));
		LS = SUNDenseLinearSolver(state, A);
		
		// Call CVDense to specify the CVDENSE dense linear solver
		if (CVDlsSetLinearSolver(cvode_mem, LS, A) < 0) {
			throw std::runtime_error(string("Error calling CVDlsSetLinearSolver in solver_setup."));
		}

		CVDlsSetJacFn(cvode_mem, Jac);

		CVodeSetMaxNumSteps(cvode_mem, 800000);

		// Call CVodeSetConstraints to initialize constraints
		N_Vector constraints = N_VNew_Serial(static_cast<long>(Nspecies));
		N_VConst(1.0, constraints); // all 1's for nonnegative solution values
		if (CVodeSetConstraints(cvode_mem, constraints) < 0) {
			throw std::runtime_error(string("Error calling CVodeSetConstraints in solver_setup."));
		}
		N_VDestroy(constraints);
	}


	solver(vector<double> paramsIn) {
		sensi = false;
		commonSetup(paramsIn);
	}

	solver(vector<double> paramsIn, array<double, Nspecies> actIn) {
		sensi = true;
		std::copy(actIn.begin(), actIn.end(), activities.begin());
		commonSetup(paramsIn);


	}

	int CVodeRun(double endT) {
		int returnVal;

		if (sensi) {
			returnVal = CVodeF(cvode_mem, endT, state, &tret, CV_NORMAL, &ncheck);
		} else {
			returnVal = CVode(cvode_mem, endT, state, &tret, CV_NORMAL);
		}
		
		if (returnVal < 0) std::cout << "CVode error in CVode. Code: " << returnVal << std::endl;

		return returnVal;
	}

	ratesS<double> getRates() {
		return ratesS<double>(params);
	}

	~solver() {
		if (sensi) {
			CVodeSensFree(cvode_mem);
			N_VDestroy_Serial(qB);
		}

		N_VDestroy_Serial(state);
		CVodeFree(&cvode_mem);
		SUNLinSolFree(LS);
		SUNMatDestroy(A);
	}
};


// fQ routine. Compute fQ(t,y)
static int fQ(double, N_Vector y, N_Vector qdot, void *ehdata) {
	solver *sMem = static_cast<solver *>(ehdata);

	NV_Ith_S(qdot, 0) = std::inner_product(sMem->activities.begin(), sMem->activities.end(), NV_DATA_S(y), 0.0);

	return 0;
}


// fB routine. Compute fB(t,y,yB). 
static int fB(double, N_Vector y, N_Vector yB, N_Vector yBdot, void *user_dataB) {
	solver *sMem = static_cast<solver *>(user_dataB);
	ratesS<double> rattes = sMem->getRates();

	std::copy(sMem->activities.begin(), sMem->activities.end(), NV_DATA_S(yBdot));

	Eigen::Map<Eigen::Matrix<double, Nspecies, 1>> yBv(NV_DATA_S(yB), Nspecies);
	Eigen::Map<Eigen::Matrix<double, Nspecies, 1>> yBdotv(NV_DATA_S(yBdot), Nspecies);

	JacMat jac;

	// Actually get the Jacobian
	fullJacobian(NV_DATA_S(y), &rattes, jac);

	// TODO: Might need to transpose

	yBdotv += jac*yBv;

	return 0;
}


// fQB routine. Compute integrand for quadratures 
static int fQB(double, N_Vector y, N_Vector yB, N_Vector qBdot, void *user_dataB) {
	using CppAD::AD;
	using Eigen::Map;
	using Eigen::Dynamic;
	using Eigen::Matrix;
	solver *sMem = static_cast<solver *>(user_dataB);

	size_t Np = sMem->params.size();

	Map<Matrix<double, Nspecies, 1>> yBdotv(NV_DATA_S(qBdot), Nspecies);

	vector<AD<double>> X(sMem->params.size());
	vector<double> xx(sMem->params.size());

	std::copy(sMem->params.begin(), sMem->params.end(), X.begin());

	CppAD::Independent(X);

	vector<AD<double>> dydt(Nspecies);

	ratesS<AD<double>> rattes = ratesS<AD<double>>(X);

	// Get the data in the right form
	fullModel(NV_DATA_S(y), &rattes, dydt.data());

	CppAD::ADFun<double> f(X, dydt);

	// Copy in values to actually use as x
	std::copy(sMem->params.begin(), sMem->params.end(), xx.begin());



	vector<double> jac(Nspecies * sMem->params.size());
	jac = f.Jacobian(xx); // Jacobian for operation sequence
	Map<Matrix<double, Dynamic, Dynamic>> jac_M(jac.data(), Nspecies, Np);

	// Jac wrt params

	// Jac_p times y

	Eigen::Map<Eigen::Matrix<double, Nspecies, 1>> yBv(NV_DATA_S(yB), Nspecies);

	yBdotv = jac_M*yBv;

	return(0);
}


static void errorHandler(int error_code, const char *module, const char *function, char *msg, void *ehdata) {
	if (error_code == CV_WARNING) return;
	solver *sMem = static_cast<solver *>(ehdata);
	ratesS<double> ratt = sMem->getRates();

	std::lock_guard<std::mutex> lock(print_mutex);

	std::cout << "Internal CVode error in " << function << ", module: " << module << ", error code: " << error_code << std::endl;
	std::cout << msg << std::endl;
	std::cout << "Parameters: ";

	for (size_t ii = 0; ii < Nparams; ii++) {
		std::cout << sMem->params[ii] << "\t";
	}
	
	ratt.print();

	if (sMem->sensi)
		std::cout << "Sensitivity enabled." << std::endl;

	std::cout << std::endl << std::endl;
}


int ewt(N_Vector y, N_Vector w, void *ehdata) {
	solver *sMem = static_cast<solver *>(ehdata);

	double tolIn = 1E-7;

	if (sMem->sensi) tolIn = 1E-3;

	for (size_t i = 0; i < Nspecies; i++) {
		NV_Ith_S(w, i) = 1.0/(fabs(NV_Ith_S(y, i))*tolIn + tolIn);
	}

	return 0;
}


int Jac(realtype, N_Vector y, N_Vector, SUNMatrix J, void *user_data, N_Vector, N_Vector, N_Vector) {
	solver *sMem = static_cast<solver *>(user_data);
	ratesS<double> rattes = sMem->getRates();

	Eigen::Map<JacMat> jac(SM_DATA_D(J));

	// Actually get the Jacobian
	fullJacobian(NV_DATA_S(y), &rattes, jac);

	jac.transposeInPlace();

	return 0;
}


int fullModelCVode(const double, const N_Vector xx, N_Vector dxxdt, void *user_data) {
	solver *sMem = static_cast<solver *>(user_data);
	ratesS<double> rattes = sMem->getRates();

	// Get the data in the right form
	fullModel(NV_DATA_S(xx), &rattes, NV_DATA_S(dxxdt));

	return 0;
}


extern "C" int runCkine (double * const tps, const size_t ntps, double * const out, const double * const rxnRatesIn, const bool sensi, double * const sensiOut, bool IL2case) {
	size_t itps = 0;

	std::vector<double> v;

	if (IL2case) {
		v = std::vector<double>(rxnRatesIn, rxnRatesIn + NIL2params);
	} else {
		v = std::vector<double>(rxnRatesIn, rxnRatesIn + Nparams);
	}

	solver sMem(v);

	if (tps[0] < std::numeric_limits<double>::epsilon()) {
		std::copy_n(NV_DATA_S(sMem.state), Nspecies, out);

		itps = 1;
	}

	for (; itps < ntps; itps++) {
		if (tps[itps] < sMem.tret) {
			std::cout << "Can't go backwards." << std::endl;
			return -1;
		}
		
		if (sMem.CVodeRun(tps[itps]) < 0) return -1;

		// Copy out result
		std::copy_n(NV_DATA_S(sMem.state), Nspecies, out + Nspecies*itps);
	}

	return 0;
}


extern "C" int runCkinePretreat (const double pret, const double tt, double * const out, const double * const rxnRatesIn, const double * const postStim) {
	solver sMem(std::vector<double>(rxnRatesIn, rxnRatesIn + Nparams));

	if (sMem.CVodeRun(pret) < 0) return -1;

	std::copy_n(postStim, Nlig, sMem.params.begin()); // Copy in stimulation ligands

	CVodeReInit(sMem.cvode_mem, pret, sMem.state);
    
	if (tt > std::numeric_limits<double>::epsilon()) {
		if (sMem.CVodeRun(pret + tt) < 0) return -1;
	}

	// Copy out result
	std::copy_n(NV_DATA_S(sMem.state), Nspecies, out);

	return 0;
}

ThreadPool pool;

extern "C" int runCkineParallel (const double * const rxnRatesIn, double tp, size_t nDoses, bool sensi, double *out, double *sensiOut) {
	int retVal = 1000;
	std::list<std::future<int>> results;

	// Actually run the simulations
	for (size_t ii = 0; ii < nDoses; ii++)
		results.push_back(pool.enqueue(runCkine, &tp, 1, out + Nspecies*ii, rxnRatesIn + ii*Nparams, sensi, sensiOut + Nspecies*Nparams*ii, false));

	// Synchronize all threads
	for (std::future<int> &th:results) retVal = std::min(th.get(), retVal);

	// Get the worst case to return
	return retVal;
}


extern "C" void fullJacobian_C(double *y_in, double, double *dydt, double *rxn_in) {
	std::vector<double> v(rxn_in, rxn_in + Nparams);
	ratesS<double> r(v);

	Eigen::Map<JacMat> out(dydt);

	fullJacobian(y_in, &r, out);
}
