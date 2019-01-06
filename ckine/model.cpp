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
#include <adept.h>
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
static int fQ(double, N_Vector y, N_Vector qdot, void *ehdata);
static int fB(double, N_Vector y, N_Vector yB, N_Vector yBdot, void *user_dataB);
static int fQB(double, N_Vector y, N_Vector yB, N_Vector qBdot, void *user_dataB);

std::mutex print_mutex; // mutex to prevent threads printing on top of each other

typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> eigenV;
typedef Eigen::Map<Eigen::Matrix<double, Nspecies, 1>> eigenVC;


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

constexpr bool debug = false;


class solver {
public:
	void *cvode_mem;
	SUNLinearSolver LS, LSB;
	N_Vector state, qB, q, yB;
	SUNMatrix A, AB;
	bool sensi;
	int ncheck, indexB;
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

		CVodeSetMaxNumSteps(cvode_mem, 8000);

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

		// CVodeQuadInit to allocate initernal memory and initialize quadrature integration
		q = N_VNew_Serial(1);
		fQ(0.0, state, q, static_cast<void *>(this));
		if (CVodeQuadInit(cvode_mem, fQ, q) < 0) {
			throw std::runtime_error(string("Error calling CVodeQuadInit in solver_setup."));
		}

		// Whether or not the quadrature variables are to be used in the step size control mechanism
		if (CVodeSetQuadErrCon(cvode_mem, SUNTRUE) < 0) {
			throw std::runtime_error(string("Error calling CVodeSetQuadErrCon in solver_setup."));
		}

		// Specify scalar relative and absolute tolerances
		if (CVodeQuadSStolerances(cvode_mem, 1.0E-6, 1.0E-6) < 0) {
			throw std::runtime_error(string("Error calling CVodeQuadSStolerances in solver_setup."));
		}

		// CVodeAdjInit to update CVODES memory block by allocting the internal memory needed for backward integration
		constexpr int steps = 1; // no. of integration steps between two consecutive ckeckpoints
		if (CVodeAdjInit(cvode_mem, steps, CV_HERMITE) < 0) {
			throw std::runtime_error(string("Error calling CVodeAdjInit in solver_setup."));
		}
	}


	void backward (double TB1) {
		indexB = 1;
		yB = N_VNew_Serial(Nspecies); // Initialize yB
		qB = N_VNew_Serial(params.size()); // Initialize qB
		N_VConst(0.0, yB);
		N_VConst(0.0, qB);

		// CVodeCreateB to specify the solution method for the backward problem
		if (CVodeCreateB(cvode_mem, CV_BDF, &indexB) < 0) {
			throw std::runtime_error(string("Error calling CVodeCreateB in solver_setup."));
		}

		// Call CVodeInitB to allocate internal memory and initialize the backward problem
		if (CVodeInitB(cvode_mem, indexB, fB, TB1, yB) < 0) {
			throw std::runtime_error(string("Error calling CVodeInitB in solver_setup."));
		}

		// Set the scalar relative and absolute tolerances
		if (CVodeSStolerancesB(cvode_mem, indexB, 1.0E-3, 1.0E-3) < 0) {
			throw std::runtime_error(string("Error calling CVodeSStolerancesB in solver_setup."));
		}

		// Attach the user data for backward problem
		if (CVodeSetUserDataB(cvode_mem, indexB, static_cast<void *>(this)) < 0) {
			throw std::runtime_error(string("Error calling CVodeSetUserDataB in solver_setup."));
		}

		AB = SUNDenseMatrix(Nspecies, Nspecies);
		LSB = SUNLinSol_Dense(yB, AB);
		
		// Call CVDense to specify the CVDENSE dense linear solver
		if (CVodeSetLinearSolverB(cvode_mem, indexB, LSB, AB) < 0) {
			throw std::runtime_error(string("Error calling CVodeSetLinearSolverB in solver_setup."));
		}

		// Set the user-supplied Jacobian routine JacB
		if (CVodeSetJacFnB(cvode_mem, indexB, NULL) < 0) {
			throw std::runtime_error(string("Error calling CVodeSetJacFnB in solver_setup."));
		}

		// Allocate internal memory and initialize backward quadrature integration
		if (CVodeQuadInitB(cvode_mem, indexB, fQB, qB) < 0) {
		 	throw std::runtime_error(string("Error calling CVodeQuadInitB in solver_setup."));
		}

		// Whether or not the quadrature variables are to be used in the step size control
		if (CVodeSetQuadErrConB(cvode_mem, indexB, true) < 0) {
			throw std::runtime_error(string("Error calling CVodeSetQuadErrConB in solver_setup."));
		}

		// Specify the scalar relative and absolute tolerances for the backward problem
		if (CVodeQuadSStolerancesB(cvode_mem, indexB, 1.0E-6, 1.0E-6) < 0) {
			throw std::runtime_error(string("Error calling CVodeQuadSStolerancesB in solver_setup."));
		}
	}

	int CVodeRun(const double endT) {
		int returnVal;

		if (sensi) {
			returnVal = CVodeF(cvode_mem, endT, state, &tret, CV_NORMAL, &ncheck);
		} else {
			returnVal = CVode(cvode_mem, endT, state, &tret, CV_NORMAL);
		}

		if (returnVal >= 0 && debug) {
			long nst;
			CVodeGetNumSteps(cvode_mem, &nst);
			cout << "Number of steps: " << nst << std::endl;
			cout << "Final time: " << tret << std::endl;
		}
		
		if (returnVal < 0) cout << "CVode error in CVode. Code: " << returnVal << std::endl;

		return returnVal;
	}

	ratesS<double> getRates() {
		return ratesS<double>(params);
	}

	~solver() {
		if (sensi) {
			CVodeSensFree(cvode_mem);
			N_VDestroy_Serial(qB);
			N_VDestroy_Serial(q);
			N_VDestroy_Serial(yB);
			SUNLinSolFree(LSB);
			SUNMatDestroy(AB);
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

	eigenVC yBv(NV_DATA_S(yB), Nspecies);
	eigenVC yBdotv(NV_DATA_S(yBdot), Nspecies);

	JacMat jac;

	// Actually get the Jacobian
	fullJacobian(NV_DATA_S(y), &rattes, jac);

	yBdotv = -yBdotv - jac.transpose()*yBv;

	return 0;
}


// fQB routine. Compute integrand for quadratures 
static int fQB(double, N_Vector y, N_Vector yB, N_Vector qBdot, void *user_dataB) {
	using adept::adouble;
	solver *sMem = static_cast<solver *>(user_dataB);

	size_t Np = sMem->params.size();

	// Wrap the vectors we'll use
	eigenV yBdotv(NV_DATA_S(qBdot), NV_LENGTH_S(qBdot));
	eigenV yBv(NV_DATA_S(yB), NV_LENGTH_S(yB));

	adept::Stack stack;

	vector<adouble> X(Np);
	adept::set_values(&X[0], Np, sMem->params.data());

	stack.new_recording();

	vector<adouble> dydt(Nspecies);

	ratesS<adouble> rattes = ratesS<adouble>(X);

	// Get the data in the right form
	fullModel(NV_DATA_S(y), &rattes, dydt.data());

	stack.independent(&X[0], Np);
	stack.dependent(&dydt[0], Nspecies);

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> jac_M(Nspecies, Np);

	stack.jacobian(jac_M.data());

	yBdotv = jac_M.transpose()*yBv;

	return(0);
}


static void errorHandler(int error_code, const char *module, const char *function, char *msg, void *ehdata) {
	if (error_code == CV_WARNING) return;
	solver *sMem = static_cast<solver *>(ehdata);
	ratesS<double> ratt = sMem->getRates();

	std::lock_guard<std::mutex> lock(print_mutex);

	cout << "Internal CVode error in " << function << ", module: " << module << ", error code: " << error_code << std::endl;
	cout << msg << std::endl;
	cout << "Parameters: ";

	for (size_t ii = 0; ii < Nparams; ii++) {
		cout << sMem->params[ii] << "\t";
	}
	
	ratt.print();

	if (sMem->sensi)
		cout << "Sensitivity enabled." << std::endl;

	cout << std::endl << std::endl;
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


extern "C" int runCkine (double * const tps, const size_t ntps, double * const out, const double * const rxnRatesIn, bool IL2case) {
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
			cout << "Can't go backwards." << std::endl;
			return -1;
		}
		
		if (sMem.CVodeRun(tps[itps]) < 0) return -1;

		// Copy out result
		std::copy_n(NV_DATA_S(sMem.state), Nspecies, out + Nspecies*itps);
	}

	return 0;
}


extern "C" int runCkineS (const double * const tps, const size_t ntps, double * const out, double * const Sout, double * const actV, const double * const rxnRatesIn, bool IL2case) {
	size_t itps = 0;

	std::vector<double> v;
	std::array<double, Nspecies> actVv;
	std::copy_n(actV, Nspecies, actVv.begin());

	v = std::vector<double>(rxnRatesIn, rxnRatesIn + Nparams);

	solver sMem(v, actVv);

	if (tps[0] < std::numeric_limits<double>::epsilon()) {
		out[0] = std::inner_product(sMem.activities.begin(), sMem.activities.end(), NV_DATA_S(sMem.state), 0.0);

		itps = 1;
	}

	for (; itps < ntps; itps++) {
		if (tps[itps] < sMem.tret) {
			cout << "Can't go backwards." << std::endl;
			return -1;
		}
		
		if (sMem.CVodeRun(tps[itps]) < 0) return -1;

		// Copy out result
		if (CVodeGetQuad(sMem.cvode_mem, &sMem.tret, sMem.q) < 0) return -100;
		out[itps] = NV_Ith_S(sMem.q, 0);
	}

	if (sMem.CVodeRun(tps[ntps-1] + 10.0) < 0) return -1;
	sMem.backward(tps[ntps-1] + 10.0);

	for (int bitps = ntps - 1; bitps >= 0; bitps--) {
		if (CVodeB(sMem.cvode_mem, tps[bitps], CV_NORMAL) < 0) {
			cout << "CVodeB error." << std::endl;
		}

		if (CVodeGetQuadB(sMem.cvode_mem, sMem.indexB, &sMem.tret, sMem.qB) < 0) {
			cout << "CVodeGetQuadB error." << std::endl;
		}

		std::copy_n(NV_DATA_S(sMem.qB), NV_LENGTH_S(sMem.qB), Sout + bitps*NV_LENGTH_S(sMem.qB));
	}

	if (CVodeB(sMem.cvode_mem, 0.0, CV_NORMAL) < 0) {
			cout << "CVodeB error at 0." << std::endl;
	}

	if (CVodeGetQuadB(sMem.cvode_mem, sMem.indexB, &sMem.tret, sMem.qB) < 0) {
		cout << "CVodeGetQuadB error at 0." << std::endl;
	}

	// Store the T=0 case
	Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> SoutV(Sout, ntps, sMem.params.size());
	Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> St0(sMem.params.size());
	std::copy_n(NV_DATA_S(sMem.qB), NV_LENGTH_S(sMem.qB), St0.data());


	// ###### NOW NEED TO ADD IN EFFECT OF INITIAL CONDITION
	using CppAD::AD;

	size_t Np = sMem.params.size();

	vector<AD<double>> X(Np);
	vector<double> xx(Np);

	std::copy(sMem.params.begin(), sMem.params.end(), X.begin());
	std::copy(sMem.params.begin(), sMem.params.end(), xx.begin());

	CppAD::Independent(X);

	ratesS<AD<double>> rattes = ratesS<AD<double>>(X);

	// Get the data in the right form
	std::array<AD<double>, Nspecies> outAD = solveAutocrine(&rattes);

	vector<AD<double>> outt(1);
	outt[0] = 0.0;

	for (size_t ii = 0; ii < outAD.size(); ii++) {
		outt[0] += outAD[ii] * sMem.activities[ii];
	}

	CppAD::ADFun<double> f(X, outt);

	vector<double> gradZ(Np);
	gradZ = f.Jacobian(xx);
	Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic>> gradZV(gradZ.data(), Np);
	
	// Add in the effect of the initial condition, then finalize output
	for (size_t ii = 0; ii < ntps; ii++) {
		SoutV.row(ii) += gradZV - St0;
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

extern "C" int runCkineParallel (const double * const rxnRatesIn, double tp, size_t nDoses, double *out) {
	int retVal = 1000;
	std::list<std::future<int>> results;

	// Actually run the simulations
	for (size_t ii = 0; ii < nDoses; ii++)
		results.push_back(pool.enqueue(runCkine, &tp, 1, out + Nspecies*ii, rxnRatesIn + ii*Nparams, false));

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
