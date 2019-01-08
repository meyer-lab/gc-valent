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
using adept::adouble;

constexpr double solveTol = 1.0E-8;

static void errorHandler(int, const char *, const char *, char *, void *);
int Jac(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector, N_Vector, N_Vector);
int fullModelCVode (const double, const N_Vector, N_Vector, void *);
static int fB(double, N_Vector y, N_Vector yB, N_Vector yBdot, void *user_dataB);

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
	N_Vector state, qB, yB;
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
		
		// Set the scalar relative and absolute tolerances
		if (CVodeSStolerances(cvode_mem, solveTol, solveTol) < 0) {
			throw std::runtime_error(string("Error calling CVodeSStolerances in solver_setup."));
		}

		A = SUNDenseMatrix(NV_LENGTH_S(state), NV_LENGTH_S(state));
		LS = SUNDenseLinearSolver(state, A);
		
		// Call CVDense to specify the CVDENSE dense linear solver
		if (CVDlsSetLinearSolver(cvode_mem, LS, A) < 0) {
			throw std::runtime_error(string("Error calling CVDlsSetLinearSolver in solver_setup."));
		}

		CVDlsSetJacFn(cvode_mem, Jac);

		CVodeSetMaxNumSteps(cvode_mem, 80000);

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

		// CVodeAdjInit to update CVODES memory block by allocting the internal memory needed for backward integration
		constexpr int steps = 10; // no. of integration steps between two consecutive ckeckpoints
		if (CVodeAdjInit(cvode_mem, steps, CV_HERMITE) < 0) {
			throw std::runtime_error(string("Error calling CVodeAdjInit in solver_setup."));
		}
	}


	void backward (double TB1) {
		indexB = 1;
		yB = N_VNew_Serial(Nspecies); // Initialize yB
		qB = N_VNew_Serial(params.size()); // Initialize qB
		std::copy_n(activities.begin(), Nspecies, NV_DATA_S(yB));
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
		if (CVodeSStolerancesB(cvode_mem, indexB, solveTol, solveTol) < 0) {
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

		CVodeSetMaxNumStepsB(cvode_mem, indexB, 80000);
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

	double getActivity() {
		return std::inner_product(activities.begin(), activities.end(), NV_DATA_S(state), 0.0);
	}

	~solver() {
		if (sensi) {
			CVodeSensFree(cvode_mem);
			N_VDestroy_Serial(qB);
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


// fB routine. Compute fB(t,y,yB). 
static int fB(double t, N_Vector y, N_Vector yB, N_Vector yBdot, void *user_data) {
	ratesS<double> rattes = static_cast<solver *>(user_data)->getRates();

	eigenVC yBv(NV_DATA_S(yB), Nspecies);
	eigenVC yBdotv(NV_DATA_S(yBdot), Nspecies);

	JacMat jac;

	// Actually get the Jacobian
	fullJacobian(NV_DATA_S(y), &rattes, jac);

	yBdotv = -jac*yBv;

	return 0;
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


int Jac(double, N_Vector y, N_Vector, SUNMatrix J, void *user_data, N_Vector, N_Vector, N_Vector) {
	ratesS<double> rattes = static_cast<solver *>(user_data)->getRates();

	Eigen::Map<JacMat> jac(SM_DATA_D(J));

	// Actually get the Jacobian
	fullJacobian(NV_DATA_S(y), &rattes, jac);

	jac.transposeInPlace();

	return 0;
}


int fullModelCVode(const double, const N_Vector xx, N_Vector dxxdt, void *user_data) {
	ratesS<double> rattes = static_cast<solver *>(user_data)->getRates();

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

typedef Eigen::Matrix<double, Nspecies, Eigen::Dynamic> x0JacM;

x0JacM xNotp (vector<double> &params) {
	size_t Np = params.size();

	adept::Stack stack;

	vector<adouble> X(Np);
	adept::set_values(&X[0], Np, params.data());

	stack.new_recording();

	ratesS<adouble> rattes = ratesS<adouble>(X);

	// Get the data in the right form
	std::array<adouble, Nspecies> outAD = solveAutocrine(&rattes);

	stack.independent(&X[0], Np);
	stack.dependent(&outAD[0], Nspecies);

	x0JacM gradZV(Nspecies, Np);
	stack.jacobian(gradZV.data());

	return gradZV;
}



extern "C" int runCkineS (const double * const tps, const size_t ntps, double * const out, double * const Sout, double * const actV, const double * const rxnRatesIn, bool IL2case) {
	size_t itps = 0;

	std::vector<double> v;
	std::array<double, Nspecies> actVv;
	std::copy_n(actV, Nspecies, actVv.begin());

	v = std::vector<double>(rxnRatesIn, rxnRatesIn + Nparams);

	solver sMem(v, actVv);

	if (tps[0] < std::numeric_limits<double>::epsilon()) {
		out[0] = sMem.getActivity();

		itps = 1;
	}

	for (; itps < ntps; itps++) {
		if (tps[itps] < sMem.tret) {
			cout << "Can't go backwards." << std::endl;
			return -1;
		}
		
		if (sMem.CVodeRun(tps[itps]) < 0) return -1;

		// Copy out result
		out[itps] = sMem.getActivity();
	}

	if (sMem.CVodeRun(tps[ntps-1] + 10.0) < 0) return -1;
	sMem.backward(tps[ntps-1] + 10.0);

	/*for (int bitps = ntps - 1; bitps >= 0; bitps--) {
		if (CVodeB(sMem.cvode_mem, tps[bitps], CV_NORMAL) < 0) {
			cout << "CVodeB error." << std::endl;
		}
	}*/

	if (CVodeB(sMem.cvode_mem, 0.0, CV_NORMAL) < 0) {
		cout << "CVodeB error at 0." << std::endl;
	}

	CVodeGetB(sMem.cvode_mem, sMem.indexB, &sMem.tret, sMem.yB);

	// Store the T=0 case
	Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> SoutV(Sout, ntps, sMem.params.size());
	Eigen::Matrix<double, Nspecies, 1> St0;
	std::copy_n(NV_DATA_S(sMem.yB), Nspecies, St0.data());

	// ###### NOW NEED TO ADD IN EFFECT OF INITIAL CONDITION
	x0JacM x0p = xNotp(sMem.params);
	
	// Add in the effect of the initial condition, then finalize output
	for (size_t ii = 0; ii < ntps; ii++) {
		SoutV.row(ii) = St0.transpose()*x0p; // gradZV 
	}

	//cout << endl << St0 << endl << endl;
	cout << endl << SoutV << endl;

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
