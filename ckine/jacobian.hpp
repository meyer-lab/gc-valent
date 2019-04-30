#include <adept.h>
#include "reaction.hpp"

typedef Eigen::Matrix<double, Nspecies, Nspecies, Eigen::RowMajor> JacMat;


template <class T>
void fullJacobian(const double * const yv, const ratesS<double> * const r, T &out) {
	adept::Stack stack;

	std::array<adept::adouble, Nspecies> y, dydt;

	adept::set_values(&y[0], Nspecies, yv);

	stack.new_recording();

	// Get the data in the right form
	fullModel(y.data(), &r, dydt.data());

	stack.independent(&y[0], Nspecies);
	stack.dependent(&dydt[0], Nspecies);

	Eigen::Matrix<double, Nspecies, Nspecies> jac_Auto, jac_Ana;
	stack.jacobian(out.data());
}
