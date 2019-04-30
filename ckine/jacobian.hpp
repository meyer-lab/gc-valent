#include "reaction.hpp"

typedef Eigen::Matrix<double, Nspecies, Nspecies, Eigen::RowMajor> JacMat;

template <class T>
void jacobian(const double * const y, const bindingRates<double> * const r, T &out, const double * const ILs) {
	// IL2 in nM
	const double IL2Ra = y[0];
	const double IL2Rb = y[1];
	const double gc = y[2];
	const double IL2_IL2Ra = y[3];
	const double IL2_IL2Rb = y[4];
	const double IL2_IL2Ra_gc = y[6];
	const double IL2_IL2Rb_gc = y[7];
	
	// IL15 in nM
	const double IL15Ra = y[9];
	const double IL15_IL15Ra = y[10];
	const double IL15_IL2Rb = y[11];
	const double IL15_IL15Ra_gc = y[13];
	const double IL15_IL2Rb_gc = y[14];

	auto complexCkine = [&out, &gc, &r, &y](const size_t ij, const size_t RaIDX, const double revOne, const double revTwo) {
		out(RaIDX, ij+1) = -r->kfwd * y[RaIDX]; // Ra with respect to IL_Rb
		out(RaIDX, ij+4) = -r->kfwd * y[RaIDX]; // Ra with respect to IL_Rb_gc

		out(2, 2) -= r->kfwd * (y[ij] + y[ij+1] + y[ij+2]); // gc with respect to gc
		out(2, ij) = - r->kfwd * gc; // gc with respect to IL_Ra
		out(2, ij+1) = - r->kfwd * gc; // gc with respect to IL_Rb
		out(2, ij+2) = - r->kfwd * gc; // gc with respect to IL_Ra_Rb

		out(1, ij) = - r->kfwd * y[1]; // Rb with respect to IL_Ra
		out(1, ij+3) = - r->kfwd * y[1]; // Rb with respect to IL_Ra_gc

		out(ij, 1) = -r->kfwd * y[ij]; // IL_Ra with respect to Rb
		out(ij, 2) = - r->kfwd * y[ij]; // IL_Ra with respect to gc

		out(ij+5, ij+2) = r->kfwd * gc; // IL_Ra_Rb_gc with respect to IL_Ra_Rb

		out(ij+4, 2) = r->kfwd * y[ij+1]; // IL_Rb_gc with respect to gc
		out(ij+4, RaIDX) = -r->kfwd * y[ij+4]; // IL_Rb_gc with respect to Ra
		out(ij+4, ij+1) = r->kfwd * gc; // IL_Rb_gc with respect to IL_Rb
		out(ij+4, ij+4) = -r->kfwd * y[RaIDX] - revOne; // IL_Rb_gc with respect to IL_Rb_gc
		out(ij+4, ij+5) = revTwo; // IL_Rb_gc with respect to IL_Ra_Rb_gc

		out(2, ij+4) = revOne; // gc with respect to IL_Rb_gc
		out(ij+1, ij+4) = revOne; // IL_Rb with respect to IL_Rb_gc
		out(RaIDX, ij+5) = revTwo; // Ra with respect to IL_Ra_Rb_gc

		out(ij+5, RaIDX) = r->kfwd * y[ij+4]; // IL_Ra_Rb_gc with respect to Ra
		out(ij+5, 1) = r->kfwd * y[ij+3]; // IL_Ra_Rb_gc with respect to Rb
		out(ij+5, 2) = r->kfwd * y[ij+2]; // IL_Ra_Rb_gc with respect to gc
		out(ij+5, ij+3) = r->kfwd * y[1]; // IL_Ra_Rb_gc with respect to IL_Ra_gc
		out(ij+5, ij+4) = r->kfwd * y[RaIDX]; // IL_Ra_Rb_gc with respect to IL_Rb_gc
		out(ij+2, RaIDX) = r->kfwd * y[ij+1]; // IL_Ra_Rb with respect to Ra
		out(ij+2, 1) = r->kfwd * y[ij]; // IL_Ra_Rb with respect to Rb
		out(ij+2, 2) = -r->kfwd * y[ij+2]; // IL_Ra_Rb with respect to gc
		out(ij+2, ij) = r->kfwd * y[1]; // IL_Ra_Rb with respect to IL_Ra
		out(ij+2, ij+1) = r->kfwd * y[RaIDX]; // IL_Ra_Rb with respect to IL_Rb

		out(ij+1, RaIDX) = -r->kfwd * y[ij+1]; // IL2_IL2Rb with respect to IL2Ra
		out(ij+1, 2) = - r->kfwd * y[ij+1]; // IL2_IL2Rb with respect to gc
		out(ij+3, 1) = -r->kfwd * y[ij+3]; // IL2_IL2Ra_gc with respect to IL2Rb
		out(ij+3, 2) = r->kfwd * y[ij]; // IL2_IL2Ra_gc with respect to gc
		out(ij+3, ij) = r->kfwd * gc; // IL2_IL2Ra_gc with respect to IL2_IL2Ra
	};

	complexCkine(3, 0, r->k5rev, k8rev(r)); // IL2
	complexCkine(10, 9, r->k17rev, k20rev(r)); // IL15
	
	// IL2Ra
	out(0, 0) = -kfbnd * ILs[0] - r->kfwd * IL2_IL2Rb_gc - r->kfwd * IL2_IL2Rb; // IL2Ra with respect to IL2Ra
	out(0, 3) = r->k1rev; // IL2Ra with respect to IL2_IL2Ra
	out(0, 5) = k12rev(r); // IL2Ra with respect to IL2_IL2Ra_IL2Rb
	
	// IL2Rb
	out(1, 1) = -kfbnd * (ILs[0] + ILs[1]) - r->kfwd * (IL2_IL2Ra_gc + IL2_IL2Ra + IL15_IL15Ra_gc + IL15_IL15Ra); // partial derivative of IL2Rb with respect to IL2Rb
	out(1, 4) = r->k2rev; // IL2Rb with respect to IL2_IL2Rb
	out(1, 5) = r->k11rev; // IL2Rb with respect to IL2_IL2Ra_IL2Rb
	out(1, 8) = k9rev(r); // IL2Rb with respect to IL2_IL2Ra_IL2Rb_gc
	out(1, 11) = r->k14rev; // IL2Rb with respect to IL15_IL2Rb
	out(1, 12) = r->k23rev; // IL2Rb with respect to IL15_IL15Ra_IL2Rb
	out(1, 15) = k21rev(r); // IL2Rb with respect to IL15_IL15Ra_IL2Rb_gc
	
	// gc
	out(2, 6) = r->k4rev; // gc with respect to IL2_IL2Ra_gc
	out(2, 8) = r->k10rev; // gc with respect to IL2_IL2Ra_IL2Rb_gc
	out(2, 13) = r->k16rev; // gc with respect to IL15_IL15Ra_gc
	out(2, 15) = r->k22rev; // gc with respect to IL15_IL15Ra_IL2Rb_gc
	
	// IL2_IL2Ra
	out(3, 0) = kfbnd * ILs[0]; // IL2_IL2Ra with respect to IL2Ra
	out(3, 3) = -r->kfwd * IL2Rb - r->kfwd * gc - r->k1rev; // IL2_IL2Ra with respect to IL2_IL2Ra
	out(3, 5) = r->k11rev; // IL2_IL2Ra with respect to IL2_IL2Ra_IL2Rb
	out(3, 6) = r->k4rev; // IL2_IL2Ra with respect to IL2_IL2Ra_gc
	
	// IL2_IL2Rb
	out(4, 1) = kfbnd * ILs[0]; // IL2_IL2Rb with respect to IL2Rb
	out(4, 4) = -r->kfwd * IL2Ra - r->kfwd * gc - r->k2rev; // IL2_IL2Rb with respect to IL2_IL2Rb
	out(4, 5) = k12rev(r); // IL2_IL2Rb with respect to IL2_IL2Ra_IL2Rb
	
	// IL2_IL2Ra_IL2Rb
	out(5, 5) = -r->kfwd * gc - r->k11rev - k12rev(r); // IL2_IL2Ra_IL2Rb with respect to IL2_IL2Ra_IL2Rb
	out(5, 8) = r->k10rev; // IL2_IL2Ra_IL2Rb with respect to IL2_IL2Ra_IL2Rb_gc
	
	// IL2_IL2Ra_gc
	out(6, 6) = -r->kfwd * IL2Rb - r->k4rev; // IL2_IL2Ra_gc with respect to IL2_IL2Ra_gc
	out(6, 8) = k9rev(r); // IL2_IL2Ra_gc with respect to IL2_IL2Ra_IL2Rb_gc
	
	// IL2_IL2Ra_IL2Rb_gc
	out(8, 8) = - k8rev(r) - k9rev(r) - r->k10rev; // IL2_IL2Ra_IL2Rb_gc with respect to IL2_IL2Ra_IL2Rb_gc
	
	// IL15Ra
	out(9, 9) = -kfbnd * ILs[1] - r->kfwd * IL15_IL2Rb_gc - r->kfwd * IL15_IL2Rb; // IL15Ra with respect to IL15Ra
	out(9, 10) = r->k13rev; // IL15Ra with respect to IL15_IL15Ra
	out(9, 12) = k24rev(r); // IL15Ra with respect to IL15_IL15Ra_IL2Rb
	
	// IL15_IL15Ra
	out(10, 9) = kfbnd * ILs[1]; // IL15_IL15Ra with respect to IL15Ra
	out(10, 10) = -r->kfwd * IL2Rb - r->kfwd * gc - r->k13rev; // IL15_IL15Ra with respect to IL15_IL15Ra
	out(10, 12) = r->k23rev; // IL15_IL15Ra with respect to IL15_IL15Ra_IL2Rb
	out(10, 13) = r->k16rev; // IL15_IL15Ra with respect to IL15_IL15Ra_gc
	
	// IL15_IL2Rb
	out(11, 1) = kfbnd * ILs[1]; // IL15_IL2Rb with respect to IL2Rb
	out(11, 11) = -r->kfwd * IL15Ra - r->kfwd * gc - r->k14rev; // IL15_IL2Rb with respect to IL15_IL2Rb
	out(11, 12) = k24rev(r); // IL15_IL2Rb with respect to IL15_IL15Ra_IL2Rb
	
	// IL15_IL15Ra_IL2Rb
	out(12, 12) = -r->kfwd * gc - r->k23rev - k24rev(r); // IL15_IL15Ra_IL2Rb with respect to IL15_IL15Ra_IL2Rb
	out(12, 15) = r->k22rev; // IL15_IL15Ra_IL2Rb with respect to IL15_IL15Ra_IL2Rb_gc
	
	// IL15_IL15Ra_gc
	out(13, 13) = -r->kfwd * IL2Rb - r->k16rev; // IL15_IL15Ra_gc with respect to IL15_IL15Ra_gc
	out(13, 15) = k21rev(r); // IL15_IL15Ra_gc with respect to IL15_IL15Ra_IL2Rb_gc
	
	// IL15_IL15Ra_IL2Rb_gc
	out(15, 15) = - k20rev(r) - k21rev(r) - r->k22rev; // IL15_IL15Ra_IL2Rb_gc with respect to IL15_IL15Ra_IL2Rb_gc


	auto simpleCkine = [&out, &gc, &r, &y](const size_t ij, const double revOne, const double revTwo, const double IL) {
		out(2, 2) -= r->kfwd * y[ij+1]; // gc with respect to gc
		out(2, ij+1) = -r->kfwd * gc; // gc with respect to Ra_IL
		out(2, ij+2) = revTwo; // gc with respect to Ra_gc_IL

		out(ij, ij  ) = -kfbnd * IL; // Ra with respect to Ra
		out(ij, ij+1) = revOne; // IL_Ra

		out(ij+1, 2) = - r->kfwd * y[ij+1]; // IL_Ra with respect to gc
		out(ij+1, ij) = kfbnd * IL; // IL_Ra with respect to Ra
		out(ij+1, ij+1) = - revOne - r->kfwd * gc; // IL_Ra with respect to IL_Ra
		out(ij+1, ij+2) = revTwo; // IL_Ra with respect to IL_Ra_gc

		out(ij+2,    2) = r->kfwd * y[ij+1]; // IL_Ra_gc with respect to gc
		out(ij+2, ij+1) = r->kfwd * gc; // IL_Ra_gc with respect to IL_Ra
		out(ij+2, ij+2) = -revTwo; // IL_Ra_gc with respect to IL_Ra_gc
	};

	simpleCkine(16, r->k25rev, r->k27rev, ILs[2]); // IL7
	simpleCkine(19, r->k29rev, r->k31rev, ILs[3]); // IL9
	simpleCkine(22, r->k32rev, r->k33rev, ILs[4]); // IL4
	simpleCkine(25, r->k34rev, r->k35rev, ILs[5]); // IL21
}


template <class T>
void fullJacobian(const double * const y, const ratesS<double> * const r, T &out) {
	
	adept::Stack stack;

	array<adouble, Nspecies> y, dydt;

	adept::set_values(&y[0], Nspecies, yv.data());

	stack.new_recording();

	// Get the data in the right form
	fullModel(y.data(), &rattes, dydt.data());

	stack.independent(&y[0], Nspecies);
	stack.dependent(&dydt[0], Nspecies);

	Eigen::Matrix<double, Nspecies, Nspecies> jac_Auto, jac_Ana;
	stack.jacobian(T.data());
}
