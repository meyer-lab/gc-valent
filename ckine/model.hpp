#include <array>
#include <string>


constexpr size_t Nlig = 6; // Number of ligands

// Measured in the literature
constexpr double kfbnd = 0.60; // Assuming on rate of 10^7 M-1 sec-1

// Literature values for k values for IL-15
constexpr double k13rev = kfbnd * 0.065; // based on the multiple papers suggesting 30-100 pM
constexpr double k14rev = kfbnd * 438; // doi:10.1038/ni.2449, 438 nM

// Literature values for IL-7
constexpr double k25rev = kfbnd * 59; // DOI:10.1111/j.1600-065X.2012.01160.x, 59 nM

// Literature value for IL-9
constexpr double k29rev = kfbnd * 0.1; // DOI:10.1073/pnas.89.12.5690, ~100 pM

// Literature value for IL-4
constexpr double k32rev = kfbnd * 1.0; // DOI: 10.1126/scisignal.aal1253 (human)

// Literature value for IL-21
constexpr double k34rev = kfbnd * 0.07; // DOI: 10.1126/scisignal.aal1253 (human)


struct bindingRates {
	double kfwd;
	double k1rev;
	double k2rev;
	double k4rev;
	double k5rev;
	double k10rev;
	double k11rev;
	double k16rev;
	double k17rev;
	double k20rev;
	double k21rev;
	double k22rev;
	double k23rev;
	double k24rev;
	double k27rev;
	double k31rev;
	double k33rev;
	double k35rev;
};


class ratesS {
public:
	std::array<double, Nlig> ILs; // IL2, 15, 7, 9, 4, 21
	bindingRates surface, endosome;
	double endo;
	double activeEndo;
	double sortF;
	double kRec;
	double kDeg;
	std::array<double, 8> Rexpr;

	explicit ratesS(const double * const rxntfR) {
		std::copy_n(rxntfR, ILs.size(), ILs.begin());
		surface.kfwd = rxntfR[6];
		surface.k1rev = kfbnd * 10; // doi:10.1016/j.jmb.2004.04.038, 10 nM
		surface.k2rev = kfbnd * 144; // doi:10.1016/j.jmb.2004.04.038, 144 nM
		surface.k4rev = rxntfR[7];
		surface.k5rev = rxntfR[8];
		surface.k16rev = rxntfR[9];
		surface.k17rev = rxntfR[10];
		surface.k22rev = rxntfR[11];
		surface.k23rev = rxntfR[12];
		surface.k27rev = rxntfR[13];
		surface.k31rev = rxntfR[14];
		surface.k33rev = rxntfR[15];
		surface.k35rev = rxntfR[16];

		// These are probably measured in the literature
		surface.k10rev = 12.0 * surface.k5rev / 1.5; // doi:10.1016/j.jmb.2004.04.038
		surface.k11rev = 63.0 * surface.k5rev / 1.5; // doi:10.1016/j.jmb.2004.04.038

		// IL15
		// To satisfy detailed balance these relationships should hold
		// _Based on initial assembly steps
		surface.k24rev = k13rev * surface.k23rev / k14rev; // loop for IL15_IL15Ra_IL2Rb still holds

		// _Based on formation of full complex
		surface.k21rev = surface.k22rev * surface.k23rev / surface.k16rev;
		surface.k20rev = surface.k22rev * surface.k24rev / surface.k17rev;

		// Set the rates
		endo = rxntfR[17];
		activeEndo = rxntfR[18];
		sortF = rxntfR[19];
		kRec = rxntfR[20];
		kDeg = rxntfR[21];

		if (sortF > 1.0) {
			throw std::runtime_error(std::string("sortF is a fraction and cannot be greater than 1.0."));
		}

		// Expression: IL2Ra, IL2Rb, gc, IL15Ra, IL7Ra, IL9R, IL4Ra, IL21Ra
		std::copy_n(rxntfR + 22, 8, Rexpr.begin());

		endosome = surface;
	}

	void print() {
		std::cout << "kfwd: " << surface.kfwd << std::endl;
		std::cout << "k4rev: " << surface.k4rev << std::endl;
		std::cout << "k5rev: " << surface.k5rev << std::endl;
		std::cout << "k10rev: " << surface.k10rev << std::endl;
		std::cout << "k11rev: " << surface.k11rev << std::endl;
		std::cout << "k16rev: " << surface.k16rev << std::endl;
		std::cout << "k17rev: " << surface.k17rev << std::endl;
		std::cout << "k20rev: " << surface.k20rev << std::endl;
		std::cout << "k21rev: " << surface.k21rev << std::endl;
		std::cout << "k22rev: " << surface.k22rev << std::endl;
		std::cout << "k23rev: " << surface.k23rev << std::endl;
		std::cout << "k24rev: " << surface.k24rev << std::endl;
		std::cout << "k27rev: " << surface.k27rev << std::endl;
		std::cout << "k31rev: " << surface.k31rev << std::endl;
		std::cout << "k33rev: " << surface.k33rev << std::endl;
		std::cout << "k35rev: " << surface.k35rev << std::endl;
		std::cout << "endo: " << endo << std::endl;
		std::cout << "activeEndo: " << activeEndo << std::endl;
		std::cout << "sortF: " << sortF << std::endl;
		std::cout << "kRec: " << kRec << std::endl;
		std::cout << "kDeg: " << kDeg << std::endl;
	}
};


constexpr double tolIn = 1.5E-7;
constexpr double internalV = 623.0; // Same as that used in TAM model
constexpr double internalFrac = 0.5; // Same as that used in TAM model

constexpr size_t Nparams = 30; // length of rxntfR vector
constexpr size_t Nspecies = 62; // number of complexes in surface + endosome + free ligand
constexpr size_t halfL = 28; // number of complexes on surface alone

extern "C" int runCkine (double *tps, size_t ntps, double *out, const double * const rxnRatesIn, const bool sensi, double *sensiOut);
extern "C" int runCkinePretreat (const double pret, const double tt, double * const out, const double * const rxnRatesIn, const double * const postStim);
