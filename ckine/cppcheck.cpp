#include <iostream>
#include <string>
#include <random>
#include <cmath>
#include <array>
#include <random>
#include <algorithm>
#include <adept.h>
#include <Eigen/Core>
#include <Eigen/Dense>

#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/XmlOutputter.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestRunner.h>

#include "model.hpp"
#include "jacobian.hpp"
#include "reaction.hpp"

using namespace std;

class interfaceTestCase : public CppUnit::TestCase {
public:
	/// Setup method
	void setUp() {
		random_device rd;
		gen = new mt19937(rd());
	}
 
	/// Teardown method
	void tearDown() {
		delete gen;
	}

	// method to create a suite of tests
	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("interfaceTestCase");

		suiteOfTests->addTest(new CppUnit::TestCaller<interfaceTestCase>("testrunCkine",
			&interfaceTestCase::testrunCkine));
		suiteOfTests->addTest(new CppUnit::TestCaller<interfaceTestCase>("testrunCkinePretreat",
			&interfaceTestCase::testrunCkinePretreat));
		suiteOfTests->addTest(new CppUnit::TestCaller<interfaceTestCase>("testJacobian",
			&interfaceTestCase::testJacobian));
		return suiteOfTests;
	}

	mt19937 *gen;

protected:
	void checkRetVal(int retVal, array<double, Nparams> &rxnRatesIn) {
		if (retVal < 0) {
			for (auto i = rxnRatesIn.begin(); i != rxnRatesIn.end(); ++i)
				std::cout << *i << ' ';

			cout << std::endl;
		}
	}

	array<double, Nparams> getParams() {
		array<double, Nparams> rxnRatesIn;
		lognormal_distribution<> dis(0.1, 0.25);

		generate(rxnRatesIn.begin(), rxnRatesIn.end(), [this, &dis]() { return dis(*this->gen); });

		rxnRatesIn[19] = tanh(rxnRatesIn[19])*0.9;

		return rxnRatesIn;
	}

	void testrunCkine() {
		array<double, 7> tps = {{0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0}};
		array<double, Nspecies*tps.size()> output;
		array<double, Nspecies*tps.size()> output2;
		array<double, Nparams> rxnRatesIn;
		array<double, Nparams*Nspecies*tps.size()> soutput;
		array<double, Nparams*Nspecies*tps.size()> soutput2;

		for (size_t ii = 0; ii < 3; ii++) {
			rxnRatesIn = getParams();

			int retVal = runCkine(tps.data(), tps.size(), output.data(), rxnRatesIn.data(), true, soutput.data(), false);

			// Run a second time to make sure we get the same thing
			int retVal2 = runCkine(tps.data(), tps.size(), output2.data(), rxnRatesIn.data(), true, soutput2.data(), false);

			checkRetVal(retVal, rxnRatesIn);
			checkRetVal(retVal2, rxnRatesIn);

			CPPUNIT_ASSERT(retVal >= 0);
			CPPUNIT_ASSERT(retVal2 >= 0);
			CPPUNIT_ASSERT(std::equal(output.begin(), output.end(), output2.begin()));
			CPPUNIT_ASSERT(std::equal(soutput.begin(), soutput.end(), soutput2.begin()));
		}
	}

	void testrunCkinePretreat() {
		lognormal_distribution<> dis(0.6, 0.25);

		array<double, 6> postStim = {{0.1, 0.1, 0.0, 1.0, 1.0, 0.1}};
		array<double, Nspecies> output;
		array<double, Nspecies> output2;
		array<double, Nparams> rxnRatesIn;

		for (size_t ii = 0; ii < 3; ii++) {
			rxnRatesIn = getParams();

			int retVal = runCkinePretreat(10.0, 10.0, output.data(), rxnRatesIn.data(), postStim.data());

			// Run a second time to make sure we get the same thing
			int retVal2 = runCkinePretreat(10.0, 10.0, output2.data(), rxnRatesIn.data(), postStim.data());

			checkRetVal(retVal, rxnRatesIn);
			checkRetVal(retVal2, rxnRatesIn);

			CPPUNIT_ASSERT(retVal >= 0);
			CPPUNIT_ASSERT(retVal2 >= 0);
			CPPUNIT_ASSERT(std::equal(output.begin(), output.end(), output2.begin()));
		}
	}

	void testJacobian() {
		using adept::adouble;

		lognormal_distribution<> dis(0.1, 0.25);
		array<double, Nspecies> yv;
		array<double, Nparams> pIn = getParams();
		vector<double> params(Nparams);
		copy(pIn.begin(), pIn.end(), params.begin());
		ratesS<double> rattes = ratesS<double>(params);
		generate(yv.begin(), yv.end(), [this, &dis]() { return dis(*this->gen); });

		adept::Stack stack;

		array<adouble, Nspecies> y, dydt;

		adept::set_values(&y[0], Nspecies, yv.data());

		stack.new_recording();

		// Get the data in the right form
		fullModel(y.data(), &rattes, dydt.data());

		stack.independent(&y[0], Nspecies);
		stack.dependent(&dydt[0], Nspecies);

		Eigen::Matrix<double, Nspecies, Nspecies> jac_Auto, jac_Ana;
		stack.jacobian(jac_Auto.data());

		fullJacobian(yv.data(), &rattes, jac_Ana);

		CPPUNIT_ASSERT((jac_Auto - jac_Ana).squaredNorm() < 0.000001);
	}
};

// the main method
int main () {
	CppUnit::TextUi::TestRunner runner;

	ofstream outputFile("testResults.xml");
	CppUnit::XmlOutputter* outputter = new CppUnit::XmlOutputter(&runner.result(), outputFile);    
	runner.setOutputter(outputter);

	runner.addTest(interfaceTestCase::suite());
	
	runner.run();

	outputFile.close();

	return 0;
}
