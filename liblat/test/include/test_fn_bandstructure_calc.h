// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_FN_BANDSTRUCTURE_CALC_
#define _TEST_FN_BANDSTRUCTURE_CALC_

#include <cppunit/extensions/HelperMacros.h>

class test_fn_bandstructure_calc: public CppUnit::TestFixture {
public:
	void test_findBandEdges();
	void test_calcBS();
	void test_calcFoldedBS();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_FN_BANDSTRUCTURE_CALC_
