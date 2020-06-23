// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_FN_WANNIER_MATCHING_
#define _TEST_FN_WANNIER_MATCHING_

#include <cppunit/extensions/HelperMacros.h>

class test_fn_wannier_matching: public CppUnit::TestFixture {
public:
	void test_matchCenters();
	void test_clusterize();
	void test_genCenters();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_FN_WANNIER_MATCHING_
