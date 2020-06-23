// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_FN_VECTOR_SUMS_
#define _TEST_FN_VECTOR_SUMS_

#include "lm_defs.h"
#include <cppunit/extensions/HelperMacros.h>

class test_fn_vector_sums: public CppUnit::TestFixture {
public:
	// tests
	void test_sum_r();
	void test_sum_c();

	// delta
#ifdef DOUBLE__
	static constexpr RE__ delta = 1e-10;
#endif
#ifdef SINGLE__
	static constexpr RE__ delta = 1e-5f;
#endif

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_FN_VECTOR_SUMS_
