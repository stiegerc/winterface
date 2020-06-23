// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_FN_R3_ONLY_
#define _TEST_FN_R3_ONLY_

#include "lm_defs.h"
#include <cppunit/extensions/HelperMacros.h>

class test_fn_R3_only: public CppUnit::TestFixture {
public:
	// tests
	void test_cross_crossn();
	void test_getR();
	

	// delta
#ifdef DOUBLE__
	static constexpr RE__ delta = 1e-6;
#endif
#ifdef SINGLE__
	static constexpr RE__ delta = 1e-3f;
#endif

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_FN_R3_ONLY_
