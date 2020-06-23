// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_FN_EIGENVALUE_COMPUTATION_
#define _TEST_FN_EIGENVALUE_COMPUTATION_

#include "lm_defs.h"
#include <cppunit/extensions/HelperMacros.h>

class test_fn_eigenvalue_computation: public CppUnit::TestFixture {
public:
	// tests
	void test_eig();
	void test_eigr();
	void test_eigl();
	void test_eigrl();
	void test_eigh();
	void test_eighv();

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

#endif // _TEST_FN_EIGENVALUE_COMPUTATION_
