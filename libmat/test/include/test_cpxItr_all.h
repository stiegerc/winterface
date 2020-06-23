// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_CPXITR_ALL_
#define _TEST_CPXITR_ALL_

#include "lm_cpxItr.h"
#include <cstddef>
#include <cppunit/extensions/HelperMacros.h>

template<ptrdiff_t s>
class test_cpxItr_all: public CppUnit::TestFixture {
public:
	// types
	typedef lm_c_cpxItr<s> c_cpxItr;
	typedef lm_cpxItr<s> cpxItr;

	// tests
	void test_c_ctor_assign();
	void test_ctor_assign();
	void test_swap_distance();
	void test_all_dereference();

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_CPXITR_ALL_
