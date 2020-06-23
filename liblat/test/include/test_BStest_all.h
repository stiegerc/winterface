// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_BSTEST_ALL_
#define _TEST_BSTEST_ALL_

#include <cppunit/extensions/HelperMacros.h>

class test_BStest_all: public CppUnit::TestFixture {
public:
	void test_exceptions();
	void test_random_expansions();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_BSTEST_ALL_
