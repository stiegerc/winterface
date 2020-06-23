// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_AUX_FOR_EACH_
#define _TEST_AUX_FOR_EACH_

#include <cppunit/extensions/HelperMacros.h>

class test_aux_for_each: public CppUnit::TestFixture {
public:
	void test_all();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_AUX_FOR_EACH_
