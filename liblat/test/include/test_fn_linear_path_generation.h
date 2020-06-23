// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_FN_LINEAR_PATH_GENERATION_
#define _TEST_FN_LINEAR_PATH_GENERATION_

#include <cppunit/extensions/HelperMacros.h>

class test_fn_linear_path_generation: public CppUnit::TestFixture {
public:
	void test_genPath();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_FN_LINEAR_PATH_GENERATION_
