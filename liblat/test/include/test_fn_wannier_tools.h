// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_FN_WANNIER_TOOLS_
#define _TEST_FN_WANNIER_TOOLS_

#include <cppunit/extensions/HelperMacros.h>

class test_fn_wannier_tools: public CppUnit::TestFixture {
public:
	void test_wiToT_TToWi();
	void test_checkWi();
	void test_checkHr();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_FN_WANNIER_TOOLS_
