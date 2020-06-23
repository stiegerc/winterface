// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_MTOL_
#define _TEST_MTOL_

#include "lm_defs.h"
#include <cppunit/extensions/HelperMacros.h>

class test_mtol: public CppUnit::TestFixture {
public:
	// tests
	void test_all();

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_MTOL_
