// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_CELL_CTOR_
#define _TEST_CELL_CTOR_

#include "ll_cell.h"
#include <cppunit/extensions/HelperMacros.h>

class test_cell_ctor: public CppUnit::TestFixture {
public:
	void test_default();
	void test_B_Ap_();
	void test_file();
	void test_copy();
	void test_move();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_CELL_CTOR_
