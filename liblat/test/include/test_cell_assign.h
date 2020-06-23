// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_CELL_ASSIGN_
#define _TEST_CELL_ASSIGN_

#include "ll_cell.h"
#include <cppunit/extensions/HelperMacros.h>

class test_cell_assign: public CppUnit::TestFixture {
public:
	// types
	typedef ll_cell unitCell;

public:
	void test_swap();
	void test_operator_equal();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_CELL_ASSIGN_
