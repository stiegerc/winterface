// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_CELL_MEMORY_MANAGEMENT_
#define _TEST_CELL_MEMORY_MANAGEMENT_

#include "ll_cell.h"
#include <cppunit/extensions/HelperMacros.h>

class test_cell_memory_management: public CppUnit::TestFixture {
public:
	// types
	typedef ll_cell cell;
	
public:
	void test_moveB();
	void test_moveAp();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_CELL_MEMORY_MANAGEMENT_
