// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_CELL_ITERATORS_
#define _TEST_CELL_ITERATORS_

#include "ll_cell.h"
#include <cppunit/extensions/HelperMacros.h>

class test_cell_iterators: public CppUnit::TestFixture {
public:
	// types
	typedef ll_cell cell;
	
public:
	void test_all();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_CELL_ITERATORS_
