// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_CELL_CONVERSION_
#define _TEST_CELL_CONVERSION_

#include "ll_cell.h"
#include <cppunit/extensions/HelperMacros.h>

class test_cell_conversion: public CppUnit::TestFixture {
public:
	// types
	typedef ll_cell cell;
	
public:
	void test_getRB();
	void test_getAp();
	void test_getcAp();
	void test_getSubCell();
	void test_getBonds();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_CELL_CONVERSION_
