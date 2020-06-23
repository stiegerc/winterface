// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_CELL_INFORMATION_
#define _TEST_CELL_INFORMATION_

#include "ll_cell.h"
#include <cppunit/extensions/HelperMacros.h>

class test_cell_information: public CppUnit::TestFixture {
public:
	// types
	typedef ll_cell cell;
	
public:
	void test_empty();
	void test_dim();
	void test_N_Nspecies();
	void test_vol_sign();
	void test_find();
	void test_lVec();
	void test_primitive();
	void test_validBasis();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_CELL_INFORMATION_
