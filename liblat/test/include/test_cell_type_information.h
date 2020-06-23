// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_CELL_TYPE_INFORMATION_
#define _TEST_CELL_TYPE_INFORMATION_

#include "ll_cell.h"
#include <cppunit/extensions/HelperMacros.h>

class test_cell_type_information: public CppUnit::TestFixture {
public:
	// types
	typedef ll_cell cell;

public:
	void test_types_inds_();
	void test_type();
	void test_id_();
	void test_Ntype();
	void test_leastFreqType_mostFreqType();
	void test_ind();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_CELL_TYPE_INFORMATION_
