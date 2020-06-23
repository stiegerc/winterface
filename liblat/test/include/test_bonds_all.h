// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_BONDS_ALL_
#define _TEST_BONDS_ALL_

#include <cppunit/extensions/HelperMacros.h>

class test_bonds_all: public CppUnit::TestFixture {
public:
	void test_i_i_R();
	void test_ctor_info_print();
	void test_mod_conv();
	void test_getBondCenters();
	void test_search();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_BONDS_ALL_
