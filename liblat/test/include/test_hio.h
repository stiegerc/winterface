// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_HIO_
#define _TEST_HIO_

#include <cppunit/extensions/HelperMacros.h>
#include "ll_types.h"

class test_hio: public CppUnit::TestFixture {
public:
	void test_writer();
	void test_getConnectedGrid();
	void test_hctor();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_HIO_
