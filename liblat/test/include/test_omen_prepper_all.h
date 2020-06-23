// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_OMEN_PREPPER_ALL_
#define _TEST_OMEN_PREPPER_ALL_

#include <cppunit/extensions/HelperMacros.h>

class test_omen_prepper_all: public CppUnit::TestFixture {
public:
	void test_ctor_exceptions();
	void test_ctor_properties();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_OMEN_PREPPER_ALL_
