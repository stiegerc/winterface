// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_PARSER_ALL_
#define _TEST_PARSER_ALL_

#include <cppunit/extensions/HelperMacros.h>

class test_parser_all: public CppUnit::TestFixture {
public:
	void test_tParser();
	void test_screenFile();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_PARSER_ALL_
