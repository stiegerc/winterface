// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_cell_memory_management.h"
#include "ll_testTools.h"
#include "aux_io.h"
#include <iostream>

using namespace ll__;
using namespace ll__::test;
using namespace aux;


void test_cell_memory_management::test_moveB() {
	auto tCell1 = genRandom();

	const auto Bptr = tCell1.B().begin();
	const auto Bc = tCell1.B();
	const auto Bm = tCell1.moveB();

	CPPUNIT_ASSERT(Bc==Bm);
	CPPUNIT_ASSERT_EQUAL(Bptr,Bm.begin());
	CPPUNIT_ASSERT(tCell1.B().empty());
}

void test_cell_memory_management::test_moveAp() {
	auto tCell1 = genRandom();

	const auto Apptr = tCell1.Ap().begin();
	const auto Apc = tCell1.Ap();
	const auto Apm = tCell1.moveAp();

	CPPUNIT_ASSERT(Apc==Apm);
	CPPUNIT_ASSERT_EQUAL(Apptr,Apm.begin());
	CPPUNIT_ASSERT(tCell1.Ap().empty());
}

const char* test_cell_memory_management::test_id() noexcept {
	return "test_cell_memory_management";
}

CppUnit::Test* test_cell_memory_management::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_cell_memory_management>(
		"test_moveB", &test_cell_memory_management::test_moveB));
	suite->addTest(new CppUnit::TestCaller<test_cell_memory_management>(
		"test_moveAp", &test_cell_memory_management::test_moveAp));
	
	return suite;
}
