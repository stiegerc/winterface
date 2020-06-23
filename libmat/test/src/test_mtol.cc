// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_mtol.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


// tests
void test_mtol::test_all() {
#ifndef NVARTOL__
	// check at least one element is on the stack
	CPPUNIT_ASSERT_EQUAL(size_t(1),mtol_depth());

	// check mtol is set to default
	CPPUNIT_ASSERT_EQUAL(RE__(MTOL__),mtol());
	CPPUNIT_ASSERT_EQUAL(size_t(1),mtol_depth());

	// set mtol to new value, check stack
	set_mtol(1e-3);
	CPPUNIT_ASSERT_EQUAL(RE__(1e-3),mtol());
	CPPUNIT_ASSERT_EQUAL(size_t(2),mtol_depth());

	// set mtol to new value, check stack
	set_mtol(1e-2);
	CPPUNIT_ASSERT_EQUAL(RE__(1e-2),mtol());
	CPPUNIT_ASSERT_EQUAL(size_t(3),mtol_depth());

	// check reset reduces stack, produces previous value
	reset_mtol();
	CPPUNIT_ASSERT_EQUAL(RE__(1e-3),mtol());
	CPPUNIT_ASSERT_EQUAL(size_t(2),mtol_depth());
	
	// check reset reduces stack, produces previous value
	reset_mtol();
	CPPUNIT_ASSERT_EQUAL(RE__(MTOL__),mtol());
	CPPUNIT_ASSERT_EQUAL(size_t(1),mtol_depth());

#ifdef _OPENMP
	// check spawning and collapsing stacks
	const size_t N = genRndST();
	spawn_mtol_stacks(N);
	CPPUNIT_ASSERT_EQUAL(size_t(N),mtol_nstacks());

	collapse_mtol_stacks();
	CPPUNIT_ASSERT_EQUAL(size_t(1),mtol_nstacks());
#endif

#else
	// check mtol is set to default
	CPPUNIT_ASSERT_EQUAL(RE__(MTOL__),mtol());
#endif
}


// test id
const char* test_mtol::test_id() noexcept {
	return "test_mtol";
}

CppUnit::Test* test_mtol::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_mtol>(
		"test_all", &test_mtol::test_all));
	
	return suite;
}
