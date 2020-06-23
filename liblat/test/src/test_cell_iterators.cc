// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_cell_iterators.h"
#include "ll_testTools.h"
#include "aux_io.h"
#include <iostream>

using namespace ll__;
using namespace ll__::test;
using namespace lm__;
using namespace lm__::test;
using namespace aux;


void test_cell_iterators::test_all() {
	const auto D = genRndST(1,5);
	const auto tCell1 = genRandom(D);

	for (aT t=0; t!=tCell1.Nspecies(); ++t) {
		auto I = tCell1.ind(t);
		CPPUNIT_ASSERT_EQUAL(tCell1.Ap().ccBegin()+I.front(),tCell1.ccBegin(t));
		CPPUNIT_ASSERT_EQUAL(tCell1.Ap().ccBegin()+I.back()+1,tCell1.ccEnd(t));
		CPPUNIT_ASSERT_EQUAL(tCell1.ccBegin(t),tCell1.cBegin(t));
		CPPUNIT_ASSERT_EQUAL(tCell1.ccEnd(t),tCell1.cEnd(t));
	}

	CPPUNIT_ASSERT_EQUAL(tCell1.Ap().cBegin(),tCell1.cBegin());
	CPPUNIT_ASSERT_EQUAL(tCell1.Ap().cEnd(),tCell1.cEnd());
	CPPUNIT_ASSERT_EQUAL(tCell1.Ap().ccBegin(),tCell1.ccBegin());
	CPPUNIT_ASSERT_EQUAL(tCell1.Ap().ccEnd(),tCell1.ccEnd());
}

const char* test_cell_iterators::test_id() noexcept {
	return "test_cell_iterators";
}

CppUnit::Test* test_cell_iterators::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_cell_iterators>(
		"test_all", &test_cell_iterators::test_all));
	
	return suite;
}
