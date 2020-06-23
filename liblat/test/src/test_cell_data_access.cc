// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_cell_data_access.h"
#include "ll_testTools.h"
#include "aux_io.h"
#include <iostream>

using namespace lm__;
using namespace ll__;
using namespace ll__::test;
using namespace lm__::test;
using namespace aux;


void test_cell_data_access::test_B_Ap() {
	const auto M = genRndST(1,5);
	const auto N = genRndST(3,20);
	
	auto B = rnd_sqb<fMat>(M,-.5,.5);
	fMat Ap;
	do Ap = rand<fMat>(M,N,0.0,.99);
	while (!cunique(Ap));

	std::sort(Ap.cBegin(),Ap.cEnd(),vcmp);
	std::vector<size_t> T = {N};

	const ll_cell tCell1(B,Ap,T);

	CPPUNIT_ASSERT(B==tCell1.B());
	CPPUNIT_ASSERT(Ap==tCell1.Ap());
}

void test_cell_data_access::test_cAt_cFront_cBack() {
	const auto M = genRndST(1,5);
	const auto N = genRndST(3,20);
	
	auto B = rnd_sqb<fMat>(M,-.5,.5);
	fMat Ap;
	do Ap = rand<fMat>(M,N,0.0,.99);
	while (!cunique(Ap));

	std::sort(Ap.cBegin(),Ap.cEnd(),vcmp);
	std::vector<size_t> T = {N};

	const ll_cell tCell1(B,Ap,T);

	const size_t ri = genRndST(0,N-1);
	CPPUNIT_ASSERT(tCell1.Ap().cAt(ri)==tCell1.cAt(ri));
	CPPUNIT_ASSERT(tCell1.Ap().cFront()==tCell1.cFront());
	CPPUNIT_ASSERT(tCell1.Ap().cBack()==tCell1.cBack());
}


const char* test_cell_data_access::test_id() noexcept {
	return "test_cell_data_access";
}

CppUnit::Test* test_cell_data_access::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_cell_data_access>(
		"test_B_Ap", &test_cell_data_access::test_B_Ap));
	suite->addTest(new CppUnit::TestCaller<test_cell_data_access>(
		"test_cAt_cFront_cBack", &test_cell_data_access::test_cAt_cFront_cBack));
	
	return suite;
}
