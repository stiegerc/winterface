// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_fn_wannier_tools.h"
#include "ll_fn.h"
#include "ll_io.h"
#include "ll_cell.h"
#include "ll_testTools.h"
#include "lm_testTools.h"


using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;


void test_fn_wannier_tools::test_wiToT_TToWi() {
	const size_t N = genRndST(50,100);
	const size_t n = genRndST(1,N);

	const auto I = genRndWi(N,n);
	
	const auto T = wiToT(I);
	CPPUNIT_ASSERT_EQUAL(N,T.size());

	for (size_t i=0; i!=T.size(); ++i)
		CPPUNIT_ASSERT(std::find(I[T[i]].begin(),I[T[i]].end(),i)!=I[T[i]].end());

	const auto J = TToWi(T);
	CPPUNIT_ASSERT(J==I);
}

void test_fn_wannier_tools::test_checkWi() {
	const size_t N = genRndST(50,100);
	const size_t n = genRndST(1,N);
	const auto I = genRndWi(N,n);
	
	CPPUNIT_ASSERT(checkWi(I));
	CPPUNIT_ASSERT(checkWi(I,N));
}

void test_fn_wannier_tools::test_checkHr() {
	const std::string ipath = "data/w90/mos2/";

	const auto Hr = readHr(ipath+"wannier90_hr.dat");
	const auto M = checkHr(Hr);

	for (size_t i=0; i!=Hr.size(); ++i) {
		CPPUNIT_ASSERT_EQUAL(max(abs(real(Hr[i]))),M(0,i));
		CPPUNIT_ASSERT_EQUAL(max(abs(imag(Hr[i]))),M(1,i));
	}
}


const char* test_fn_wannier_tools::test_id() noexcept {
	return "test_fn_wannier_tools";
}

CppUnit::Test* test_fn_wannier_tools::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_fn_wannier_tools>(
		"test_wiToT_TToWi", &test_fn_wannier_tools::test_wiToT_TToWi));
	suite->addTest(new CppUnit::TestCaller<test_fn_wannier_tools>(
		"test_checkWi", &test_fn_wannier_tools::test_checkWi));
	suite->addTest(new CppUnit::TestCaller<test_fn_wannier_tools>(
		"test_checkHr", &test_fn_wannier_tools::test_checkHr));
	
	return suite;
}
