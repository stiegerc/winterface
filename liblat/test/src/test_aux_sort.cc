// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_aux_sort.h"
#include "aux_sort.h"
#include "libmat.h"
#include "lm_testTools.h"
#include <tuple>
#include <cstddef>
#include <vector>
#include <iostream>

using namespace lm__;
using namespace lm__::test;


void test_aux_sort::test_all() {
	const size_t M = genRndST(30,50);
	const size_t N = genRndST(30,50);

	// sorted_order, reorder for size_t
	{
		std::vector<size_t> S(N);
		std::iota(S.begin(),S.end(),size_t(0));
		std::random_shuffle(S.begin(),S.end());

		const auto I = aux::sorted_order(S.cbegin(),S.cend());
		
		CPPUNIT_ASSERT(std::is_permutation(S.begin(),S.end(),I.begin()));
		for (size_t i=0; i!=N; ++i)
			CPPUNIT_ASSERT_EQUAL(S[i],I[i]);

		aux::reorder(S.begin(),I);
		for (size_t i=0; i!=N; ++i)
			CPPUNIT_ASSERT_EQUAL(i,S[i]);
	}

	// sorted order, reorder for double
	{
		auto tMat1 = rnd<fMat>(N,1,-10.0,10.0);

		const auto I = aux::sorted_order(tMat1.cbegin(),tMat1.cend());
		aux::reorder(tMat1.begin(),I);
		CPPUNIT_ASSERT(std::is_sorted(tMat1.cbegin(),tMat1.cend()));
	}

	// sorted order, reorder for columns of integer matrix
	{
		auto tMat1 = round(rnd<fMat>(M,N,-3.0,3.0));
		
		const auto I = aux::sorted_order(tMat1.ccBegin(),tMat1.ccEnd(),vcmp);
		aux::reorder(tMat1.cBegin(),I);
		CPPUNIT_ASSERT(std::is_sorted(tMat1.ccBegin(),tMat1.ccEnd(),vcmp));
	}
}


const char* test_aux_sort::test_id() noexcept {
	return "test_aux_sort";
}

CppUnit::Test* test_aux_sort::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_aux_sort>(
		"test_all", &test_aux_sort::test_all));
	
	return suite;
}
