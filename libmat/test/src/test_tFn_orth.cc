// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tFn_orth.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


template<class TT, class FT, class CT>
void test_tFn_orth<TT,FT,CT>::test_gsorth() {
	
	const size_t M = genRndST();

	auto tMat1 = rnd_b<tMat>(M,M);
	gsorth(tMat1);
	
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(M,tMat1.N());

	CPPUNIT_ASSERT_DELTA(std::abs(det(tMat1)),RE__(1.0),delta);
}

template<class TT, class FT, class CT>
void test_tFn_orth<TT,FT,CT>::test_complement() {
	auto S = genRndMEST();
	if (S.N>S.M) std::swap(S.M,S.N);

	// check complement to square is empty
	{
		const auto tMat1 = rnd_sqb<tMat>(S.M);
		const auto tMat2 = complement(tMat1);
		
		CPPUNIT_ASSERT(tMat2.empty());
		CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),tMat2.N());
	}

	const auto tMat1 = rnd_b<tMat>(S.M,S.N);
	const auto tMat2 = complement(tMat1);

	CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.M());
	CPPUNIT_ASSERT_EQUAL(tMat1.M()-tMat1.N(),tMat2.N());
	CPPUNIT_ASSERT(zeros<tMat>(tMat1.N(),tMat2.N())==T(tMat1).prod(tMat2));
	CPPUNIT_ASSERT(mnorm(tMat2)==TT(1.0));

	const auto tMat3 = ncat(tMat1,tMat2);
	CPPUNIT_ASSERT(tMat3.square());
	CPPUNIT_ASSERT_EQUAL(S.M,rank(tMat3));
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tFn_orth<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tFn_orth>(
		"test_gsorth", &test_tFn_orth<TT,FT,CT>::test_gsorth));
	suite->addTest(new CppUnit::TestCaller<test_tFn_orth>(
		"test_complement", &test_tFn_orth<TT,FT,CT>::test_complement));
	
	return suite;
}
