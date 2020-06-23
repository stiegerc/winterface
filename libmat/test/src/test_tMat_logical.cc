// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_logical.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


template<class TT, class FT, class CT>
void test_tMat_logical<TT,FT,CT>::test_logical() {

	const size_t M = genRndST();
	const size_t N = genRndST();
	const auto tMat1 = rnd<tMat>(M,N);
	
	CPPUNIT_ASSERT(!tMat1.logical());
	CPPUNIT_ASSERT((tMat1.lt(TT(.5))).logical());
}
template<class TT, class FT, class CT>
void test_tMat_logical<TT,FT,CT>::test_permutation() {
	{
		const auto S = genRndMEST();
			
		const auto tMat1 = rnd<tMat>(S.M,S.N);
		CPPUNIT_ASSERT(!tMat1.permutation());
	}
	{
		const size_t M = genRndST();
		const auto tMat1 = rnd<tMat>(M,M);

		CPPUNIT_ASSERT(!tMat1.permutation());
	}
	{
		const size_t M = genRndST();
		auto tMat1 = eye<tMat>(M);
			
		for (size_t i=0; i<M; ++i)
			swap(tMat1.cAt(i),tMat1.cAt(genRndST(0,M-1)));

		CPPUNIT_ASSERT(tMat1.permutation());
	}
}
template<class TT, class FT, class CT>
void test_tMat_logical<TT,FT,CT>::test_operator_not() {

	const size_t M = genRndST();
	const size_t N = genRndST();

	const auto tMat1 = rnd<tMat>(M,N);
	const auto tMat2 = (tMat1.gt(0.0));
	const auto tMat3 = (tMat1.leq(0.0));

	CPPUNIT_ASSERT(tMat2==(~tMat3));
}
template<class TT, class FT, class CT>
void test_tMat_logical<TT,FT,CT>::test_operator_and() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		const auto tMat1 = eye<tMat>(M,N);
		auto tMat2 = tMat1;

		tMat2 = tMat1&FT(1.0);
		CPPUNIT_ASSERT(tMat1==tMat2);
	}
	{
		auto tMat1 = rnd<tMat>(M,N);
		auto tMat2 = fMat(tMat1);

		tMat1 = (tMat1.lt(TT(.5)));
		tMat2 = (tMat2.geq(FT(0.5)));

		const auto tMat3 = tMat1&tMat2;
		CPPUNIT_ASSERT(all(~tMat3));
	}
	
	{
		const auto tMat1 = eye<tMat>(M,N);
		auto tMat2 = tMat1;

		tMat2 = tMat2&CT(1.0);
		CPPUNIT_ASSERT(tMat1==tMat2);
	}
	{
		auto tMat1 = rnd<tMat>(M,N);
		auto tMat2 = cMat(tMat1);

		tMat1 = (tMat1.lt(TT(.5)));
		tMat2 = (tMat2.geq(CT(0.5)));

		const auto tMat3 = tMat1&tMat2;
		CPPUNIT_ASSERT(all(~tMat3));
	}
}
template<class TT, class FT, class CT>
void test_tMat_logical<TT,FT,CT>::test_operator_or() {
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		const auto tMat1 = eye<tMat>(M,N);
		auto tMat2 = tMat1;

		tMat2 = tMat1|FT(0.0);
		CPPUNIT_ASSERT(tMat1==tMat2);
	}
	{
		auto tMat1 = rnd<tMat>(M,N);
		auto tMat2 = fMat(tMat1);

		tMat1 = (tMat1.lt(TT(.5)));
		tMat2 = (tMat2.geq(FT(0.5)));

		const auto tMat3 = tMat1|tMat2;
		CPPUNIT_ASSERT(all(tMat3));
	}
	
	{
		const auto tMat1 = eye<tMat>(M,N);
		auto tMat2 = tMat1;

		tMat2 = tMat2|CT(0.0);
		CPPUNIT_ASSERT(tMat1==tMat2);
	}
	{
		auto tMat1 = rnd<tMat>(M,N);
		auto tMat2 = cMat(tMat1);

		tMat1 = (tMat1.lt(TT(.5)));
		tMat2 = (tMat2.geq(CT(0.5)));

		const auto tMat3 = tMat1|tMat2;
		CPPUNIT_ASSERT(all(tMat3));
	}
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tMat_logical<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tMat_logical>(
		"test_logical", &test_tMat_logical<TT,FT,CT>::test_logical));
	suite->addTest(new CppUnit::TestCaller<test_tMat_logical>(
		"test_permutation", &test_tMat_logical<TT,FT,CT>::test_permutation));
	suite->addTest(new CppUnit::TestCaller<test_tMat_logical>(
		"test_operator_not", &test_tMat_logical<TT,FT,CT>::test_operator_not));
	suite->addTest(new CppUnit::TestCaller<test_tMat_logical>(
		"test_operator_and", &test_tMat_logical<TT,FT,CT>::test_operator_and));
	suite->addTest(new CppUnit::TestCaller<test_tMat_logical>(
		"test_operator_or", &test_tMat_logical<TT,FT,CT>::test_operator_or));
	
	return suite;
}
