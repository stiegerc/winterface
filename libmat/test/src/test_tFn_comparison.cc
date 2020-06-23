// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tFn_comparison.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


template<class TT, class FT, class CT>
void test_tFn_comparison<TT,FT,CT>::test_all() {
	const size_t M = genRndST();
	const size_t N = genRndST();
	auto tMat1 = rnd<tMat>(M,N);

	{
		tMat1 = zeros<tMat>(M,N);
		CPPUNIT_ASSERT(!all(tMat1));
		
		tMat1 = ones<tMat>(M,N);
		CPPUNIT_ASSERT(all(tMat1));
		
		tMat1 = eye<tMat>(M,N);
		CPPUNIT_ASSERT(!all(tMat1));
	}
	{
		tMat1 = zeros<tMat>(M,N);
		for (auto r=tMat1.crBegin(), e=tMat1.crEnd(); r!=e; ++r)
			CPPUNIT_ASSERT(!all(*r));
		
		tMat1 = ones<tMat>(M,N);
		for (auto r=tMat1.crBegin(), e=tMat1.crEnd(); r!=e; ++r)
			CPPUNIT_ASSERT(all(*r));
		
		tMat1 = eye<tMat>(M,N);
		for (auto r=tMat1.crBegin(), e=tMat1.crEnd(); r!=e; ++r)
			CPPUNIT_ASSERT(!all(*r));
	}
	{
		tMat1 = zeros<tMat>(M,N);
		for (auto c=tMat1.ccBegin(), e=tMat1.ccEnd(); c!=e; ++c)
			CPPUNIT_ASSERT(!all(*c));
		
		tMat1 = ones<tMat>(M,N);
		for (auto c=tMat1.ccBegin(), e=tMat1.ccEnd(); c!=e; ++c)
			CPPUNIT_ASSERT(all(*c));
		
		tMat1 = eye<tMat>(M,N);
		for (auto c=tMat1.crBegin(), e=tMat1.crEnd(); c!=e; ++c)
			CPPUNIT_ASSERT(!all(*c));
	}
}

template<class TT, class FT, class CT>
void test_tFn_comparison<TT,FT,CT>::test_mall() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	auto tMat1 = rnd<tMat>(M,N);

	{
		tMat1 = zeros<tMat>(0,N);
		const auto act = mall(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(0),act.M());
		CPPUNIT_ASSERT_EQUAL(N,act.N());
		
	}
	{
		tMat1 = zeros<tMat>(M,N);
		const auto act = mall(tMat1);
		
		CPPUNIT_ASSERT_EQUAL(size_t(1),act.M());
		CPPUNIT_ASSERT_EQUAL(N,act.N());
		CPPUNIT_ASSERT(act.logical());
		CPPUNIT_ASSERT(!any(act));
	}
	{
		tMat1 = ones<tMat>(M,N);
		const auto act = mall(tMat1);
		
		CPPUNIT_ASSERT_EQUAL(size_t(1),act.M());
		CPPUNIT_ASSERT_EQUAL(N,act.N());
		CPPUNIT_ASSERT(act.logical());
		CPPUNIT_ASSERT(all(act));
	}
	{
		tMat1 = eye<tMat>(M,N);
		const auto act = mall(tMat1);
		
		CPPUNIT_ASSERT_EQUAL(size_t(1),act.M());
		CPPUNIT_ASSERT_EQUAL(N,act.N());
		CPPUNIT_ASSERT(act.logical());
		CPPUNIT_ASSERT(!any(act));
	}
}

template<class TT, class FT, class CT>
void test_tFn_comparison<TT,FT,CT>::test_nall() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	auto tMat1 = rnd<tMat>(M,N);

	{
		tMat1 = zeros<tMat>(M,0);
		const auto act = nall(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,act.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),act.N());
		
	}
	{
		tMat1 = zeros<tMat>(M,N);
		const auto act = nall(tMat1);
		
		CPPUNIT_ASSERT_EQUAL(M,act.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),act.N());
		CPPUNIT_ASSERT(act.logical());
		CPPUNIT_ASSERT(!any(act));
	}
	{
		tMat1 = ones<tMat>(M,N);
		const auto act = nall(tMat1);
		
		CPPUNIT_ASSERT_EQUAL(M,act.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),act.N());
		CPPUNIT_ASSERT(act.logical());
		CPPUNIT_ASSERT(all(act));
	}
	{
		tMat1 = eye<tMat>(M,N);
		const auto act = nall(tMat1);
		
		CPPUNIT_ASSERT_EQUAL(M,act.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),act.N());
		CPPUNIT_ASSERT(act.logical());
		CPPUNIT_ASSERT(!any(act));
	}
}

template<class TT, class FT, class CT>
void test_tFn_comparison<TT,FT,CT>::test_any() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		const auto tMat1 = zeros<tMat>(M,N);
		CPPUNIT_ASSERT(!any(tMat1));
		
		const auto tMat2 = ones<tMat>(M,N);
		CPPUNIT_ASSERT(any(tMat2));
		
		const auto tMat3 = eye<tMat>(M,N);
		CPPUNIT_ASSERT(any(tMat3));
	}
	{
		const auto tMat1 = zeros<tMat>(M,N);
		for (auto r=tMat1.crBegin(), e=tMat1.crEnd(); r!=e; ++r)
			CPPUNIT_ASSERT(!any(*r));
		
		const auto tMat2 = ones<tMat>(M,N);
		for (auto r=tMat2.crBegin(), e=tMat2.crEnd(); r!=e; ++r)
			CPPUNIT_ASSERT(any(*r));
		
		const auto tMat3 = eye<tMat>(M,M);
		for (auto r=tMat3.crBegin(), e=tMat3.crEnd(); r!=e; ++r)
			CPPUNIT_ASSERT(any(*r));
	}
	{
		const auto tMat1 = zeros<tMat>(M,N);
		for (auto c=tMat1.ccBegin(), e=tMat1.ccEnd(); c!=e; ++c)
			CPPUNIT_ASSERT(!any(*c));
		
		const auto tMat2 = ones<tMat>(M,N);
		for (auto c=tMat2.ccBegin(), e=tMat2.ccEnd(); c!=e; ++c)
			CPPUNIT_ASSERT(any(*c));
		
		const auto tMat3 = eye<tMat>(M,M);
		for (auto c=tMat3.crBegin(), e=tMat3.crEnd(); c!=e; ++c)
			CPPUNIT_ASSERT(any(*c));
	}
}

template<class TT, class FT, class CT>
void test_tFn_comparison<TT,FT,CT>::test_many() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	auto tMat1 = rnd<tMat>(M,N);

	{
		tMat1 = zeros<tMat>(0,N);
		const auto act = many(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(0),act.M());
		CPPUNIT_ASSERT_EQUAL(N,act.N());
		
	}
	{
		tMat1 = zeros<tMat>(M,N);
		const auto act = many(tMat1);
		
		CPPUNIT_ASSERT_EQUAL(size_t(1),act.M());
		CPPUNIT_ASSERT_EQUAL(N,act.N());
		CPPUNIT_ASSERT(act.logical());
		CPPUNIT_ASSERT(!any(act));
	}
	{
		tMat1 = ones<tMat>(M,N);
		const auto act = many(tMat1);
		
		CPPUNIT_ASSERT_EQUAL(size_t(1),act.M());
		CPPUNIT_ASSERT_EQUAL(N,act.N());
		CPPUNIT_ASSERT(act.logical());
		CPPUNIT_ASSERT(any(act));
	}
	{
		tMat1 = eye<tMat>(M,M);
		const auto act = many(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(1),act.M());
		CPPUNIT_ASSERT_EQUAL(M,act.N());
		CPPUNIT_ASSERT(act.logical());
		CPPUNIT_ASSERT(any(act));
	}
}

template<class TT, class FT, class CT>
void test_tFn_comparison<TT,FT,CT>::test_nany() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	auto tMat1 = rnd<tMat>(M,N);

	{
		tMat1 = zeros<tMat>(M,0);
		const auto act = nany(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,act.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),act.N());
	}
	{
		tMat1 = zeros<tMat>(M,N);
		const auto act = nany(tMat1);
		
		CPPUNIT_ASSERT_EQUAL(M,act.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),act.N());
		CPPUNIT_ASSERT(act.logical());
		CPPUNIT_ASSERT(!any(act));
	}
	{
		tMat1 = ones<tMat>(M,N);
		const auto act = nany(tMat1);
		
		CPPUNIT_ASSERT_EQUAL(M,act.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),act.N());
		CPPUNIT_ASSERT(act.logical());
		CPPUNIT_ASSERT(any(act));
	}
	{
		tMat1 = eye<tMat>(M,M);
		const auto act = nany(tMat1);
		
		CPPUNIT_ASSERT_EQUAL(M,act.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),act.N());
		CPPUNIT_ASSERT(act.logical());
		CPPUNIT_ASSERT(any(act));
	}
}

template<class TT, class FT, class CT>
void test_tFn_comparison<TT,FT,CT>::test_find() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	auto tMat1 = (rnd<tMat>(M,N).gt(TT(.5)));

	std::vector<size_t> ref;
	for (size_t i=0; i!=tMat1.L(); ++i)
		if (ops::nz(tMat1[i])) ref.push_back(i);
	const auto act = find(tMat1);

	CPPUNIT_ASSERT(ref==act);
}

template<class TT, class FT, class CT>
void test_tFn_comparison<TT,FT,CT>::test_nnz() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	const auto tMat1 = rand<tMat>(M,N,-5.0*mtol(),5.0*mtol());

	CPPUNIT_ASSERT_EQUAL(size_t(sum(tMat1.neq(0.0))),nnz(tMat1));
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tFn_comparison<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tFn_comparison>(
		"test_all", &test_tFn_comparison<TT,FT,CT>::test_all));
	suite->addTest(new CppUnit::TestCaller<test_tFn_comparison>(
		"test_mall", &test_tFn_comparison<TT,FT,CT>::test_mall));
	suite->addTest(new CppUnit::TestCaller<test_tFn_comparison>(
		"test_nall", &test_tFn_comparison<TT,FT,CT>::test_nall));
	suite->addTest(new CppUnit::TestCaller<test_tFn_comparison>(
		"test_any", &test_tFn_comparison<TT,FT,CT>::test_any));
	suite->addTest(new CppUnit::TestCaller<test_tFn_comparison>(
		"test_many", &test_tFn_comparison<TT,FT,CT>::test_many));
	suite->addTest(new CppUnit::TestCaller<test_tFn_comparison>(
		"test_nany", &test_tFn_comparison<TT,FT,CT>::test_nany));
	suite->addTest(new CppUnit::TestCaller<test_tFn_comparison>(
		"test_find", &test_tFn_comparison<TT,FT,CT>::test_find));
	suite->addTest(new CppUnit::TestCaller<test_tFn_comparison>(
		"test_nnz", &test_tFn_comparison<TT,FT,CT>::test_nnz));
	
	return suite;
}
