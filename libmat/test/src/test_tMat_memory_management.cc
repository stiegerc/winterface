// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_memory_management.h"
#include "testTools.h"
#include <iostream>

using namespace lm__::test;

template<class TT, class FT, class CT>
void test_tMat_memory_management<TT,FT,CT>::test_reserve() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M,N);
	
	tMat1.reserve(N);
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
	CPPUNIT_ASSERT_EQUAL(M*N,tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.ccap());

	tMat1.reserve(N+1);
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
	CPPUNIT_ASSERT_EQUAL(M*(N+1),tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(N+1,tMat1.ccap());
}

template<class TT, class FT, class CT>
void test_tMat_memory_management<TT,FT,CT>::test_shrink_to_fit() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M,N);
	
	tMat1.reserve(N+3);
	tMat1.shrink_to_fit();
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
	CPPUNIT_ASSERT_EQUAL(M*N,tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.ccap());
}

template<class TT, class FT, class CT>
void test_tMat_memory_management<TT,FT,CT>::test_operator_shift() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M,N);
	
	tMat1 << N;
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
	CPPUNIT_ASSERT_EQUAL(M*N,tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.ccap());

	tMat1 << (N+1);
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
	CPPUNIT_ASSERT_EQUAL(M*(N+1),tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(N+1,tMat1.ccap());
}

template<class TT, class FT, class CT>
void test_tMat_memory_management<TT,FT,CT>::test_move() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M,N);
	const auto ptr = tMat1.data();

	const auto data_ = tMat1.move();
	CPPUNIT_ASSERT_EQUAL(size_t(0),tMat1.N());
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(size_t(0),tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(size_t(0),tMat1.ccap());
	CPPUNIT_ASSERT_EQUAL(ptr,data_);
	CPPUNIT_ASSERT_EQUAL(long(nullptr),long(tMat1.data()));

	delete[] data_;
}

template<class TT, class FT, class CT>
CppUnit::Test* test_tMat_memory_management<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tMat_memory_management>(
		"test_reserve", &test_tMat_memory_management<TT,FT,CT>::test_reserve));
	suite->addTest(new CppUnit::TestCaller<test_tMat_memory_management>(
		"test_shrink_to_fit", &test_tMat_memory_management<TT,FT,CT>::test_shrink_to_fit));
	suite->addTest(new CppUnit::TestCaller<test_tMat_memory_management>(
		"test_operator_shift", &test_tMat_memory_management<TT,FT,CT>::test_operator_shift));
	suite->addTest(new CppUnit::TestCaller<test_tMat_memory_management>(
		"test_move", &test_tMat_memory_management<TT,FT,CT>::test_move));
	
	return suite;
}
