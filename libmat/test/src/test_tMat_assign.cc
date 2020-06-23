// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_assign.h"
#include "testTools.h"
#include <iostream>

using namespace lm__::test;

template<class TT, class FT, class CT>
void test_tMat_assign<TT,FT,CT>::test_swap() {
	const size_t M1 = genRndST();
	const size_t N1 = genRndST();
	const size_t M2 = genRndST();
	const size_t N2 = genRndST();

	auto tMat1 = rnd<tMat>(M1,N1);
	auto tMat2 = rnd<tMat>(M2,N2);
	
	const auto ptr1 = tMat1.begin();
	const auto ptr2 = tMat2.begin();
	const size_t lcap1 = tMat1.lcap();
	const size_t lcap2 = tMat2.lcap();
	const size_t ccap1 = tMat1.ccap();
	const size_t ccap2 = tMat2.ccap();

	swap(tMat1,tMat2);

	CPPUNIT_ASSERT_EQUAL(M2,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(N2,tMat1.N());
	CPPUNIT_ASSERT_EQUAL(M1,tMat2.M());
	CPPUNIT_ASSERT_EQUAL(N1,tMat2.N());
	CPPUNIT_ASSERT_EQUAL(lcap2,tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(lcap1,tMat2.lcap());
	CPPUNIT_ASSERT_EQUAL(ccap2,tMat1.ccap());
	CPPUNIT_ASSERT_EQUAL(ccap1,tMat2.ccap());
	CPPUNIT_ASSERT_EQUAL(ptr2,tMat1.begin());
	CPPUNIT_ASSERT_EQUAL(ptr1,tMat2.begin());
}

template<class TT, class FT, class CT>
void test_tMat_assign<TT,FT,CT>::test_operator_equal_fArray() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	// from fMat, same size
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		const auto ptr = tMat1.begin();
		const auto lcap = tMat1.lcap();
		const auto ccap = tMat1.ccap();

		tMat1 = tMat2;

		CPPUNIT_ASSERT_EQUAL(lcap,tMat1.lcap());
		CPPUNIT_ASSERT_EQUAL(ccap,tMat1.ccap());
		CPPUNIT_ASSERT_EQUAL(ptr,tMat1.begin());

		CPPUNIT_ASSERT(tMat2==tMat1);
	}
	
	// from fMat, larger size
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<fMat>(2*M,2*N);

		const auto ptr = tMat1.begin();

		tMat1 = tMat2;

		CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
		CPPUNIT_ASSERT_EQUAL(tMat2.ccap(),tMat1.ccap());
		CPPUNIT_ASSERT(ptr!=tMat1.begin());

		CPPUNIT_ASSERT(tMat2==tMat1);
	}

	// from fRow
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		const size_t m = genRndST(0,M-1);
		auto j = tMat2.crBegin()+m;
		
		tMat1 = *j;

		CPPUNIT_ASSERT(*j==tMat1);
	}
	
	// from fCol
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		const size_t n = genRndST(0,N-1);
		auto j = tMat2.ccBegin()+n;
		
		tMat1 = *j;

		CPPUNIT_ASSERT(*j==tMat1);
	}
}

template<class TT, class FT, class CT>
void test_tMat_assign<TT,FT,CT>::test_operator_equal_cArray() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	// from cMat, same size
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		const auto ptr = tMat1.begin();
		const auto lcap = tMat1.lcap();
		const auto ccap = tMat1.ccap();

		tMat1 = tMat2;

		CPPUNIT_ASSERT_EQUAL(lcap,tMat1.lcap());
		CPPUNIT_ASSERT_EQUAL(ccap,tMat1.ccap());
		CPPUNIT_ASSERT_EQUAL(ptr,tMat1.begin());

		CPPUNIT_ASSERT(tMat2==tMat1);
	}
	
	// from cMat, larger size
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<cMat>(2*M,2*N);

		const auto ptr = tMat1.begin();

		tMat1 = tMat2;

		CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
		CPPUNIT_ASSERT_EQUAL(tMat2.ccap(),tMat1.ccap());
		CPPUNIT_ASSERT(ptr!=tMat1.begin());

		CPPUNIT_ASSERT(tMat2==tMat1);
	}

	// from cRow
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		const size_t m = genRndST(0,M-1);
		auto j = tMat2.crBegin()+m;
		
		tMat1 = *j;

		CPPUNIT_ASSERT(*j==tMat1);
	}
	
	// from cCol
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		const size_t n = genRndST(0,N-1);
		auto j = tMat2.ccBegin()+n;
		
		tMat1 = *j;

		CPPUNIT_ASSERT(*j==tMat1);
	}
}

template<class TT, class FT, class CT>
void test_tMat_assign<TT,FT,CT>::test_operator_equal_cpy() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M,N);
	const auto tMat2 = rnd<tMat>(M,N);

	const auto ptr = tMat1.begin();
	const auto lcap = tMat1.lcap();
	const auto ccap = tMat1.ccap();

	tMat1 = tMat2;

	CPPUNIT_ASSERT_EQUAL(lcap,tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(ccap,tMat1.ccap());
	CPPUNIT_ASSERT_EQUAL(ptr,tMat1.begin());

	CPPUNIT_ASSERT(tMat2==tMat1);
}

template<class TT, class FT, class CT>
void test_tMat_assign<TT,FT,CT>::test_operator_equal_move() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M,N);
	auto tMat2 = rnd<tMat>(M/2,N/2);

	TT* data;
	size_t lcap, ccap;
	{
		auto tMat3 = tMat2;
		
		data = tMat3.data();
		lcap = tMat3.lcap();
		ccap = tMat3.ccap();
		
		tMat1 = std::move(tMat3);
	}

	CPPUNIT_ASSERT_EQUAL(lcap,tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(ccap,tMat1.ccap());
	CPPUNIT_ASSERT_EQUAL(data,tMat1.data());

	CPPUNIT_ASSERT(tMat2==tMat1);
}

template<class TT, class FT, class CT>
CppUnit::Test* test_tMat_assign<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tMat_assign>(
		"test_swap", &test_tMat_assign<TT,FT,CT>::test_swap));
	suite->addTest(new CppUnit::TestCaller<test_tMat_assign>(
		"test_operator_equal_fArray", &test_tMat_assign<TT,FT,CT>::test_operator_equal_fArray));
	suite->addTest(new CppUnit::TestCaller<test_tMat_assign>(
		"test_operator_equal_cArray", &test_tMat_assign<TT,FT,CT>::test_operator_equal_cArray));
	suite->addTest(new CppUnit::TestCaller<test_tMat_assign>(
		"test_operator_equal_cpy", &test_tMat_assign<TT,FT,CT>::test_operator_equal_cpy));
	suite->addTest(new CppUnit::TestCaller<test_tMat_assign>(
		"test_operator_equal_move", &test_tMat_assign<TT,FT,CT>::test_operator_equal_move));
	
	return suite;
}
