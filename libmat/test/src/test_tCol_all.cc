// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tCol_all.h"
#include "testTools.h"
#include "lm_tMat.h"
#include "lm_lambda.h"
#include <iostream>
#include <type_traits>

using namespace lm__;
using namespace lm__::test;

template<class TT, class FT, class CT>
void test_tCol_all<TT,FT,CT>::test_swap() {
	
	const size_t M = genRndST();
	const size_t N1 = genRndST();
	const size_t N2 = genRndST();

	auto tMat1 = rnd<tMat>(M,N1);
	auto tMat2 = rnd<tMat>(M,N2);

	const size_t n1 = genRndST(0,N1-1);
	const size_t n2 = genRndST(0,N2-1);

	const auto tMat3 = tMat1;
	const auto tMat4 = tMat2;

	swap(tMat1.cAt(n1),tMat2.cAt(n2));

	for (size_t m=0; m!=M; ++m)
		CPPUNIT_ASSERT_EQUAL(tMat4(m,n2),tMat1(m,n1));
	for (size_t m=0; m!=M; ++m)
		CPPUNIT_ASSERT_EQUAL(tMat3(m,n1),tMat2(m,n2));
}

template<class TT, class FT, class CT>
void test_tCol_all<TT,FT,CT>::test_data_access() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		auto tMat1 = rnd<tMat>(M,N);
		const size_t n = genRndST(0,N-1);
	
		auto j = tMat1.cBegin()+n;
		CPPUNIT_ASSERT_EQUAL(tMat1.data()+n*M,j->data());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const size_t n = genRndST(0,N-1);
	
		auto j = tMat1.cBegin()+n;
		CPPUNIT_ASSERT_EQUAL(tMat1.data()+n*M,j->data());
	}
}

template<class TT, class FT, class CT>
void test_tCol_all<TT,FT,CT>::test_information() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M,N);
	
	const size_t n = genRndST(0,N-1);
	auto j = tMat1.cBegin()+n;

	CPPUNIT_ASSERT_EQUAL(ptrdiff_t(n),j->i());
	CPPUNIT_ASSERT_EQUAL(ptrdiff_t(N),j->iL());
	CPPUNIT_ASSERT_EQUAL(M,j->M());
	CPPUNIT_ASSERT_EQUAL(size_t(1),j->N());
	CPPUNIT_ASSERT_EQUAL(size_t(1),j->incr());
}

template<class TT, class FT, class CT>
void test_tCol_all<TT,FT,CT>::test_conversion() {
	const size_t M = genRndST();
	
	auto tMat1 = rnd<tMat>(M,1);
	auto j = tMat1.cBegin();
	
	CPPUNIT_ASSERT(tMat1==j->copy());

	fMat tMat2(tMat1);
	CPPUNIT_ASSERT(tMat2==j->fcopy());
	
	cMat tMat3(tMat1);
	CPPUNIT_ASSERT(tMat3==j->ccopy());
}

template<class TT, class FT, class CT>
void test_tCol_all<TT,FT,CT>::test_copy_algo_() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M,N);
	auto ptr = tMat1.ptr();

	// copy holds own data, temporary does not
	{
		// temporary
		auto i = tMat1.cBegin();
		CPPUNIT_ASSERT_EQUAL(ptr,i->ptr());
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),i->i());

		// copy
		const auto c = *i;
		CPPUNIT_ASSERT(c.ptr()!=ptr);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(PTRDIFF_MAX),c.i());

		// overwrite and compare
		const auto tMat2 = tMat1;

		*i = *(i+1);
		CPPUNIT_ASSERT(c!=*i);

		*i = c;
		CPPUNIT_ASSERT(tMat2==tMat1);
	}

	// try sorting
	{
		std::sort(tMat1.cBegin(),tMat1.cEnd(),vcmp);
		CPPUNIT_ASSERT(std::is_sorted(tMat1.ccBegin(),tMat1.ccEnd(),vcmp));
	}

	// sorting on integer matrix
	{
		auto tMat2 = round(rnd<tMat>(M,N,FT(-1.0),FT(1.0)));
		std::sort(tMat2.cBegin(),tMat2.cEnd(),vcmp);
		std::is_sorted(tMat2.ccBegin(),tMat2.ccEnd(),vcmp);
	}
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tCol_all<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tCol_all>(
		"test_assign", &test_tCol_all<TT,FT,CT>::test_assign));
	suite->addTest(new CppUnit::TestCaller<test_tCol_all>(
		"test_swap", &test_tCol_all<TT,FT,CT>::test_swap));
	suite->addTest(new CppUnit::TestCaller<test_tCol_all>(
		"test_data_access", &test_tCol_all<TT,FT,CT>::test_data_access));
	suite->addTest(new CppUnit::TestCaller<test_tCol_all>(
		"test_information", &test_tCol_all<TT,FT,CT>::test_information));
	suite->addTest(new CppUnit::TestCaller<test_tCol_all>(
		"test_conversion", &test_tCol_all<TT,FT,CT>::test_conversion));
	suite->addTest(new CppUnit::TestCaller<test_tCol_all>(
		"test_matrix_arithmetic", &test_tCol_all<TT,FT,CT>::test_matrix_arithmetic));
	suite->addTest(new CppUnit::TestCaller<test_tCol_all>(
		"test_copy_algo_", &test_tCol_all<TT,FT,CT>::test_copy_algo_));
	
	return suite;
}
