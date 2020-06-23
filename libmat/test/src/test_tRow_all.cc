// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tRow_all.h"
#include "testTools.h"
#include "lm_tMat.h"
#include "lm_lambda.h"
#include <iostream>


using namespace lm__;
using namespace lm__::test;

template<class TT, class FT, class CT>
void test_tRow_all<TT,FT,CT>::test_swap() {
	
	const size_t M1 = genRndST();
	const size_t M2 = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M1,N);
	auto tMat2 = rnd<tMat>(M2,N);

	const size_t m1 = genRndST(0,M1-1);
	const size_t m2 = genRndST(0,M2-1);

	const auto tMat3 = tMat1;
	const auto tMat4 = tMat2;

	swap(tMat1.rAt(m1),tMat2.rAt(m2));

	for (size_t n=0; n!=N; ++n)
		CPPUNIT_ASSERT_EQUAL(tMat4(m2,n),tMat1(m1,n));
	for (size_t n=0; n!=N; ++n)
		CPPUNIT_ASSERT_EQUAL(tMat3(m1,n),tMat2(m2,n));
}

template<class TT, class FT, class CT>
void test_tRow_all<TT,FT,CT>::test_data_access() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		auto tMat1 = rnd<tMat>(M,N);
		const size_t m = genRndST(0,M-1);
	
		auto j = tMat1.rBegin()+m;
		CPPUNIT_ASSERT_EQUAL(tMat1.data()+m,j->data());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const size_t m = genRndST(0,M-1);
	
		auto j = tMat1.rBegin()+m;
		CPPUNIT_ASSERT_EQUAL(tMat1.data()+m,j->data());
	}
}

template<class TT, class FT, class CT>
void test_tRow_all<TT,FT,CT>::test_information() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M,N);
	
	const size_t m = genRndST(0,M-1);
	auto j = tMat1.rBegin()+m;

	CPPUNIT_ASSERT_EQUAL(ptrdiff_t(m),j->i());
	CPPUNIT_ASSERT_EQUAL(ptrdiff_t(M),j->iL());
	CPPUNIT_ASSERT_EQUAL(N,j->N());
	CPPUNIT_ASSERT_EQUAL(size_t(1),j->M());
	CPPUNIT_ASSERT_EQUAL(M,j->incr());
}

template<class TT, class FT, class CT>
void test_tRow_all<TT,FT,CT>::test_conversion() {
	const size_t N = genRndST();
	
	auto tMat1 = rnd<tMat>(1,N);
	auto j = tMat1.rBegin();
	
	CPPUNIT_ASSERT(tMat1==j->copy());

	fMat tMat2(tMat1);
	CPPUNIT_ASSERT(tMat2==j->fcopy());
	
	cMat tMat3(tMat1);
	CPPUNIT_ASSERT(tMat3==j->ccopy());
}

template<class TT, class FT, class CT>
void test_tRow_all<TT,FT,CT>::test_copy_algo_() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M,N);
	auto ptr = tMat1.ptr();

	// copy holds own data, temporary does not
	{
		// temporary
		auto i = tMat1.rBegin();
		CPPUNIT_ASSERT_EQUAL(ptr,i->ptr());
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),i->i());

		// copy
		const auto r = *i;
		CPPUNIT_ASSERT(r.ptr()!=ptr);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(PTRDIFF_MAX),r.i());

		// overwrite and compare
		const auto tMat2 = tMat1;

		*i = *(i+1);
		CPPUNIT_ASSERT(r!=*i);

		*i = r;
		CPPUNIT_ASSERT(tMat2==tMat1);
	}

	// try sorting
	{
		std::sort(tMat1.rBegin(),tMat1.rEnd(),vcmp);
		CPPUNIT_ASSERT(std::is_sorted(tMat1.crBegin(),tMat1.crEnd(),vcmp));
	}

	// sorting on integer matrix
	{
		auto tMat2 = round(rnd<tMat>(M,N,FT(-1.0),FT(1.0)));
		std::sort(tMat2.rBegin(),tMat2.rEnd(),vcmp);
		std::is_sorted(tMat2.crBegin(),tMat2.crEnd(),vcmp);
	}
}

template<class TT, class FT, class CT>
CppUnit::Test* test_tRow_all<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tRow_all>(
		"test_assign", &test_tRow_all<TT,FT,CT>::test_assign));
	suite->addTest(new CppUnit::TestCaller<test_tRow_all>(
		"test_swap", &test_tRow_all<TT,FT,CT>::test_swap));
	suite->addTest(new CppUnit::TestCaller<test_tRow_all>(
		"test_data_access", &test_tRow_all<TT,FT,CT>::test_data_access));
	suite->addTest(new CppUnit::TestCaller<test_tRow_all>(
		"test_information", &test_tRow_all<TT,FT,CT>::test_information));
	suite->addTest(new CppUnit::TestCaller<test_tRow_all>(
		"test_conversion", &test_tRow_all<TT,FT,CT>::test_conversion));
	suite->addTest(new CppUnit::TestCaller<test_tRow_all>(
		"test_matrix_arithmetic", &test_tRow_all<TT,FT,CT>::test_matrix_arithmetic));
	suite->addTest(new CppUnit::TestCaller<test_tRow_all>(
		"test_copy_algo_", &test_tRow_all<TT,FT,CT>::test_copy_algo_));
	
	return suite;
}
