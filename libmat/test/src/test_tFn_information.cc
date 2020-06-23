// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tFn_information.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


// tests
template<class TT, class FT, class CT>
void test_tFn_information<TT,FT,CT>::test_isnan() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = isnan(tMat1);

		CPPUNIT_ASSERT(size(tMat1)==size(tMat2));
		CPPUNIT_ASSERT(!tMat2.cpx());
		CPPUNIT_ASSERT(tMat2.logical());
		CPPUNIT_ASSERT(!any(tMat2));
	}
	{
		auto tMat1 = rnd<tMat>(M,N);
		for (size_t n=0; n!=N; ++n)
			for (size_t m=0; m!=M; ++m)
				if (m==n) tMat1(m,n) = TT(std::numeric_limits<FT>::quiet_NaN());
		const auto tMat2 = isnan(tMat1);

		CPPUNIT_ASSERT(size(tMat1)==size(tMat2));
		CPPUNIT_ASSERT(!tMat2.cpx());
		CPPUNIT_ASSERT(tMat2.logical());
		CPPUNIT_ASSERT(eye<fMat>(M,N)==tMat2);
	}
	
	{
		auto tMat1 = rnd<tMat>(M,N);
		for (size_t n=0; n!=N; ++n)
			for (size_t m=0; m!=M; ++m)
				if (m==n) tMat1(m,n) = TT(std::numeric_limits<FT>::quiet_NaN());
		
		auto r = tMat1.crBegin();
		const size_t lim = (M<N)?M:N;
		for (size_t m=0; m!=lim; ++m, ++r) {
			const auto tMat2 = isnan(*r);

			CPPUNIT_ASSERT_EQUAL(size_t(1),tMat2.M());
			CPPUNIT_ASSERT_EQUAL(N,tMat2.N());
			CPPUNIT_ASSERT(!tMat2.cpx());
			CPPUNIT_ASSERT(tMat2.logical());
			CPPUNIT_ASSERT(rId<fMat>(N,m)==tMat2);
		}
	}
	{
		auto tMat1 = rnd<tMat>(M,N);
		for (size_t n=0; n!=N; ++n)
			for (size_t m=0; m!=M; ++m)
				if (m==n) tMat1(m,n) = TT(std::numeric_limits<FT>::quiet_NaN());
		
		auto c = tMat1.ccBegin();
		const size_t lim = (N<M)?N:M;
		for (size_t n=0; n!=lim; ++n, ++c) {
			const auto tMat2 = isnan(*c);

			CPPUNIT_ASSERT_EQUAL(M,tMat2.M());
			CPPUNIT_ASSERT_EQUAL(size_t(1),tMat2.N());
			CPPUNIT_ASSERT(!tMat2.cpx());
			CPPUNIT_ASSERT(tMat2.logical());
			CPPUNIT_ASSERT(cId<fMat>(M,n)==tMat2);
		}
	}
}

template<class TT, class FT, class CT>
void test_tFn_information<TT,FT,CT>::test_length() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	const auto tMat1 = rnd<tMat>(M,N);

	const auto ref = (M>N)?M:N;
	const auto act = length(tMat1);

	CPPUNIT_ASSERT_EQUAL(ref,act);
}

template<class TT, class FT, class CT>
void test_tFn_information<TT,FT,CT>::test_rank() {
	
	size_t M = genRndST(5,10);
	size_t N = genRndST(5,10);
	if (M>N) std::swap(M,N);

	{
		tMat tMat1(M,N);
		std::fill(tMat1.begin(),tMat1.end(),rnd<tMat>(1,1,.1,.2)[0]);
		
		const auto tMat2 = rnd<tMat>(1,M,.3,.9);
		std::copy(tMat2.cbegin(),tMat2.cend(),tMat1.dbegin());
		if (rank(tMat1)!=M)

		for (size_t m=0; m!=M; ++m) {
			const size_t ref = M-m;
			const size_t act = rank(tMat1);
			CPPUNIT_ASSERT_EQUAL(ref,act);
			tMat1.rAt(m)=FT(0.0);
		}
	}
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tFn_information<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tFn_information<TT,FT,CT>>(
		"test_isnan", &test_tFn_information<TT,FT,CT>::test_isnan));
	suite->addTest(new CppUnit::TestCaller<test_tFn_information<TT,FT,CT>>(
		"test_length", &test_tFn_information<TT,FT,CT>::test_length));
	suite->addTest(new CppUnit::TestCaller<test_tFn_information<TT,FT,CT>>(
		"test_rank", &test_tFn_information<TT,FT,CT>::test_rank));
	
	return suite;
}
