// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_basic_modification.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


template<class TT, class FT, class CT>
void test_tMat_basic_modification<TT,FT,CT>::test_R_C_T() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	// R
	{
		const auto tMat1 = rnd<tMat>(M,N);
		auto tMat2 = tMat1;

		tMat2.R();
		CPPUNIT_ASSERT(tMat2.row());
		CPPUNIT_ASSERT_EQUAL(tMat2.L(),tMat1.L());
		CPPUNIT_ASSERT_EQUAL(tMat2.L(),tMat2.ccap());
		CPPUNIT_ASSERT_EQUAL(tMat2.L(),tMat2.lcap());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat2[i],tMat1[i]);
	}

	// C
	{
		const auto tMat1 = rnd<tMat>(M,N);
		auto tMat2 = tMat1;
		tMat2.C();
		CPPUNIT_ASSERT(tMat2.col());
		CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat2.L());
		CPPUNIT_ASSERT_EQUAL(size_t(1),tMat2.ccap());
		CPPUNIT_ASSERT_EQUAL(tMat2.L(),tMat2.lcap());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat2[i],tMat1[i]);
	}

	// T
	{
		const size_t M = genRndST();
		const auto tMat1 = rnd<tMat>(M,1);
		auto tMat2=tMat1;
		tMat2.T();
			
		CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.N());
		CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat2.M());
			
		for (size_t n=0; n!=tMat1.N(); ++n)
			for (size_t m=0; m!=tMat1.M(); ++m)
				CPPUNIT_ASSERT_EQUAL(CT(tMat1(m,n)),CT(std::conj(tMat2(n,m))));
	

		const auto tMat3 = rnd<tMat>(M,M);
		auto tMat4=tMat3;
		
		tMat4.T();
		for (size_t n=0; n!=tMat1.N(); ++n)
			for (size_t m=0; m!=tMat1.M(); ++m)
				CPPUNIT_ASSERT_EQUAL(CT(tMat3(m,n)),CT(std::conj(tMat4(n,m))));

		const auto S = genRndMEST();
		const auto tMat5 = rnd<tMat>(S.M,S.N);
		auto tMat6=tMat5;
		tMat6.T();
			
		CPPUNIT_ASSERT_EQUAL(tMat6.M(),tMat5.N());
		CPPUNIT_ASSERT_EQUAL(tMat6.N(),tMat5.M());
			
		for (size_t n=0; n!=tMat5.N(); ++n)
			for (size_t m=0; m!=tMat5.M(); ++m)
				CPPUNIT_ASSERT_EQUAL(CT(tMat5(m,n)),CT(std::conj(tMat6(n,m))));
	}
}

template<class TT, class FT, class CT>
void test_tMat_basic_modification<TT,FT,CT>::test_clear_resize() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
		
	auto tMat1 = rnd<tMat>(M,N);
	const auto tMat2 = tMat1;

	// clear
	{
		auto tMat3 = tMat1;
		tMat3.clear();
		CPPUNIT_ASSERT(tMat3.empty());
		CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat3.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),tMat3.N());
		CPPUNIT_ASSERT_EQUAL(tMat1.lcap(),tMat3.lcap());
	}

	// resize
	{
		const auto pre = tMat1.begin();
		tMat1.resize(2*N);
		const auto post = tMat1.begin();
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(2*N,tMat1.N());
		CPPUNIT_ASSERT_EQUAL(2*N,tMat1.ccap());
		CPPUNIT_ASSERT_EQUAL(2*M*N,tMat1.lcap());
		CPPUNIT_ASSERT(pre!=post);
		CPPUNIT_ASSERT(tMat2==tMat1.get(0,0,M,N));
	}
	{
		const auto pre = tMat1.begin();
		tMat1.resize(N);
		const auto post = tMat1.begin();
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
		CPPUNIT_ASSERT_EQUAL(2*N,tMat1.ccap());
		CPPUNIT_ASSERT_EQUAL(2*M*N,tMat1.lcap());
		CPPUNIT_ASSERT_EQUAL(pre,post);
		CPPUNIT_ASSERT(tMat1==tMat2);
	}
}

template<class TT, class FT, class CT>
void test_tMat_basic_modification<TT,FT,CT>::test_push_back_shift() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// pushBack
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = tMat1;
		const auto tMat3 = rnd<tMat>(M,1);

		{
			const auto pre = tMat1.begin();
			tMat1.push_back(tMat3);
			const auto post = tMat1.begin();
			
			CPPUNIT_ASSERT(pre!=post);
			CPPUNIT_ASSERT(tMat1==ncat(tMat2,tMat3));
			CPPUNIT_ASSERT_EQUAL(N+1,tMat1.ccap());
			CPPUNIT_ASSERT_EQUAL((N+1)*M,tMat1.lcap());
		}
		{
			tMat1.reserve(N+2);
				
			const auto pre = tMat1.begin();
			tMat1.push_back(tMat3);
			const auto post = tMat1.begin();

			CPPUNIT_ASSERT(pre==post);
			CPPUNIT_ASSERT(tMat1==ncat(ncat(tMat2,tMat3),tMat3));
			CPPUNIT_ASSERT_EQUAL(N+2,tMat1.ccap());
			CPPUNIT_ASSERT_EQUAL((N+2)*M,tMat1.lcap());
		}
		{
			tMat tMat4(M,0); tMat4.reserve(N);
			for (size_t i=0; i<tMat1.N(); ++i)
				tMat4.push_back(tMat1.cGet(i));
			CPPUNIT_ASSERT(tMat4==tMat1);
		}
		{
			auto tMat4 = rnd<tMat>(M,N);
			const auto tMat5 = tMat4;
			tMat4.push_back(tMat1);
			CPPUNIT_ASSERT(tMat4==ncat(tMat5,tMat1));
		}
	}

	// operator<<
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<tMat>(M,N/2);

		{
			auto tMat3 = tMat1;
			const auto pre = tMat3.begin();

			tMat3 << tMat2;
			const auto post = tMat3.begin();

			CPPUNIT_ASSERT(pre!=post);
			CPPUNIT_ASSERT(tMat3==ncat(tMat1,tMat2));
			CPPUNIT_ASSERT_EQUAL(tMat1.N()+tMat2.N(),tMat3.ccap());
			CPPUNIT_ASSERT_EQUAL((tMat1.N()+tMat2.N())*M,tMat3.lcap());
		}
		{
			auto tMat3 = tMat1;
			tMat3.reserve(tMat1.N()+tMat2.N());
			const auto pre = tMat3.begin();
	
			tMat3 << tMat2;
			const auto post = tMat3.begin();

			CPPUNIT_ASSERT(pre==post);
			CPPUNIT_ASSERT(tMat3==ncat(tMat1,tMat2));
			CPPUNIT_ASSERT_EQUAL(tMat1.N()+tMat2.N(),tMat3.ccap());
			CPPUNIT_ASSERT_EQUAL((tMat1.N()+tMat2.N())*M,tMat3.lcap());
		}
		{
			auto tMat3 = tMat1;
			
			tMat3 << (tMat3.N()+2*tMat2.N()) << tMat2 << tMat2;
			
			CPPUNIT_ASSERT(tMat3==ncat(tMat1,ncat(tMat2,tMat2)));
			CPPUNIT_ASSERT_EQUAL(tMat1.N()+2*tMat2.N(),tMat3.ccap());
			CPPUNIT_ASSERT_EQUAL((tMat1.N()+2*tMat2.N())*M,tMat3.lcap());
		}
	}
}

template<class TT, class FT, class CT>
void test_tMat_basic_modification<TT,FT,CT>::test_pop_back() {
	const size_t M = genRndST();
	const size_t N = genRndST();
		
	auto tMat1 = rnd<tMat>(M,N);
	const auto tMat2 = tMat1;

	const auto pre = tMat1.begin();
	tMat1.pop_back();
	const auto post = tMat1.begin();
		
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(N-1,tMat1.N());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.ccap());
	CPPUNIT_ASSERT_EQUAL(M*N,tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(pre,post);
	CPPUNIT_ASSERT(tMat1==tMat2.get(0,0,M,N-1));
}

template<class TT, class FT, class CT>
CppUnit::Test* test_tMat_basic_modification<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tMat_basic_modification>(
		"test_R_C_T", &test_tMat_basic_modification<TT,FT,CT>::test_R_C_T));
	suite->addTest(new CppUnit::TestCaller<test_tMat_basic_modification>(
		"test_clear_resize", &test_tMat_basic_modification<TT,FT,CT>::test_clear_resize));
	suite->addTest(new CppUnit::TestCaller<test_tMat_basic_modification>(
		"test_push_back_shift", &test_tMat_basic_modification<TT,FT,CT>::test_push_back_shift));
	suite->addTest(new CppUnit::TestCaller<test_tMat_basic_modification>(
		"test_pop_back", &test_tMat_basic_modification<TT,FT,CT>::test_pop_back));
	
	return suite;
}
