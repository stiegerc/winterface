// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tFn_mat_gen.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;



template<class TT, class FT, class CT>
void test_tFn_mat_gen<TT,FT,CT>::test_ones() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		const auto tMat1 = ones<fMat>(0,0);
		CPPUNIT_ASSERT(tMat1.empty());
	}
	{
		const auto tMat1 = ones<tMat>(M);
		for (auto& i: tMat1)
			CPPUNIT_ASSERT(i==TT(1.0));
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(M,tMat1.N());
	}
	{
		const auto tMat1 = ones<tMat>(M,N);
		for (auto& i: tMat1)
			CPPUNIT_ASSERT(i==TT(1.0));
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
	}
}

template<class TT, class FT, class CT>
void test_tFn_mat_gen<TT,FT,CT>::test_zeros() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		const auto tMat1 = zeros<fMat>(0,0);
		CPPUNIT_ASSERT(tMat1.empty());
	}
	{
		const auto tMat1 = zeros<tMat>(M);
		for (auto& i: tMat1)
			CPPUNIT_ASSERT(i==TT(0.0));
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(M,tMat1.N());
	}
	{
		const auto tMat1 = zeros<tMat>(M,N);
		for (auto& i: tMat1)
			CPPUNIT_ASSERT(i==TT(0.0));
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
	}
}

template<class TT, class FT, class CT>
void test_tFn_mat_gen<TT,FT,CT>::test_eye() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		const auto tMat1 = eye<fMat>(0,0);
		CPPUNIT_ASSERT(tMat1.empty());
	}
	{
		const auto tMat1 = eye<tMat>(M);
		for (size_t n=0; n!=tMat1.N(); ++n)
			for (size_t m=0; m!=tMat1.M(); ++m)
				CPPUNIT_ASSERT_EQUAL(m==n ? TT(1.0): TT(0.0),tMat1(m,n));
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(M,tMat1.N());
	}
	{
		const auto tMat1 = eye<tMat>(M,N);
		for (size_t n=0; n!=tMat1.N(); ++n)
			for (size_t m=0; m!=tMat1.M(); ++m)
				CPPUNIT_ASSERT_EQUAL(m==n ? TT(1.0): TT(0.0),tMat1(m,n));
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
	}
}

template<class TT, class FT, class CT>
void test_tFn_mat_gen<TT,FT,CT>::test_rId_cId() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	for (size_t i=0; i!=N; ++i) {
		const auto cr = rId<tMat>(N,i);
		CPPUNIT_ASSERT(cr.row());
		CPPUNIT_ASSERT(cr.L()==N);
		for (size_t j=0; j!=N; ++j)
			if (j==i) CPPUNIT_ASSERT_EQUAL(TT(1.0),cr[j]);
			else CPPUNIT_ASSERT_EQUAL(TT(0.0),cr[j]);
	}
	for (size_t i=0; i!=M; ++i) {
		const auto cc = cId<tMat>(M,i);
		CPPUNIT_ASSERT(cc.col());
		CPPUNIT_ASSERT(cc.L()==M);
		for (size_t j=0; j!=M; ++j)
			if (j==i) CPPUNIT_ASSERT_EQUAL(TT(1.0),cc[j]);
			else CPPUNIT_ASSERT_EQUAL(TT(0.0),cc[j]);
	}
}

template<class TT, class FT, class CT>
void test_tFn_mat_gen<TT,FT,CT>::test_rand() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		const auto tMat1 = rand<fMat>(0,0);
		CPPUNIT_ASSERT(tMat1.empty());
	}
	{
		const auto tMat1 = rand<tMat>(M);
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(M,tMat1.N());
	}
	{
		const RE__ l = 0.0;
		const RE__ u = 1.0;
		const auto tMat1 = rand<tMat>(M,N,l,u);
		
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
		for (const auto i: tMat1)
			CPPUNIT_ASSERT(i>=l && i<=u);
	}
	{
		const RE__ ll_ = -1.0;
		const RE__ ul_ = 0.0;
		const RE__ lu_ = 1.0;
		const RE__ uu_ = 2.0;
	
		const auto l = rand<fMat>(M,N,ll_,ul_);
		const auto u = rand<fMat>(M,N,lu_,uu_);
		
		const auto tMat1 = rand<tMat>(l,u);
		
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT(tMat1[i]>=l[i] && tMat1[i]<=u[i]);
	}
}

template<class TT, class FT, class CT>
void test_tFn_mat_gen<TT,FT,CT>::test_randi() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = randi<fMat>(0,0);
		CPPUNIT_ASSERT(tMat1.empty());
	}
	{
		const auto tMat1 = randi<tMat>(M);
		
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(M,tMat1.N());
		for (const auto i: tMat1)
			CPPUNIT_ASSERT_EQUAL(std::round(i),i);
	}
	{
		const long l=0;
		const long u=100;
		const auto tMat1 = randi<tMat>(M,N);
		
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
		for (const auto i: tMat1) {
			CPPUNIT_ASSERT_EQUAL(std::round(i),i);
			CPPUNIT_ASSERT(i>=l && i<=u);
		}
	}
	{
		const long ll_ = -100;
		const long ul_ = 0;
		const long lu_ = 1;
		const long uu_ = 100;
	
		const auto l = randi<fMat>(M,N,ll_,ul_);
		const auto u = randi<fMat>(M,N,lu_,uu_);
		
		const auto tMat1 = randi<tMat>(l,u);
		
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
		for (size_t i=0; i!=tMat1.L(); ++i) {
			CPPUNIT_ASSERT_EQUAL(std::round(tMat1[i]),tMat1[i]);
			CPPUNIT_ASSERT(tMat1[i]>=l[i] && tMat1[i]<=u[i]);
		}
	}
}

template<class TT, class FT, class CT>
void test_tFn_mat_gen<TT,FT,CT>::test_gsorth() {
	
	const size_t M = genRndST();

	auto tMat1 = rnd_b<tMat>(M,M);
	gsorth(tMat1);
	
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(M,tMat1.N());

	CPPUNIT_ASSERT_DELTA(std::abs(det(tMat1)),RE__(1.0),delta);
}

template<class TT, class FT, class CT>
void test_tFn_mat_gen<TT,FT,CT>::test_mcat() {
	
	const size_t M1 = genRndST();
	const size_t M2 = genRndST();
	const size_t N = genRndST();

	const auto tMat1 = rnd<tMat>(M1,N);
	const auto tMat2 = rnd<tMat>(M2,N);
	
	{
		const auto tMat3 = mcat(tMat1,tMat2);
		CPPUNIT_ASSERT_EQUAL(M1+M2,tMat3.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat3.N());

		for (size_t n=0; n!=N; ++n) {
			for (size_t m=0; m!=M1; ++m)
				CPPUNIT_ASSERT_EQUAL(tMat1(m,n),tMat3(m,n));

			for (size_t m=0; m!=M2; ++m)
				CPPUNIT_ASSERT_EQUAL(tMat2(m,n),tMat3(m+M1,n));
		}
	}
	{
		const auto tMat3 = mcat(tMat1,*tMat2.crBegin());
		CPPUNIT_ASSERT_EQUAL(M1+1,tMat3.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat3.N());

		for (size_t n=0; n!=N; ++n) {
			for (size_t m=0; m!=M1; ++m)
				CPPUNIT_ASSERT_EQUAL(tMat1(m,n),tMat3(m,n));

			CPPUNIT_ASSERT_EQUAL(tMat2(0,n),tMat3(M1,n));
		}
	}
	{
		const auto tMat3 = mcat(*tMat1.crBegin(),tMat2);
		CPPUNIT_ASSERT_EQUAL(M2+1,tMat3.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat3.N());

		for (size_t n=0; n!=N; ++n) {
			CPPUNIT_ASSERT_EQUAL(tMat1(0,n),tMat3(0,n));
			
			for (size_t m=0; m!=M2; ++m)
				CPPUNIT_ASSERT_EQUAL(tMat2(m,n),tMat3(m+1,n));
		}
	}
	{
		const auto tMat3 = mcat(*tMat1.crBegin(),*tMat2.crBegin());
		CPPUNIT_ASSERT_EQUAL(size_t(2),tMat3.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat3.N());

		for (size_t n=0; n!=N; ++n) {
			CPPUNIT_ASSERT_EQUAL(tMat1(0,n),tMat3(0,n));
			CPPUNIT_ASSERT_EQUAL(tMat2(0,n),tMat3(1,n));
		}
	}
}

template<class TT, class FT, class CT>
void test_tFn_mat_gen<TT,FT,CT>::test_ncat() {
	
	const size_t M = genRndST();
	const size_t N1 = genRndST();
	const size_t N2 = genRndST();

	const auto tMat1 = rnd<tMat>(M,N1);
	const auto tMat2 = rnd<tMat>(M,N2);
	
	{
		const auto tMat3 = ncat(tMat1,tMat2);
		CPPUNIT_ASSERT_EQUAL(M,tMat3.M());
		CPPUNIT_ASSERT_EQUAL(N1+N2,tMat3.N());

		for (size_t m=0; m!=M; ++m) {
			for (size_t n=0; n!=N1; ++n)
				CPPUNIT_ASSERT_EQUAL(tMat1(m,n),tMat3(m,n));

			for (size_t n=0; n!=N2; ++n)
				CPPUNIT_ASSERT_EQUAL(tMat2(m,n),tMat3(m,n+N1));
		}
	}
	{
		const auto tMat3 = ncat(tMat1,*tMat2.ccBegin());
		CPPUNIT_ASSERT_EQUAL(M,tMat3.M());
		CPPUNIT_ASSERT_EQUAL(N1+1,tMat3.N());

		for (size_t m=0; m!=M; ++m) {
			for (size_t n=0; n!=N1; ++n)
				CPPUNIT_ASSERT_EQUAL(tMat1(m,n),tMat3(m,n));

			CPPUNIT_ASSERT_EQUAL(tMat2(m,0),tMat3(m,N1));
		}
	}
	{
		const auto tMat3 = ncat(*tMat1.ccBegin(),tMat2);
		CPPUNIT_ASSERT_EQUAL(M,tMat3.M());
		CPPUNIT_ASSERT_EQUAL(N2+1,tMat3.N());

		for (size_t m=0; m!=M; ++m) {
			CPPUNIT_ASSERT_EQUAL(tMat1(m,0),tMat3(m,0));
			
			for (size_t n=0; n!=N2; ++n)
				CPPUNIT_ASSERT_EQUAL(tMat2(m,n),tMat3(m,n+1));
		}
	}
	{
		const auto tMat3 = ncat(*tMat1.ccBegin(),*tMat2.ccBegin());
		CPPUNIT_ASSERT_EQUAL(M,tMat3.M());
		CPPUNIT_ASSERT_EQUAL(size_t(2),tMat3.N());

		for (size_t m=0; m!=M; ++m) {
			CPPUNIT_ASSERT_EQUAL(tMat1(m,0),tMat3(m,0));
			CPPUNIT_ASSERT_EQUAL(tMat2(m,0),tMat3(m,1));
		}
	}
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tFn_mat_gen<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tFn_mat_gen>(
		"test_ones", &test_tFn_mat_gen<TT,FT,CT>::test_ones));
	suite->addTest(new CppUnit::TestCaller<test_tFn_mat_gen>(
		"test_zeros", &test_tFn_mat_gen<TT,FT,CT>::test_zeros));
	suite->addTest(new CppUnit::TestCaller<test_tFn_mat_gen>(
		"test_eye", &test_tFn_mat_gen<TT,FT,CT>::test_eye));
	suite->addTest(new CppUnit::TestCaller<test_tFn_mat_gen>(
		"test_rId_cId", &test_tFn_mat_gen<TT,FT,CT>::test_rId_cId));
	suite->addTest(new CppUnit::TestCaller<test_tFn_mat_gen>(
		"test_rand", &test_tFn_mat_gen<TT,FT,CT>::test_rand));
	suite->addTest(new CppUnit::TestCaller<test_tFn_mat_gen>(
		"test_randi", &test_tFn_mat_gen<TT,FT,CT>::test_randi));
	suite->addTest(new CppUnit::TestCaller<test_tFn_mat_gen>(
		"test_gsorth", &test_tFn_mat_gen<TT,FT,CT>::test_gsorth));
	suite->addTest(new CppUnit::TestCaller<test_tFn_mat_gen>(
		"test_mcat", &test_tFn_mat_gen<TT,FT,CT>::test_mcat));
	suite->addTest(new CppUnit::TestCaller<test_tFn_mat_gen>(
		"test_ncat", &test_tFn_mat_gen<TT,FT,CT>::test_ncat));
	
	return suite;
}
