// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_conversion.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


template<class TT, class FT, class CT>
void test_tMat_conversion<TT,FT,CT>::test_copy_fcopy_ccopy() {

	const size_t M = genRndST();
	const size_t N = genRndST();
	const auto tMat1 = rnd<tMat>(M,N);

	// copy
	{
		const tMat tMat2(tMat1);
		CPPUNIT_ASSERT(tMat2==tMat1.copy());
	}
	
	// fcopy
	{
		const fMat tMat2(tMat1);
		CPPUNIT_ASSERT(tMat2==tMat1.fcopy());
	}

	// ccopy
	{
		const cMat tMat2(tMat1);
		CPPUNIT_ASSERT(tMat2==tMat1.ccopy());
	}
}

template<class TT, class FT, class CT>
void test_tMat_conversion<TT,FT,CT>::test_get_getl() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	const auto tMat1 = rnd<tMat>(M,N);
	
	const auto m = genRndI(1,M/4,0,M/2);
	const auto n = genRndI(1,N/4,0,N/2);
				
	// generic case
	{
		auto tMat2 = tMat1.get(m,n);
		CPPUNIT_ASSERT_EQUAL(m.size(),tMat2.M());
		CPPUNIT_ASSERT_EQUAL(n.size(),tMat2.N());
		for (size_t nv=0; nv<n.size(); ++nv)
			for (size_t mv=0; mv<m.size(); ++mv)
				CPPUNIT_ASSERT_EQUAL(tMat1(m[mv],n[nv]),tMat2(mv,nv));
	}
	
	// full m, partial n
	{
		auto tMat2 = tMat1.get({},n);
		CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.M());
		CPPUNIT_ASSERT_EQUAL(n.size(),tMat2.N());
		for (size_t nv=0; nv<n.size(); ++nv)
			for (size_t mv=0; mv<tMat1.M(); ++mv)
				CPPUNIT_ASSERT_EQUAL(tMat1(mv,n[nv]),tMat2(mv,nv));
	}
	
	// full n, partial m
	{
		auto tMat2 = tMat1.get(m,{});
		CPPUNIT_ASSERT_EQUAL(m.size(),tMat2.M());
		CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat2.N());
		for (size_t nv=0; nv<tMat1.N(); ++nv)
			for (size_t mv=0; mv<m.size(); ++mv)
				CPPUNIT_ASSERT_EQUAL(tMat1(m[mv],nv),tMat2(mv,nv));
	}
	
	// full m, full n
	{
		auto tMat2 = tMat1.get({},{});
		CPPUNIT_ASSERT(tMat1==tMat2);
	}

	// 4 size_t input version
	{
		const size_t m = genRndST(0,M-1);
		const size_t n = genRndST(0,N-1);
		const size_t lm = genRndST(0,M);
		const size_t ln = genRndST(0,N);
			
		{
			const auto tMat2 = tMat1.get(m,n);
			CPPUNIT_ASSERT_EQUAL(tMat1.M()-m,tMat2.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.N()-n,tMat2.N());
			for (size_t i=m; i<tMat1.M(); ++i)
				for (size_t j=n; j<tMat1.N(); ++j)
					CPPUNIT_ASSERT_EQUAL(tMat2(i-m,j-n),tMat1(i,j));
		}
		{
			const auto tMat2 = tMat1.get(m,n,lm,ln);
			CPPUNIT_ASSERT(tMat2.M()==(tMat1.M()-m<lm ? tMat1.M()-m: lm));
			CPPUNIT_ASSERT(tMat2.N()==(tMat1.N()-n<ln ? tMat1.N()-n: ln));
			for (size_t i=0; i<tMat2.M(); ++i)
				for (size_t j=0; j<tMat2.N(); ++j)
					CPPUNIT_ASSERT_EQUAL(tMat2(i,j),tMat1(m+i,n+j));
		}
	}

	// getl (logical indexing version)
	{
		// matrix version
		{
			const fMat I = tMat1.lt(mean(tMat1));
			const auto tMat2 = tMat1.getl(I);

			const size_t cnt = sum(I);
			CPPUNIT_ASSERT_EQUAL(cnt,tMat2.M());
			CPPUNIT_ASSERT_EQUAL(size_t(1),tMat2.N());

			auto j = tMat2.cbegin();
			for (size_t i=0; i!=tMat2.L(); ++i)
				if (I[i]) CPPUNIT_ASSERT_EQUAL(*j++,tMat1[i]);
		}

		// col, row version
		{
			const auto ml = nsum(tMat1).lt(mean(nsum(tMat1)));
			const auto nl = msum(tMat1).lt(mean(msum(tMat1)));
				
			std::vector<size_t> mlck; mlck.reserve(size_t(sum(ml)));
			for (size_t i=0; i<ml.L(); ++i)
				if (bool(ml[i])) mlck.push_back(i);
			
			std::vector<size_t> nlck; nlck.reserve(size_t(sum(nl)));
			for (size_t i=0; i<nl.L(); ++i)
				if (bool(nl[i])) nlck.push_back(i);
			
			CPPUNIT_ASSERT(tMat1.getl(ml,nl)==tMat1.get(mlck,nlck));
		}
	}
}

template<class TT, class FT, class CT>
void test_tMat_conversion<TT,FT,CT>::test_rGet_cGet() {
	// rGet
	{
		const auto tMat1 = rnd<tMat>(genRndST(),genRndST());
	
		// default 1 row	
		for (size_t m=0; m!=tMat1.M(); ++m) {
			const auto row = tMat1.rGet(m);
			CPPUNIT_ASSERT_EQUAL(size_t(1),row.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.N(),row.N());
			for (size_t n=0; n!=tMat1.N(); ++n)
				CPPUNIT_ASSERT_EQUAL(tMat1(m,n),row[n]);
		}

		// random many rows
		for (size_t m=0; m!=tMat1.M(); ++m) {
			const size_t l = genRndST(0,tMat1.M()-m);
			const auto row = tMat1.rGet(m,l);
			for (size_t n=0; n!=tMat1.N(); ++n)
				for (size_t cl=0; cl!=l; ++cl)
					CPPUNIT_ASSERT_EQUAL(tMat1(m+cl,n),row(cl,n));
		}
	}
	
	// cGet
	{
		const auto tMat1 = rnd<tMat>(genRndST(),genRndST());

		// default 1 col
		for (size_t n=0; n!=tMat1.N(); ++n) {
			auto col = tMat1.cGet(n);
			CPPUNIT_ASSERT_EQUAL(tMat1.M(),col.M());
			CPPUNIT_ASSERT_EQUAL(size_t(1),col.N());
			for (size_t m=0; m!=tMat1.M(); ++m)
				CPPUNIT_ASSERT_EQUAL(tMat1(m,n),col[m]);
		}
	
		// random many rows
		for (size_t n=0; n!=tMat1.N(); ++n) {
			const size_t l = genRndST(0,tMat1.N()-n);
			const auto col = tMat1.cGet(n,l);
			for (size_t m=0; m!=tMat1.M(); ++m)
				for (size_t cl=0; cl!=l; ++cl)
					CPPUNIT_ASSERT_EQUAL(tMat1(m,n+cl),col(m,cl));
		}
	}
}

template<class TT, class FT, class CT>
void test_tMat_conversion<TT,FT,CT>::test_rWOGet_cWOGet() {
	
	// rWOGet
	{
		const auto tMat1 = rnd<tMat>(genRndST(),genRndST());
		
		for (size_t i=0; i<tMat1.M(); ++i) {
			const auto tMat2 = tMat1.rWOGet(i);
				
			CPPUNIT_ASSERT_EQUAL(tMat1.M()-1,tMat2.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat2.N());
			CPPUNIT_ASSERT_EQUAL(tMat1.L()-tMat1.N(),tMat2.L());
			for (size_t m=0; m<tMat1.M(); ++m) {
				if (m==i) continue;
				for (size_t n=0; n<tMat1.N(); ++n)
					CPPUNIT_ASSERT_EQUAL(tMat1(m,n),tMat2(m-(m>i),n));
			}
		}
	}
	
	// cWOGet
	{
		const auto tMat1 = rnd<tMat>(genRndST(),genRndST());
		
		for (size_t i=0; i<tMat1.N(); ++i) {
			const auto tMat2 = tMat1.cWOGet(i);
				
			CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.N()-1,tMat2.N());
			CPPUNIT_ASSERT_EQUAL(tMat1.L()-tMat1.M(),tMat2.L());
			for (size_t n=0; n<tMat1.N(); ++n) {
				if (n==i) continue;
				for (size_t m=0; m<tMat1.M(); ++m)
					CPPUNIT_ASSERT_EQUAL(tMat1(m,n),tMat2(m,n-(n>i)));
			}
		}
	}
}

template<class TT, class FT, class CT>
void test_tMat_conversion<TT,FT,CT>::test_lower() {
	// empty and single element mat
	{
		const tMat tMat1(0,0);
		const tMat tMat2(1,1);
		CPPUNIT_ASSERT(tMat1.lower().empty());
		CPPUNIT_ASSERT(tMat2.lower().empty());
	}

	// M==N
	{
		const size_t M = genRndST();
	
		const auto tMat1 = rnd<tMat>(M,M);

		const auto u = tMat1.lower();
		CPPUNIT_ASSERT_EQUAL(size_t(1),u.N());
		CPPUNIT_ASSERT_EQUAL((M-1)*M/2,u.L());

		size_t i=0;
		for (size_t n=0; n<M-1; ++n)
			for (size_t m=n+1; m<M; ++m)
				CPPUNIT_ASSERT_EQUAL(tMat1(m,n),u[i++]);
	}
	
	// M>N
	{
		auto S = genRndMEST();
		if (S.M<S.N) std::swap(S.M,S.N);
			
		const auto tMat1 = rnd<tMat>(S.M,S.N);

		const auto u = tMat1.lower();
		CPPUNIT_ASSERT_EQUAL(size_t(1),u.N());
		CPPUNIT_ASSERT_EQUAL((S.N-1)*S.N/2 + (S.M-S.N)*S.N,u.L());
			
		size_t i=0;
		for (size_t n=0; n<S.N; ++n)
			for (size_t m=n+1; m<S.M; ++m)
				CPPUNIT_ASSERT_EQUAL(tMat1(m,n),u[i++]);
	}
	
	// M<N
	{
		auto S = genRndMEST();
		if (S.M>S.N) std::swap(S.M,S.N);
			
		const auto tMat1 = rnd<tMat>(S.M,S.N);	
			
		const auto u = tMat1.lower();
		CPPUNIT_ASSERT_EQUAL(size_t(1),u.N());
		CPPUNIT_ASSERT_EQUAL((S.M-1)*S.M/2,u.L());
			
		size_t i=0;
		for (size_t n=0; n<S.M-1; ++n)
			for (size_t m=n+1; m<S.M; ++m)
				CPPUNIT_ASSERT_EQUAL(tMat1(m,n),u[i++]);
	}
}

template<class TT, class FT, class CT>
void test_tMat_conversion<TT,FT,CT>::test_upper() {
		// empty and single element mat
		{
			const tMat tMat1(0,0);
			const tMat tMat2(1,1);
			CPPUNIT_ASSERT(tMat1.upper().empty());
			CPPUNIT_ASSERT(tMat2.upper().empty());
		}

		// M==N
		{
			const size_t M = genRndST();
			const auto tMat1 = rnd<tMat>(M,M);

			const auto u = tMat1.upper();
			CPPUNIT_ASSERT_EQUAL(size_t(1),u.N());
			CPPUNIT_ASSERT_EQUAL((M-1)*M/2,u.L());

			size_t i=0;
			for (size_t n=1; n<M; ++n)
				for (size_t m=0; m<n; ++m)
					CPPUNIT_ASSERT_EQUAL(tMat1(m,n),u[i++]);
		}
		
		// M>N
		{
			auto S = genRndMEST();
			if (S.M<S.N) std::swap(S.M,S.N);
		
			const auto tMat1 = rnd<tMat>(S.M,S.N);	

			const auto u = tMat1.upper();
			CPPUNIT_ASSERT_EQUAL(size_t(1),u.N());
			CPPUNIT_ASSERT_EQUAL((S.N-1)*S.N/2,u.L());
			
			size_t i=0;
			for (size_t n=1; n<S.N; ++n)
				for (size_t m=0; m<n; ++m)
					CPPUNIT_ASSERT_EQUAL(tMat1(m,n),u[i++]);

			
		}
		
		// M<N
		{
			auto S = genRndMEST();
			if (S.M>S.N) std::swap(S.M,S.N);
			
			const auto tMat1 = rnd<tMat>(S.M,S.N);	

			const auto u = tMat1.upper();
			CPPUNIT_ASSERT_EQUAL(size_t(1),u.N());
			CPPUNIT_ASSERT_EQUAL(((S.M-1)*S.M/2 + (S.N-S.M)*S.M),u.L());
			
			size_t i=0;
			for (size_t n=1; n<S.M; ++n)
				for (size_t m=0; m<n; ++m)
					CPPUNIT_ASSERT_EQUAL(tMat1(m,n),u[i++]);
			for (size_t n=S.M; n<S.N; ++n)
				for (size_t m=0; m<S.M; ++m)
					CPPUNIT_ASSERT_EQUAL(tMat1(m,n),u[i++]);
		}
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tMat_conversion<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tMat_conversion>(
		"test_copy_fcopy_ccopy", &test_tMat_conversion<TT,FT,CT>::test_copy_fcopy_ccopy));
	suite->addTest(new CppUnit::TestCaller<test_tMat_conversion>(
		"test_get_getl", &test_tMat_conversion<TT,FT,CT>::test_get_getl));
	suite->addTest(new CppUnit::TestCaller<test_tMat_conversion>(
		"test_rGet_cGet", &test_tMat_conversion<TT,FT,CT>::test_rGet_cGet));
	suite->addTest(new CppUnit::TestCaller<test_tMat_conversion>(
		"test_rWOGet_cWOGet", &test_tMat_conversion<TT,FT,CT>::test_rWOGet_cWOGet));
	suite->addTest(new CppUnit::TestCaller<test_tMat_conversion>(
		"test_lower", &test_tMat_conversion<TT,FT,CT>::test_lower));
	suite->addTest(new CppUnit::TestCaller<test_tMat_conversion>(
		"test_upper", &test_tMat_conversion<TT,FT,CT>::test_upper));
	
	return suite;
}
