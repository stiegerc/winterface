// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_modification.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


template<class TT, class FT, class CT>
void test_tMat_modification<TT,FT,CT>::test_rShift_cShift() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// rShift by 0
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = tMat1;

		const size_t m = genRndST(0,M-1);

		tMat1.rShift(m,0);
		CPPUNIT_ASSERT(tMat2==tMat1);
	}

	// rShift random, no prealloc
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = tMat1;

		const size_t m = genRndST(0,M-1);
		const size_t s = genRndST();
		
		tMat1.rShift(m,s);
		CPPUNIT_ASSERT_EQUAL(M+s,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
		CPPUNIT_ASSERT(tMat2.rGet(0,m)==tMat1.rGet(0,m));
		CPPUNIT_ASSERT(tMat2.rGet(m,tMat2.M()-m)==tMat1.rGet(m+s,tMat2.M()-m));
	}
	
	// rShift random, no prealloc
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = tMat1;

		const size_t m = genRndST(0,M-1);
		const size_t s = genRndST();
	
		tMat1.reserve(((M+s)*N)/M+1);
		tMat1.rShift(m,s);
		CPPUNIT_ASSERT_EQUAL(M+s,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
		CPPUNIT_ASSERT(tMat2.rGet(0,m)==tMat1.rGet(0,m));
		CPPUNIT_ASSERT(tMat2.rGet(m,tMat2.M()-m)==tMat1.rGet(m+s,tMat2.M()-m));
	}
	
	// cShift by 0
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = tMat1;

		const size_t n = genRndST(0,N-1);

		tMat1.cShift(n,0);
		CPPUNIT_ASSERT(tMat2==tMat1);
	}

	// cShift random, no prealloc
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = tMat1;

		const size_t n = genRndST(0,N-1);
		const size_t s = genRndST();
		
		tMat1.cShift(n,s);
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N+s,tMat1.N());
		CPPUNIT_ASSERT(tMat2.cGet(0,n)==tMat1.cGet(0,n));
		CPPUNIT_ASSERT(tMat2.cGet(n,tMat2.N()-n)==tMat1.cGet(n+s,tMat2.N()-n));
	}
	
	// rShift random, no prealloc
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = tMat1;

		const size_t n = genRndST(0,N-1);
		const size_t s = genRndST();
	
		tMat1.reserve(N+s);
		tMat1.cShift(n,s);
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N+s,tMat1.N());
		CPPUNIT_ASSERT(tMat2.cGet(0,n)==tMat1.cGet(0,n));
		CPPUNIT_ASSERT(tMat2.cGet(n,tMat2.N()-n)==tMat1.cGet(n+s,tMat2.N()-n));
	}
}

template<class TT, class FT, class CT>
void test_tMat_modification<TT,FT,CT>::test_rInsert() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// insert matrix, no reserve
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<tMat>(genRndST(),N);

		const size_t m_ = genRndST(0,M-1);
		const auto tMat3 = mcat(mcat(tMat1.rGet(0,m_),tMat2),tMat1.rGet(m_,tMat1.M()-m_));
		tMat1.rInsert(m_,tMat2);

		CPPUNIT_ASSERT(tMat3==tMat1);
	}

	// insert matrix, reserve
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<tMat>(genRndST(),N);

		const size_t m_ = genRndST(0,M-1);
		const auto tMat3 = mcat(mcat(tMat1.rGet(0,m_),tMat2),tMat1.rGet(m_,tMat1.M()-m_));
		tMat1.reserve(((M+tMat2.M())*N)/M+1);
		tMat1.rInsert(m_,tMat2);

		CPPUNIT_ASSERT(tMat3==tMat1);
	}

	// insert row
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<tMat>(genRndST(),N);

		const size_t m_ = genRndST(0,M-1);
		const auto i = tMat2.crBegin();

		const auto tMat3 = mcat(mcat(tMat1.rGet(0,m_),*i),tMat1.rGet(m_,tMat1.M()-m_));
		tMat1.rInsert(m_,*i);

		CPPUNIT_ASSERT(tMat3==tMat1);
	}

	// insert matrix into itself
	{
		auto tMat1 = rnd<tMat>(M,N);

		const size_t m_ = genRndST(0,M-1);
		const auto tMat3 = mcat(mcat(tMat1.rGet(0,m_),tMat1),tMat1.rGet(m_,tMat1.M()-m_));
		tMat1.rInsert(m_,tMat1);

		CPPUNIT_ASSERT(tMat3==tMat1);
	}
	
	// insert row, same matrix
	{
		auto tMat1 = rnd<tMat>(M,N);

		const size_t m_ = genRndST(0,M-1);
		const auto i = tMat1.crBegin();

		const auto tMat3 = mcat(mcat(tMat1.rGet(0,m_),*i),tMat1.rGet(m_,tMat1.M()-m_));
		tMat1.rInsert(m_,*i);

		CPPUNIT_ASSERT(tMat3==tMat1);
	}
}

template<class TT, class FT, class CT>
void test_tMat_modification<TT,FT,CT>::test_cInsert() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// insert matrix, no reserve
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<tMat>(M,genRndST());


		const size_t n_ = genRndST(0,N-1);
		const auto tMat3 = ncat(ncat(tMat1.cGet(0,n_),tMat2),tMat1.cGet(n_,tMat1.N()-n_));
		
		tMat1.cInsert(n_,tMat2);

		CPPUNIT_ASSERT(tMat3==tMat1);
	}

	// insert matrix, reserve
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<tMat>(M,genRndST());

		const size_t n_ = genRndST(0,N-1);
		const auto tMat3 = ncat(ncat(tMat1.cGet(0,n_),tMat2),tMat1.cGet(n_,tMat1.N()-n_));
		tMat1.reserve(N+tMat2.N());
		tMat1.cInsert(n_,tMat2);

		CPPUNIT_ASSERT(tMat3==tMat1);
	}

	// insert row
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<tMat>(M,genRndST());

		const size_t n_ = genRndST(0,N-1);
		const auto i = tMat2.ccBegin();

		const auto tMat3 = ncat(ncat(tMat1.cGet(0,n_),*i),tMat1.cGet(n_,tMat1.N()-n_));
		tMat1.cInsert(n_,*i);

		CPPUNIT_ASSERT(tMat3==tMat1);
	}

	// insert matrix into itself
	{
		auto tMat1 = rnd<tMat>(M,N);

		const size_t n_ = genRndST(0,N-1);
		const auto tMat3 = ncat(ncat(tMat1.cGet(0,n_),tMat1),tMat1.cGet(n_,tMat1.N()-n_));
		tMat1.cInsert(n_,tMat1);

		CPPUNIT_ASSERT(tMat3==tMat1);
	}
	
	// insert row, same matrix
	{
		auto tMat1 = rnd<tMat>(M,N);

		const size_t n_ = genRndST(0,N-1);
		const auto i = tMat1.ccBegin();

		const auto tMat3 = ncat(ncat(tMat1.cGet(0,n_),*i),tMat1.cGet(n_,tMat1.N()-n_));
		tMat1.cInsert(n_,*i);

		CPPUNIT_ASSERT(tMat3==tMat1);
	}
}

template<class TT, class FT, class CT>
void test_tMat_modification<TT,FT,CT>::test_set_setl() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M,N);
	const auto m = genRndI(1,M/4,0,M/2);
	const auto n = genRndI(1,N/4,0,N/2);
		
	{
		const auto tMat2 = rnd<tMat>(m.size(),n.size());
		tMat1.set(tMat2,m,n);
		for (size_t i=0; i!=n.size(); ++i)
			for (size_t j=0; j!=m.size(); ++j)
				CPPUNIT_ASSERT_EQUAL(tMat2(j,i),tMat1(m[j],n[i]));
	}
	{
		const auto tMat2 = rnd<tMat>(M,n.size());
		tMat1.set(tMat2,{},n);
		for (size_t i=0; i!=n.size(); ++i)
			for (size_t j=0; j!=tMat1.M(); ++j)
				CPPUNIT_ASSERT_EQUAL(tMat2(j,i),tMat1(j,n[i]));
	}
	{
		const auto tMat2 = rnd<tMat>(m.size(),N);
		tMat1.set(tMat2,m,{});
		for (size_t i=0; i!=tMat1.N(); ++i)
			for (size_t j=0; j!=m.size(); ++j)
				CPPUNIT_ASSERT_EQUAL(tMat2(j,i),tMat1(m[j],i));
	}
	{
		const auto tMat2 = rnd<tMat>(M,N);
		tMat1.set(tMat2,{},{});
		CPPUNIT_ASSERT(tMat1==tMat2);
	}
	{
		const auto tMat2 = rnd<tMat>(M/2,N/2);
		const size_t ms = genRndST(0,M/2-1);
		const size_t ns = genRndST(0,N/2-1);
		
		tMat1.set(tMat2,ms,ns);
		for (size_t n=0; n!=tMat2.N(); ++n)
			for (size_t m=0; m!=tMat2.M(); ++m)
				CPPUNIT_ASSERT_EQUAL(tMat2(m,n),tMat1(m+ms,n+ns));
	}		
	{
		const auto ml = nsum(tMat1).lt(mean(nsum(tMat1)));
		const auto nl = msum(tMat1).lt(mean(msum(tMat1)));
	
		const auto tMat3 = rnd<tMat>(size_t(sum(ml)),size_t(sum(nl)));
	
		std::vector<size_t> mlck; mlck.reserve(size_t(sum(ml)));
		for (size_t i=0; i!=ml.L(); ++i)
			if (bool(ml[i])) mlck.push_back(i);
		
		std::vector<size_t> nlck; nlck.reserve(size_t(sum(nl)));
		for (size_t i=0; i!=nl.L(); ++i)
			if (bool(nl[i])) nlck.push_back(i);
	
		auto ck = tMat1.copy().setl(tMat3,ml,nl)==tMat1.copy().set(tMat3,mlck,nlck);
		CPPUNIT_ASSERT(ck);
	}
}

template<class TT, class FT, class CT>
void test_tMat_modification<TT,FT,CT>::test_rRm() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	const auto tMat1 = rnd<tMat>(M,N);
	
	// one row
	{
		for (size_t i=0; i!=tMat1.M(); ++i) {
			auto tMat2=tMat1;
			tMat2.rRm(i);
			
			CPPUNIT_ASSERT_EQUAL(tMat1.M()-1,tMat2.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat2.N());
			CPPUNIT_ASSERT_EQUAL(tMat1.L()-tMat1.N(),tMat2.L());
			for (size_t m=0; m!=tMat1.M(); ++m) {
				if (m==i) continue;
				for (size_t n=0; n<tMat1.N(); ++n)
					CPPUNIT_ASSERT_EQUAL(tMat1(m,n),tMat2(m-(m>i),n));
			}
			CPPUNIT_ASSERT_EQUAL(tMat1.lcap(),tMat2.lcap());
		}
	}

	// multiple rows
	{
		// empty vector
		{
			auto tMat2=tMat1;
			tMat2.rRm(std::vector<size_t>());
			CPPUNIT_ASSERT(tMat1==tMat2);
			CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
		}
			
		// remove all rows
		{
			auto tMat2=tMat1;
			std::vector<size_t> m(tMat1.M());
			for (size_t i=0; i!=m.size(); ++i)
				m[i]=i;
			
			tMat2.rRm(m);
			CPPUNIT_ASSERT(tMat2.empty());
			CPPUNIT_ASSERT(!tMat2.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat2.N());
			CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
		}
			
		// including edges
		{
			std::vector<size_t> m = {0, M/4, M/2, 2*M/3, M-1};
			auto tMat2=tMat1;

			tMat2.rRm(m);
			
			CPPUNIT_ASSERT_EQUAL(tMat1.M()-m.size(),tMat2.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat2.N());
			CPPUNIT_ASSERT_EQUAL(tMat1.L()-m.size()*tMat1.N(),tMat2.L());
		
			for (auto i: m)
				CPPUNIT_ASSERT_EQUAL(tMat2.crEnd(),
				std::find(tMat2.crBegin(),tMat2.crEnd(),tMat1.rAt(i)));
				
			std::vector<size_t> pos(tMat2.M());
			for (size_t i=0; i<tMat2.M(); ++i)
				pos[i] = std::find(tMat1.crBegin(),tMat1.crEnd(),tMat2.rAt(i))-tMat1.crBegin();
			
			for (size_t i=1; i!=pos.size(); ++i)
				CPPUNIT_ASSERT(pos[i]>pos[i-1]);
				
			CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
		}
			
		// excluding edges
		{
			std::vector<size_t> m = {M/5, M/3, M/2, 3*M/4};
			auto tMat2=tMat1;
			tMat2.rRm(m);
				
			CPPUNIT_ASSERT_EQUAL(tMat1.M()-m.size(),tMat2.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat2.N());
			CPPUNIT_ASSERT_EQUAL(tMat1.L()-m.size()*tMat1.N(),tMat2.L());
				
			for (auto i: m)
				CPPUNIT_ASSERT_EQUAL(tMat2.crEnd(),
				std::find(tMat2.crBegin(),tMat2.crEnd(),tMat1.rAt(i)));
				
			std::vector<size_t> pos(tMat2.M());
			for (size_t i=0; i<tMat2.M(); ++i)
				pos[i] = std::find(tMat1.crBegin(),tMat1.crEnd(),tMat2.rAt(i))-tMat1.crBegin();
			
			for (size_t i=1; i!=pos.size(); ++i)
				CPPUNIT_ASSERT(pos[i]>pos[i-1]);
				
			CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
		}
	}
}

template<class TT, class FT, class CT>
void test_tMat_modification<TT,FT,CT>::test_cRm() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	const auto tMat1 = rnd<tMat>(M,N);
	
	// one column
	{	
		for (size_t i=0; i<tMat1.N(); ++i) {
			auto tMat2=tMat1;
			tMat2.cRm(i);
				
			CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.N()-1,tMat2.N());
			CPPUNIT_ASSERT_EQUAL(tMat1.L()-tMat1.M(),tMat2.L());
			for (size_t n=0; n<tMat1.N(); ++n) {
				if (n==i) continue;
				for (size_t m=0; m!=tMat1.M(); ++m)
					CPPUNIT_ASSERT_EQUAL(tMat1(m,n),tMat2(m,n-(n>i)));
			}
			CPPUNIT_ASSERT_EQUAL(tMat1.lcap(),tMat2.lcap());
		}
	}
	
	// multiple columns
	{
		// empty vector
		{
			auto tMat2=tMat1;
			tMat2.cRm(std::vector<size_t>());
			CPPUNIT_ASSERT(tMat2==tMat1);
			CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
		}
			
		// remove all
		{
			auto tMat2=tMat1;
			std::vector<size_t> n(tMat1.N());
			std::iota(n.begin(),n.end(),0);
				
			tMat2.cRm(n);
			CPPUNIT_ASSERT(tMat2.empty());
			CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.M());
			CPPUNIT_ASSERT(!tMat2.N());
			CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
		}
			
		// including edges
		{
			std::vector<size_t> n = {0, N/4, N/2, N-1};
			auto tMat2=tMat1;
			tMat2.cRm(n);
				
			CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.N()-n.size(),tMat2.N());
			CPPUNIT_ASSERT_EQUAL(tMat1.L()-n.size()*tMat1.M(),tMat2.L());
				
			for (auto i: n)
				CPPUNIT_ASSERT_EQUAL(tMat2.ccEnd(),
				std::find(tMat2.ccBegin(),tMat2.ccEnd(),tMat1.cAt(i)));
				
			std::vector<size_t> pos(tMat2.N());
			for (size_t i=0; i!=tMat2.N(); ++i)
				pos[i] = std::find(tMat1.ccBegin(),tMat1.ccEnd(),tMat2.cAt(i))-tMat1.ccBegin();
					
			for (size_t i=1; i!=pos.size(); ++i)
				CPPUNIT_ASSERT(pos[i]>pos[i-1]);
				
			CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
		}
			
		// excluding edges
		{
			std::vector<size_t> n = {N/7, N/4, N/2};
			auto tMat2=tMat1;
			tMat2.cRm(n);
				
			CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.N()-n.size(),tMat2.N());
			CPPUNIT_ASSERT_EQUAL(tMat1.L()-n.size()*tMat1.M(),tMat2.L());
				
			for (auto i: n)
				CPPUNIT_ASSERT_EQUAL(tMat2.ccEnd(),
				std::find(tMat2.ccBegin(),tMat2.ccEnd(),tMat1.cAt(i)));
				
			std::vector<size_t> pos(tMat2.N());
			for (size_t i=0; i!=tMat2.N(); ++i)
				pos[i] = std::find(tMat1.ccBegin(),tMat1.ccEnd(),tMat2.cAt(i))-tMat1.ccBegin();
					
			for (size_t i=1; i!=pos.size(); ++i)
				CPPUNIT_ASSERT(pos[i]>pos[i-1]);				
				
			CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
		}
	}
}

template<class TT, class FT, class CT>
void test_tMat_modification<TT,FT,CT>::test_dRm() {
	const size_t M = genRndST();
	const size_t N = genRndST();
	const auto tMat1 = rnd<tMat>(M,N);
	
	// M!=N
	{
		const auto S = genRndMEST();
		auto tMat1 = rnd<tMat>(S.M,S.N);

		auto tMat2=tMat1;
		tMat1.dRm();
			
		if (S.N>S.M) { // add low
			CPPUNIT_ASSERT_EQUAL(tMat2.N()-1,tMat1.N());
			CPPUNIT_ASSERT_EQUAL(tMat2.M(),tMat1.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.M()*tMat1.N(),tMat1.L());
			CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());

			// lower part
			for (size_t n=0; n!=tMat1.M()-1; ++n)
				for (size_t m=n+1; m!=tMat1.M(); ++m)
					CPPUNIT_ASSERT_EQUAL(tMat2(m,n),tMat1(m,n));

			// upper part
			for (size_t n=0; n!=tMat1.M()-1; ++n)
				for (size_t m=0; m<=n; ++m)
					CPPUNIT_ASSERT_EQUAL(tMat2(m,n+1),tMat1(m,n));
				
			// right block
			for (size_t n=tMat1.M()-1; n!=tMat1.N(); ++n)
				for (size_t m=0; m!=tMat1.M(); ++m)
					CPPUNIT_ASSERT_EQUAL(tMat2(m,n+1),tMat1(m,n));
					
			} else { // add left
				CPPUNIT_ASSERT_EQUAL(tMat2.N(),tMat1.N());
				CPPUNIT_ASSERT_EQUAL(tMat2.M()-1,tMat1.M());
				CPPUNIT_ASSERT_EQUAL(tMat1.M()*tMat1.N(),tMat1.L());
				CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
			
				// lower part
				for (size_t n=0; n!=tMat1.N()-1; ++n)
					for (size_t m=n; m!=tMat1.N()-1; ++m)
						CPPUNIT_ASSERT_EQUAL(tMat2(m+1,n),tMat1(m,n));

				// upper part
				for (size_t n=1; n!=tMat1.N(); ++n)
					for (size_t m=0; m!=n; ++m)
						CPPUNIT_ASSERT_EQUAL(tMat2(m,n),tMat1(m,n));

				// low block
				for (size_t n=0; n!=tMat1.N(); ++n)
					for (size_t m=tMat1.N()-1; m!=tMat1.M(); ++m)
						CPPUNIT_ASSERT_EQUAL(tMat2(m+1,n),tMat1(m,n));
			}
	}
		
	// M==N;
	{
		const size_t M = genRndST();
		auto tMat1 = rnd<tMat>(M,M);

		// low
		{
			auto tMat2 = tMat1;
			tMat2.dRm(true);
				
			CPPUNIT_ASSERT_EQUAL(tMat2.N(),tMat1.N()-1);
			CPPUNIT_ASSERT_EQUAL(tMat2.M(),tMat1.M());
			CPPUNIT_ASSERT_EQUAL(tMat2.L(),tMat2.M()*tMat2.N());
			CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());

			// lower part
			for (size_t n=0; n!=tMat1.M()-1; ++n)
				for (size_t m=n+1; m!=tMat1.M(); ++m)
					CPPUNIT_ASSERT_EQUAL(tMat1(m,n),tMat2(m,n));

			// upper part
			for (size_t n=0; n!=tMat1.M()-1; ++n)
				for (size_t m=0; m<=n; ++m)
					CPPUNIT_ASSERT_EQUAL(tMat1(m,n+1),tMat2(m,n));
		}
	
		// left
		{
			auto tMat2 = tMat1;
			tMat2.dRm(false);

			CPPUNIT_ASSERT_EQUAL(tMat2.N(),tMat1.N());
			CPPUNIT_ASSERT_EQUAL(tMat2.M(),tMat1.M()-1);
			CPPUNIT_ASSERT_EQUAL(tMat2.M()*tMat2.N(),tMat2.L());
			CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
			
			// lower part
			for (size_t n=0; n!=tMat1.N()-1; ++n)
				for (size_t m=n; m!=tMat1.N()-1; ++m)
					CPPUNIT_ASSERT_EQUAL(tMat1(m+1,n),tMat2(m,n));

			// upper part
			for (size_t n=1; n!=tMat1.N(); ++n)
				for (size_t m=0; m!=n; ++m)
					CPPUNIT_ASSERT_EQUAL(tMat2(m,n),tMat1(m,n));
		}
	}

	// row
	{
		const size_t N = genRndST();
		auto tMat1 = rnd<tMat>(1,N);

		auto tMat2=tMat1;
		tMat1.dRm();

		CPPUNIT_ASSERT(tMat1.row());
		CPPUNIT_ASSERT_EQUAL(tMat2.N()-1,tMat1.N());
		CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());

		for (size_t i=0; i<tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat2[i+1],tMat1[i]);
	}
		
	// col
	{
		const size_t M = genRndST();
		auto tMat1 = rnd<tMat>(M,1);

		auto tMat2=tMat1;
		tMat1.dRm();

		CPPUNIT_ASSERT(tMat1.col());
		CPPUNIT_ASSERT_EQUAL(tMat2.M()-1,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat2[i+1],tMat1[i]);
	}

	// 1x1
	{
		tMat tMat1(1);
		tMat1.dRm();
				
		CPPUNIT_ASSERT(tMat1.empty());
		CPPUNIT_ASSERT(!tMat1.M());
		CPPUNIT_ASSERT(!tMat1.N());
		CPPUNIT_ASSERT(!tMat1.L());
		CPPUNIT_ASSERT_EQUAL(size_t(1),tMat1.lcap());
	}
}


template<class TT, class FT, class CT>
void test_tMat_modification<TT,FT,CT>::test_inv() {
	const size_t M = genRndST();

	const auto tMat1 = rnd_b<tMat>(M,M);
	auto tMat2 = eye<tMat>(M);
	auto tMat3 = tMat1;
	
	tMat3.inv();
	const auto tMat4 = tMat1.prod(tMat3);

	for (size_t i=0; i!=tMat1.L(); ++i)
		CPPUNIT_ASSERT_DELTA(tMat2[i],tMat4[i],delta);
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tMat_modification<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tMat_modification>(
		"test_rShift_cShift", &test_tMat_modification<TT,FT,CT>::test_rShift_cShift));
	suite->addTest(new CppUnit::TestCaller<test_tMat_modification>(
		"test_rInsert", &test_tMat_modification<TT,FT,CT>::test_rInsert));
	suite->addTest(new CppUnit::TestCaller<test_tMat_modification>(
		"test_cInsert", &test_tMat_modification<TT,FT,CT>::test_cInsert));
	suite->addTest(new CppUnit::TestCaller<test_tMat_modification>(
		"test_set_setl", &test_tMat_modification<TT,FT,CT>::test_set_setl));
	suite->addTest(new CppUnit::TestCaller<test_tMat_modification>(
		"test_rRm", &test_tMat_modification<TT,FT,CT>::test_rRm));
	suite->addTest(new CppUnit::TestCaller<test_tMat_modification>(
		"test_cRm", &test_tMat_modification<TT,FT,CT>::test_cRm));
	suite->addTest(new CppUnit::TestCaller<test_tMat_modification>(
		"test_dRm", &test_tMat_modification<TT,FT,CT>::test_dRm));
	suite->addTest(new CppUnit::TestCaller<test_tMat_modification>(
		"test_inv", &test_tMat_modification<TT,FT,CT>::test_inv));
	
	return suite;
}
