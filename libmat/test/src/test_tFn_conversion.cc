// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tFn_conversion.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_msum() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(0,N);
		const auto s = msum(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(0),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto s = msum(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());

		auto c = tMat1.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c)
			CPPUNIT_ASSERT_DELTA(sum(*c),s[n],delta);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_nsum() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,0);
		const auto s = nsum(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),s.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto s = nsum(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.N());

		auto r = tMat1.crBegin();
		for (size_t m=0; m!=M; ++m, ++r)
			CPPUNIT_ASSERT_DELTA(sum(*r),s[m],delta);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_mprod() {
	const size_t M = genRndST(2,4);
	const size_t N = genRndST(2,4);

	{
		const auto tMat1 = rnd<tMat>(0,N,1.0,1.01);
		const auto s = mprod(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(0),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N,1.0,1.01);
		const auto s = mprod(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());

		auto c = tMat1.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c)
			CPPUNIT_ASSERT_DELTA(prod(*c),s[n],delta);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_nprod() {
	const size_t M = genRndST(2,4);
	const size_t N = genRndST(2,4);

	{
		const auto tMat1 = rnd<tMat>(M,0,1.0,1.01);
		const auto s = nprod(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),s.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N,1.0,1.01);
		const auto s = nprod(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.N());

		auto r = tMat1.crBegin();
		for (size_t m=0; m!=M; ++m, ++r)
			CPPUNIT_ASSERT_DELTA(prod(*r),s[m],delta);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_mmin() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(0,N);
		const auto mm = mmin(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(0),mm.mat.M());
		CPPUNIT_ASSERT_EQUAL(N,mm.mat.N());
		CPPUNIT_ASSERT(mm.pos.empty());
	}
	{
		const auto tMat1 = FT(1000.0)*rnd<tMat>(M,N);
		const auto mm = mmin(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(1),mm.mat.M());
		CPPUNIT_ASSERT_EQUAL(N,mm.mat.N());
		CPPUNIT_ASSERT_EQUAL(size_t(N),mm.pos.size());

		auto c = tMat1.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c) {
			const auto ref = std::min_element(c->cbegin(),c->cend(),
					[](const TT& a, const TT& b){return ops::lt_s(a,b);});
			CPPUNIT_ASSERT_EQUAL(*ref,mm.mat[n]);
			CPPUNIT_ASSERT_EQUAL(size_t(ref-c->cbegin()),mm.pos[n]);
		}
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_nmin() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,0);
		const auto mm = nmin(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,mm.mat.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),mm.mat.N());
		CPPUNIT_ASSERT(mm.pos.empty());
	}
	{
		const auto tMat1 = FT(1000.0)*rnd<tMat>(M,N);
		const auto mm = nmin(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,mm.mat.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),mm.mat.N());
		CPPUNIT_ASSERT_EQUAL(size_t(M),mm.pos.size());

		auto r = tMat1.crBegin();
		for (size_t m=0; m!=M; ++m, ++r) {
			const auto ref = std::min_element(r->cbegin(),r->cend(),
					[](const TT& a, const TT& b){return ops::lt_s(a,b);});
			CPPUNIT_ASSERT_EQUAL(*ref,mm.mat[m]);
			CPPUNIT_ASSERT_EQUAL(size_t(ref-r->cbegin()),mm.pos[m]);
		}
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_mmax() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(0,N);
		const auto mm = mmax(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(0),mm.mat.M());
		CPPUNIT_ASSERT_EQUAL(N,mm.mat.N());
		CPPUNIT_ASSERT(mm.pos.empty());
	}
	{
		const auto tMat1 = FT(1000.0)*rnd<tMat>(M,N);
		const auto mm = mmax(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(1),mm.mat.M());
		CPPUNIT_ASSERT_EQUAL(N,mm.mat.N());
		CPPUNIT_ASSERT_EQUAL(size_t(N),mm.pos.size());

		auto c = tMat1.cBegin();
		for (size_t n=0; n!=N; ++n, ++c) {
			const auto ref = std::max_element(c->cbegin(),c->cend(),
					[](const TT& a, const TT& b){return ops::lt_s(a,b);});
			CPPUNIT_ASSERT_EQUAL(*ref,mm.mat[n]);
			CPPUNIT_ASSERT_EQUAL(size_t(ref-c->cbegin()),mm.pos[n]);
		}
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_nmax() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,0);
		const auto mm = nmax(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,mm.mat.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),mm.mat.N());
		CPPUNIT_ASSERT(mm.pos.empty());
	}
	{
		const auto tMat1 = FT(1000.0)*rnd<tMat>(M,N);
		const auto mm = nmax(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,mm.mat.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),mm.mat.N());
		CPPUNIT_ASSERT_EQUAL(size_t(M),mm.pos.size());

		auto r = tMat1.crBegin();
		for (size_t m=0; m!=M; ++m, ++r) {
			const auto ref = std::max_element(r->cbegin(),r->cend(),
					[](const TT& a, const TT& b){return ops::lt_s(a,b);});
			CPPUNIT_ASSERT_EQUAL(*ref,mm.mat[m]);
			CPPUNIT_ASSERT_EQUAL(size_t(ref-r->cbegin()),mm.pos[m]);
		}
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_mmean() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(0,N);
		const auto mm = mmean(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(0),mm.M());
		CPPUNIT_ASSERT_EQUAL(N,mm.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto mm = mmean(tMat1);
		const auto ref = msum(tMat1)/FT(M);
		CPPUNIT_ASSERT(size(mm)==size(ref));

		for (size_t n=0; n!=N; ++n)
			CPPUNIT_ASSERT_DELTA(ref[n],mm[n],delta);	
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_nmean() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,0);
		const auto mm = nmean(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,mm.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),mm.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto mm = nmean(tMat1);
		const auto ref = nsum(tMat1)/FT(N);
		CPPUNIT_ASSERT(size(mm)==size(ref));

		for (size_t m=0; m!=M; ++m)
			CPPUNIT_ASSERT_DELTA(ref[m],mm[m],delta);	
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_mnormsq() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(0,N);
		const auto s = mnormsq(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(0),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto s = mnormsq(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());

		auto c = tMat1.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c)
			CPPUNIT_ASSERT_DELTA(normsq(*c),s[n],delta);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_nnormsq() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,0);
		const auto s = nnormsq(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),s.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto s = nnormsq(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.N());

		auto r = tMat1.crBegin();
		for (size_t m=0; m!=M; ++m, ++r)
			CPPUNIT_ASSERT_DELTA(normsq(*r),s[m],delta);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_mnorm() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(0,N);
		const auto s = mnorm(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(0),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto s = mnorm(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());

		auto c = tMat1.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c)
			CPPUNIT_ASSERT_DELTA(norm(*c),s[n],delta);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_nnorm() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,0);
		const auto s = nnorm(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),s.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto s = nnorm(tMat1);
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.N());

		auto r = tMat1.crBegin();
		for (size_t m=0; m!=M; ++m, ++r)
			CPPUNIT_ASSERT_DELTA(norm(*r),s[m],delta);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_mdot() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(0,N);
		const auto ftMat = rnd<tMat>(0,N);
		const auto s = mdot(tMat1,ftMat);
		CPPUNIT_ASSERT(tMat1.cpx() == s.cpx());
		CPPUNIT_ASSERT_EQUAL(size_t(0),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ftMat = rnd<fMat>(M,N);
		const auto s = mdot(tMat1,ftMat);
		CPPUNIT_ASSERT(tMat1.cpx() == s.cpx());
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());

		auto c = tMat1.ccBegin();
		auto fc = ftMat.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c, ++fc)
			CPPUNIT_ASSERT_DELTA(dot(*c,*fc),s[n],delta);
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ctMat = rnd<cMat>(M,N);
		const auto s = mdot(tMat1,ctMat);
		CPPUNIT_ASSERT(s.cpx());
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());

		auto c = tMat1.ccBegin();
		auto cc = ctMat.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c, ++cc)
			CPPUNIT_ASSERT_DELTA(dot(*c,*cc),s[n],delta);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_ndot() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,0);
		const auto ftMat = rnd<tMat>(M,0);
		const auto s = ndot(tMat1,ftMat);
		CPPUNIT_ASSERT(tMat1.cpx() == s.cpx());
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),s.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ftMat = rnd<fMat>(M,N);
		const auto s = ndot(tMat1,ftMat);
		CPPUNIT_ASSERT(tMat1.cpx() == s.cpx());
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.N());

		auto r = tMat1.crBegin();
		auto fr = ftMat.crBegin();
		for (size_t m=0; m!=M; ++m, ++r, ++fr)
			CPPUNIT_ASSERT_DELTA(dot(*r,*fr),s[m],delta);
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ctMat = rnd<cMat>(M,N);
		const auto s = ndot(tMat1,ctMat);
		CPPUNIT_ASSERT(s.cpx());
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.N());

		auto r = tMat1.crBegin();
		auto cr = ctMat.crBegin();
		for (size_t m=0; m!=M; ++m, ++r, ++cr)
			CPPUNIT_ASSERT_DELTA(dot(*r,*cr),s[m],delta);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_mdotu() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(0,N);
		const auto ftMat = rnd<tMat>(0,N);
		const auto s = mdotu(tMat1,ftMat);
		CPPUNIT_ASSERT(tMat1.cpx() == s.cpx());
		CPPUNIT_ASSERT_EQUAL(size_t(0),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ftMat = rnd<fMat>(M,N);
		const auto s = mdotu(tMat1,ftMat);
		CPPUNIT_ASSERT(tMat1.cpx() == s.cpx());
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());

		auto c = tMat1.ccBegin();
		auto fc = ftMat.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c, ++fc)
			CPPUNIT_ASSERT_DELTA(dotu(*c,*fc),s[n],delta);
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ctMat = rnd<cMat>(M,N);
		const auto s = mdotu(tMat1,ctMat);
		CPPUNIT_ASSERT(s.cpx());
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.M());
		CPPUNIT_ASSERT_EQUAL(N,s.N());

		auto c = tMat1.ccBegin();
		auto cc = ctMat.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c, ++cc)
			CPPUNIT_ASSERT_DELTA(dotu(*c,*cc),s[n],delta);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_ndotu() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,0);
		const auto ftMat = rnd<tMat>(M,0);
		const auto s = ndotu(tMat1,ftMat);
		CPPUNIT_ASSERT(tMat1.cpx() == s.cpx());
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(0),s.N());
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ftMat = rnd<fMat>(M,N);
		const auto s = ndotu(tMat1,ftMat);
		CPPUNIT_ASSERT(tMat1.cpx() == s.cpx());
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.N());

		auto r = tMat1.crBegin();
		auto fr = ftMat.crBegin();
		for (size_t m=0; m!=M; ++m, ++r, ++fr)
			CPPUNIT_ASSERT_DELTA(dotu(*r,*fr),s[m],delta);
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ctMat = rnd<cMat>(M,N);
		const auto s = ndotu(tMat1,ctMat);
		CPPUNIT_ASSERT(s.cpx());
		CPPUNIT_ASSERT_EQUAL(M,s.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),s.N());

		auto r = tMat1.crBegin();
		auto cr = ctMat.crBegin();
		for (size_t m=0; m!=M; ++m, ++r, ++cr)
			CPPUNIT_ASSERT_DELTA(dotu(*r,*cr),s[m],delta);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_inv() {
	const size_t M = genRndST(3,3);

	auto tMat1 = rnd<tMat>(M,M);
	while (std::abs(det(tMat1))<FT(.5))
		tMat1 = rnd<tMat>(M,M);

	const auto itMat1 = inv(tMat1);
	const auto act = tMat1.prod(itMat1);
	const auto ref = eye<tMat>(M,M);

	for (size_t i=0; i!=M*M; ++i)
		CPPUNIT_ASSERT_DELTA(ref[i],act[i],delta);
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_diag() {
	const size_t M = genRndST();
	const size_t N = genRndST();
		
	// check creating matrix from row
	{
		const auto tMat1 = rnd<tMat>(1,N);
			
		// default os (=0)
		{
			const auto tMat2 = diag(tMat1);
			CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat2.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat2.N());
			for (size_t i=0; i<tMat1.L(); ++i)
				CPPUNIT_ASSERT_EQUAL(tMat1[i],tMat2(i,i));
		}
		
		// check a small range of os along m, from matrix
		{
			for (size_t os=1; os<4; ++os) {
				const auto tMat2 = diag(tMat1,os,true);
				
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.M());
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.N());
			
				auto i = tMat2.cdbegin(os,true);
				auto ie = tMat2.cdend(os,true);
				CPPUNIT_ASSERT_EQUAL(size_t(ie-i),tMat1.L());
			
				auto j1 = tMat1.cbegin();
				for (auto j2=tMat2.cbegin(),je=tMat2.cend(); j2!=je; ++j2)
					if (j2==i && i!=ie) CPPUNIT_ASSERT_EQUAL(*j1++,*i++);
					else CPPUNIT_ASSERT_EQUAL(TT(0.0),*j2);
			}
		}
		
		// check a small range of os along n, from matrix
		{
			for (size_t os=1; os<4; ++os) {
				const auto tMat2 = diag(tMat1,os,false);
				
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.M());
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.N());
			
				auto i = tMat2.cdbegin(os,false);
				auto ie = tMat2.cdend(os,false);
				CPPUNIT_ASSERT_EQUAL(size_t(ie-i),tMat1.L());
				
				auto j1 = tMat1.cbegin();
				for (auto j2=tMat2.cbegin(),je=tMat2.cend(); j2!=je; ++j2)
					if (j2==i && i!=ie) CPPUNIT_ASSERT_EQUAL(*j1++,*i++);
					else CPPUNIT_ASSERT_EQUAL(TT(0.0),*j2);
			}
		}

		// check a small range of os along m, from matrix
		{
			for (size_t os=1; os<4; ++os) {
				const auto tMat2 = diag(*tMat1.crBegin(),os,true);
				
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.M());
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.N());
			
				auto i = tMat2.cdbegin(os,true);
				auto ie = tMat2.cdend(os,true);
				CPPUNIT_ASSERT_EQUAL(size_t(ie-i),tMat1.L());
				
				auto j1 = tMat1.cbegin();
				for (auto j2=tMat2.cbegin(),je=tMat2.cend(); j2!=je; ++j2)
					if (j2==i && i!=ie) CPPUNIT_ASSERT_EQUAL(*j1++,*i++);
					else CPPUNIT_ASSERT_EQUAL(TT(0.0),*j2);
			}
		}
		
		// check a small range of os along n, from matrix
		{
			for (size_t os=1; os<4; ++os) {
				const auto tMat2 = diag(*tMat1.crBegin(),os,false);
				
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.M());
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.N());
			
				auto i = tMat2.cdbegin(os,false);
				auto ie = tMat2.cdend(os,false);
				CPPUNIT_ASSERT_EQUAL(size_t(ie-i),tMat1.L());
				
				auto j1 = tMat1.cbegin();
				for (auto j2=tMat2.cbegin(),je=tMat2.cend(); j2!=je; ++j2)
					if (j2==i && i!=ie) CPPUNIT_ASSERT_EQUAL(*j1++,*i++);
					else CPPUNIT_ASSERT_EQUAL(TT(0.0),*j2);
			}
		}
	}

	// check creating matrix from col
	{
		const auto tMat1 = rnd<tMat>(M,1);
			
		// default os (=0)
		{
			const auto tMat2 = diag(tMat1);
			CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat2.M());
			CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat2.N());
			for (size_t i=0; i<tMat1.L(); ++i)
				CPPUNIT_ASSERT_EQUAL(tMat1[i],tMat2(i,i));
		}
		
		// check a small range of os along m, from matrix
		{
			for (size_t os=1; os<4; ++os) {
				const auto tMat2 = diag(tMat1,os,true);
				
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.M());
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.N());
			
				auto i = tMat2.cdbegin(os,true);
				auto ie = tMat2.cdend(os,true);
				CPPUNIT_ASSERT_EQUAL(size_t(ie-i),tMat1.L());
				
				auto j1 = tMat1.cbegin();
				for (auto j2=tMat2.cbegin(),je=tMat2.cend(); j2!=je; ++j2)
					if (j2==i && i!=ie) CPPUNIT_ASSERT_EQUAL(*j1++,*i++);
					else CPPUNIT_ASSERT_EQUAL(TT(0.0),*j2);
			}
		}
		
		// check a small range of os along n, from matrix
		{
			for (size_t os=1; os<4; ++os) {
				const auto tMat2 = diag(tMat1,os,false);
				
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.M());
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.N());
			
				auto i = tMat2.cdbegin(os,false);
				auto ie = tMat2.cdend(os,false);
				CPPUNIT_ASSERT_EQUAL(size_t(ie-i),tMat1.L());
				
				auto j1 = tMat1.cbegin();
				for (auto j2=tMat2.cbegin(),je=tMat2.cend(); j2!=je; ++j2)
					if (j2==i && i!=ie) CPPUNIT_ASSERT_EQUAL(*j1++,*i++);
					else CPPUNIT_ASSERT_EQUAL(TT(0.0),*j2);
			}
		}

		// check a small range of os along m, from matrix
		{
			for (size_t os=1; os<4; ++os) {
				const auto tMat2 = diag(*tMat1.ccBegin(),os,true);
				
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.M());
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.N());
			
				auto i = tMat2.cdbegin(os,true);
				auto ie = tMat2.cdend(os,true);
				CPPUNIT_ASSERT_EQUAL(size_t(ie-i),tMat1.L());
				
				auto j1 = tMat1.cbegin();
				for (auto j2=tMat2.cbegin(),je=tMat2.cend(); j2!=je; ++j2)
					if (j2==i && i!=ie) CPPUNIT_ASSERT_EQUAL(*j1++,*i++);
					else CPPUNIT_ASSERT_EQUAL(TT(0.0),*j2);
			}
		}
		
		// check a small range of os along n, from matrix
		{
			for (size_t os=1; os<4; ++os) {
				const auto tMat2 = diag(*tMat1.ccBegin(),os,false);
				
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.M());
				CPPUNIT_ASSERT_EQUAL(tMat1.L()+os,tMat2.N());
			
				auto i = tMat2.cdbegin(os,false);
				auto ie = tMat2.cdend(os,false);
				CPPUNIT_ASSERT_EQUAL(size_t(ie-i),tMat1.L());
				
				auto j1 = tMat1.cbegin();
				for (auto j2=tMat2.cbegin(),je=tMat2.cend(); j2!=je; ++j2)
					if (j2==i && i!=ie) CPPUNIT_ASSERT_EQUAL(*j1++,*i++);
					else CPPUNIT_ASSERT_EQUAL(TT(0.0),*j2);
			}
		}
	}

	// check creating col of diagonal from matrix
	{
		const auto tMat1 = rnd<tMat>(M,N);
		
		// default os(=0)
		{
			const auto tMat2 = diag(tMat1);
			
			const size_t S = std::min(M,N);
			CPPUNIT_ASSERT_EQUAL(S,tMat2.L());
		
			for (size_t m=0; m!=S; ++m)
				CPPUNIT_ASSERT_EQUAL(tMat1(m,m),tMat2[m]);
		}

		// check s small range of os along m
		{
			for (size_t os=1; os<4; ++os) {
				const auto tMat2 = diag(tMat1,os,true);

				auto i = tMat1.cdbegin(os,true);
				auto ie = tMat1.cdend(os,true);

				CPPUNIT_ASSERT_EQUAL(size_t(ie-i),tMat2.L());
				for (auto j: tMat2)
					CPPUNIT_ASSERT_EQUAL(*i++,j);
			}
		}
		
		// check s small range of os along m
		{
			for (size_t os=1; os<4; ++os) {
				const auto tMat2 = diag(tMat1,os,false);

				auto i = tMat1.cdbegin(os,false);
				auto ie = tMat1.cdend(os,false);

				CPPUNIT_ASSERT_EQUAL(size_t(ie-i),tMat2.L());
				for (auto j: tMat2)
					CPPUNIT_ASSERT_EQUAL(*i++,j);
			}
		}
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_R_C() {
	const size_t M = genRndST();
	const size_t N = genRndST();
	const auto tMat1 = rnd<tMat>(M,N);

	{
		const auto tMat2 = R(tMat1);
		CPPUNIT_ASSERT_EQUAL(size_t(1),tMat2.M());
		CPPUNIT_ASSERT_EQUAL(M*N,tMat2.N());
		CPPUNIT_ASSERT_EQUAL(tMat1.lcap(),tMat2.lcap());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i],tMat2[i]);
	}
	{
		const auto tMat2 = C(tMat1);
		CPPUNIT_ASSERT_EQUAL(M*N,tMat2.M());
		CPPUNIT_ASSERT_EQUAL(size_t(1),tMat2.N());
		CPPUNIT_ASSERT_EQUAL(tMat1.lcap(),tMat2.lcap());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i],tMat2[i]);
	}
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_T() {
	const size_t M = genRndST();
	const auto tMat1 = rnd<tMat>(M,1);
	const auto tMat2 = T(tMat1);
	
			
	CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.N());
	CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat2.M());
			
	for (size_t n=0; n<tMat1.N(); ++n)
		for (size_t m=0; m<tMat1.M(); ++m)
			CPPUNIT_ASSERT_EQUAL(CT(tMat1(m,n)),CT(std::conj(tMat2(n,m))));
	

	const auto tMat3 = rnd<tMat>(M,M);
	const auto tMat4 = T(tMat3);
	for (size_t n=0; n<tMat1.N(); ++n)
		for (size_t m=0; m<tMat1.M(); ++m)
			CPPUNIT_ASSERT_EQUAL(CT(tMat1(m,n)),CT(std::conj(tMat2(n,m))));
		
	const auto S = genRndMEST();
	const auto tMat5 = rnd<tMat>(S.M,S.N);
	const auto tMat6 = T(tMat5);
			
	CPPUNIT_ASSERT_EQUAL(tMat6.M(),tMat5.N());
	CPPUNIT_ASSERT_EQUAL(tMat6.N(),tMat5.M());
			
	for (size_t n=0; n<tMat5.N(); ++n)
		for (size_t m=0; m<tMat5.M(); ++m)
			CPPUNIT_ASSERT_EQUAL(CT(tMat5(m,n)),CT(std::conj(tMat6(n,m))));
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_Iz_Inz() {
	const size_t M = genRndST();
	const size_t N = genRndST();
	auto tMat1 = rnd<tMat>(M,N);
	tMat1 = (tMat1.gt(mean(tMat1)));

	const auto I_z = Iz(tMat1);
	CPPUNIT_ASSERT(std::all_of(I_z.cbegin(),I_z.cend(),
			[&tMat1](const Ipair& i)->bool{return i.m<tMat1.M();}));
	CPPUNIT_ASSERT(std::all_of(I_z.cbegin(),I_z.cend(),
			[&tMat1](const Ipair& i)->bool{return i.n<tMat1.N();}));
	
	const auto I_nz = Inz(tMat1);
	CPPUNIT_ASSERT(std::all_of(I_nz.cbegin(),I_nz.cend(),
			[&tMat1](const Ipair& i)->bool{return i.m<tMat1.M();}));
	CPPUNIT_ASSERT(std::all_of(I_nz.cbegin(),I_nz.cend(),
			[&tMat1](const Ipair& i)->bool{return i.n<tMat1.N();}));

	for (size_t i=0; i!=I_z.size(); ++i)
		CPPUNIT_ASSERT_EQUAL(TT(0.0),tMat1(I_z[i].m,I_z[i].n));
	for (size_t i=0; i!=I_nz.size(); ++i)
		CPPUNIT_ASSERT(TT(0.0)!=tMat1(I_nz[i].m,I_nz[i].n));

	// check I_z and I_nz are mutually exclusive
	CPPUNIT_ASSERT(std::all_of(I_z.cbegin(),I_z.cend(),[&I_nz](const Ipair& i)->bool{
		return std::find_if(I_nz.cbegin(),I_nz.cend(),[&i](const Ipair& j)->bool{
				return j.m==i.m && j.n==i.n;
			})==I_nz.cend();
	}));

	// check all indices are included in I_z or I_nz
	CPPUNIT_ASSERT_EQUAL(tMat1.L(),I_z.size()+I_nz.size());
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_conj() {
	const size_t M = genRndST();
	const size_t N = genRndST();
	const auto tMat1 = rnd<tMat>(M,N);

	auto tMat2 = tMat1; 
	std::for_each(tMat2.begin(),tMat2.end(),[](TT& i){lm__::ops::assign(i,std::conj(i));});
	const auto tMat3 = conj(tMat1);

	CPPUNIT_ASSERT(size(tMat2)==size(tMat3));
	CPPUNIT_ASSERT(tMat2==tMat3);
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_real() {
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	const auto re = rnd<fMat>(M,N);
	const auto im = rnd<fMat>(M,N);
	const auto tMat1 = tMat(re,im);

	CPPUNIT_ASSERT(re==real(tMat1));
}

template<class TT, class FT, class CT>
void test_tFn_conversion<TT,FT,CT>::test_imag() {
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	const auto re = rnd<fMat>(M,N);
	const auto im = rnd<fMat>(M,N);
	const auto tMat1 = tMat(re,im);

	if (!tMat1.cpx()) CPPUNIT_ASSERT(zeros<fMat>(M,N)==imag(tMat1));
	else CPPUNIT_ASSERT(im==imag(tMat1));
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tFn_conversion<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_msum", &test_tFn_conversion<TT,FT,CT>::test_msum));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_nsum", &test_tFn_conversion<TT,FT,CT>::test_nsum));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_mprod", &test_tFn_conversion<TT,FT,CT>::test_mprod));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_nprod", &test_tFn_conversion<TT,FT,CT>::test_nprod));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_mmin", &test_tFn_conversion<TT,FT,CT>::test_mmin));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_nmin", &test_tFn_conversion<TT,FT,CT>::test_nmin));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_mmax", &test_tFn_conversion<TT,FT,CT>::test_mmax));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_nmax", &test_tFn_conversion<TT,FT,CT>::test_nmax));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_mmean", &test_tFn_conversion<TT,FT,CT>::test_mmean));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_nmean", &test_tFn_conversion<TT,FT,CT>::test_nmean));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_mnormsq", &test_tFn_conversion<TT,FT,CT>::test_mnormsq));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_nnormsq", &test_tFn_conversion<TT,FT,CT>::test_nnormsq));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_mnorm", &test_tFn_conversion<TT,FT,CT>::test_mnorm));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_nnorm", &test_tFn_conversion<TT,FT,CT>::test_nnorm));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_mdot", &test_tFn_conversion<TT,FT,CT>::test_mdot));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_ndot", &test_tFn_conversion<TT,FT,CT>::test_ndot));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_mdotu", &test_tFn_conversion<TT,FT,CT>::test_mdotu));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_ndotu", &test_tFn_conversion<TT,FT,CT>::test_ndotu));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_inv", &test_tFn_conversion<TT,FT,CT>::test_inv));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_diag", &test_tFn_conversion<TT,FT,CT>::test_diag));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_R_C", &test_tFn_conversion<TT,FT,CT>::test_R_C));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_T", &test_tFn_conversion<TT,FT,CT>::test_T));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_Iz_Inz", &test_tFn_conversion<TT,FT,CT>::test_Iz_Inz));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_conj", &test_tFn_conversion<TT,FT,CT>::test_conj));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_real", &test_tFn_conversion<TT,FT,CT>::test_real));
	suite->addTest(new CppUnit::TestCaller<test_tFn_conversion>(
		"test_imag", &test_tFn_conversion<TT,FT,CT>::test_imag));
	
	return suite;
}
