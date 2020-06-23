// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_comparison.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


template<class TT, class FT, class CT>
void test_tMat_comparison<TT,FT,CT>::test_msize_function() {

	const auto S = genRndMEST();
	const auto tMat1 = rnd<tMat>(S.M,S.N);
	
	CPPUNIT_ASSERT(msize(tMat1)==msize(tMat1));
	
	const auto tfMat1 = rnd<fMat>(S.M,S.N);
	const auto tfMat2 = rnd<fMat>(S.N,S.M);
	CPPUNIT_ASSERT(msize(tMat1)==msize(tfMat1));
	CPPUNIT_ASSERT(msize(tMat1)!=msize(tfMat2));

	const auto tcMat1 = rnd<cMat>(S.M,S.N);
	const auto tcMat2 = rnd<cMat>(S.N,S.M);
	CPPUNIT_ASSERT(msize(tMat1)==msize(tcMat1));
	CPPUNIT_ASSERT(msize(tMat1)!=msize(tcMat2));
}

template<class TT, class FT, class CT>
void test_tMat_comparison<TT,FT,CT>::test_operator_equal_unequal_() {
	
	const auto S = genRndMEST();
	const auto tMat1 = rnd<tMat>(S.M,S.N);
	
	{
		fMat tfMat1 = rnd<fMat>(S.M,S.N);
		CPPUNIT_ASSERT(tMat1!=tfMat1);

		tfMat1 = tMat1;
		CPPUNIT_ASSERT(real(tMat1)==tfMat1);

		CPPUNIT_ASSERT(tMat1!=tfMat1[0]);
		auto tMat2 = tfMat1;
		for (auto& i: tMat2) i=tMat2[0];
		CPPUNIT_ASSERT(tMat2==tMat2[0]);
	}
	{
		cMat tcMat1 = rnd<fMat>(S.M,S.N);
		CPPUNIT_ASSERT(tMat1!=tcMat1);

		tcMat1 = tMat1;
		CPPUNIT_ASSERT(tMat1==tcMat1);

		CPPUNIT_ASSERT(tMat1!=tcMat1[0]);
		auto tMat2 = tcMat1;
		for (auto& i: tMat2) i=tMat2[0];
		CPPUNIT_ASSERT(tMat2==tMat2[0]);
	}
	{
		fMat tfMat1 = rnd<fMat>(S.M,S.N);
		auto j = tfMat1.begin();
		for (auto i=tMat1.cbegin(),e=tMat1.cend(); i!=e; ++i,++j)
			if (lm__::ops::lt(*i,*j)) *j = std::real(*i - FT(1.0));

		CPPUNIT_ASSERT(tMat1>tfMat1 || tMat1>=tfMat1);
		CPPUNIT_ASSERT(!(tMat1<tfMat1) || !(tMat1==tfMat1));
	}
	{
		fMat tcMat1 = rnd<cMat>(S.M,S.N);
		auto j = tcMat1.begin();
		for (auto i=tMat1.cbegin(),e=tMat1.cend(); i!=e; ++i,++j)
			if (lm__::ops::lt(*i,*j)) *j = std::real(*i - FT(1.0));
		
		CPPUNIT_ASSERT(tMat1>tcMat1 || tMat1>=tcMat1);
		CPPUNIT_ASSERT(!(tMat1<tcMat1) || !(tMat1==tcMat1));
	}


#ifndef NTOLERANT__
	{
		const auto ttMat1 = tMat(real(tMat1));
		
		fMat tfMat1 = tMat1+mtol()*FT(.5);
		CPPUNIT_ASSERT(ttMat1==tfMat1);
		
		fMat tfMat2 = tMat1+mtol()*FT(2.0);
		CPPUNIT_ASSERT(ttMat1!=tfMat2);

		auto tMat2 = tfMat1;
		for (auto& i: tMat2) i=tMat2[0]+mtol()*FT(.5);
		CPPUNIT_ASSERT(tMat2==tMat2[0]);
		
		auto tMat3 = tfMat1;
		for (auto& i: tMat2) i=tMat3[0]+mtol()*FT(2.0);
		CPPUNIT_ASSERT(tMat2!=tMat3[0]);
	}
	{
		cMat tcMat1 = tMat1+mtol()*FT(.5);
		CPPUNIT_ASSERT(tMat1==tcMat1);
		
		cMat tcMat2 = tMat1+mtol()*FT(2.0);
		CPPUNIT_ASSERT(tMat1!=tcMat2);

		auto tMat2 = tcMat1;
		for (auto& i: tMat2) i=tMat2[0]+mtol()*FT(.5);
		CPPUNIT_ASSERT(tMat2==tMat2[0]);
		
		auto tMat3 = tcMat1;
		for (auto& i: tMat2) i=tMat3[0]+mtol()*FT(2.0);
		CPPUNIT_ASSERT(tMat2!=tMat3[0]);
	}
#endif
}

template<class TT, class FT, class CT>
void test_tMat_comparison<TT,FT,CT>::test_eq_neq_leq_() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// eq, neq
	auto tMat1 = rnd<tMat>(M,N); tMat1[1]=tMat1[0]+TT(1.0);	
	{
		const auto tMat2 = tMat1.eq(tMat1);
		CPPUNIT_ASSERT(size(tMat2)==size(tMat1));
		CPPUNIT_ASSERT(all(tMat2));
		
		const auto tMat3 = tMat1.neq(tMat1);
		CPPUNIT_ASSERT(size(tMat2)==size(tMat1));
		CPPUNIT_ASSERT(all(~tMat3));

		const auto tMat4 = ncat(tMat1,-rnd<tMat>(M,N)-1.0);
		const auto tMat5 = ncat(tMat1,tMat1);
		
		const auto tMat6 = tMat4.eq(tMat5);
		const auto sol6 = ncat(ones<tMat>(M,N),zeros<tMat>(M,N));
		CPPUNIT_ASSERT(tMat6==sol6);
		
		const auto tMat7 = tMat4.neq(tMat5);
		const auto sol7 = ncat(zeros<tMat>(M,N),ones<tMat>(M,N));
		CPPUNIT_ASSERT(tMat7==sol7);
	
		tMat tMat8 = real(tMat1);
		for (auto& i: tMat8)
			i=std::real(tMat1[0]);
		const auto tMat9 = tMat8.eq(std::real(tMat1[0]));
		const auto tMat10 = tMat8.eq(std::real(tMat1[1]));
		CPPUNIT_ASSERT(all(tMat9));
		CPPUNIT_ASSERT(all(~tMat10));

		tMat tMat11 = tMat1;
		for (auto& i: tMat11)
			i=tMat1[0];
		const auto tMat12 = tMat11.eq(tMat1[0]);
		const auto tMat13 = tMat11.eq(tMat1[1]);
		CPPUNIT_ASSERT(all(tMat12));
		CPPUNIT_ASSERT(all(~tMat13));
	}
#ifndef NTOLERANT__
	{
		CPPUNIT_ASSERT(all(tMat1.eq(tMat1+mtol()*FT(.5))));
		CPPUNIT_ASSERT(all(tMat1.neq(tMat1+mtol()*FT(2.0))));
	}
#endif


	// <,>
	{
		const auto tMat1 = rnd<tMat>(M,N,.6,1.0);
		const auto ftMat1 = rnd<fMat>(M,N,.0,.5);
		CPPUNIT_ASSERT(all(tMat1.gt(ftMat1)));
		CPPUNIT_ASSERT(all(tMat1.gt(FT(.5))));
		
		const auto tMat2 = rnd<tMat>(M,N,.0,.5);
		const auto ftMat2 = rnd<fMat>(M,N,.51,1.0);
		CPPUNIT_ASSERT(all(tMat2.lt(ftMat2)));
		CPPUNIT_ASSERT(all(tMat2.lt(FT(1.0))));
	}
	{
		const auto tMat1 = rnd<tMat>(M,N,.6,1.0);
		const auto ctMat1 = rnd<cMat>(M,N,.0,.5);
		CPPUNIT_ASSERT(all(tMat1.gt(ctMat1)));
		CPPUNIT_ASSERT(all(tMat1.gt(CT(.5))));
		
		const auto tMat2 = rnd<tMat>(M,N,.0,.5);
		const auto ctMat2 = rnd<cMat>(M,N,.51,1.0);
		CPPUNIT_ASSERT(all(tMat2.lt(ctMat2)));
		CPPUNIT_ASSERT(all(tMat2.lt(CT(1.0))));
	}
#ifndef NTOLERANT__
	{
		const auto tMat1 = rnd<tMat>(M,N);
		CPPUNIT_ASSERT(all((tMat1-mtol()*2.0).lt(tMat1)));
		CPPUNIT_ASSERT(all((~(tMat1-mtol()*.5).lt(tMat1))));
		CPPUNIT_ASSERT(all((tMat1+mtol()*2.0).gt(tMat1)));
		CPPUNIT_ASSERT(all((~(tMat1+mtol()*.5).gt(tMat1))));
	}
#endif


	// <=,>=
	{
		auto tMat1 = rnd<tMat>(M,N,.5,1.0); tMat1[0]=TT(.5);
		auto ftMat1 = rnd<fMat>(M,N,.0,.5); ftMat1[0]=FT(.5);
		CPPUNIT_ASSERT(all(tMat1.geq(ftMat1)));
		CPPUNIT_ASSERT(all(tMat1.geq(FT(.5))));
		
		auto tMat2 = rnd<tMat>(M,N,.0,.5); tMat2[0]=TT(.5);
		auto ftMat2 = rnd<fMat>(M,N,.5,1.0); ftMat2[0]=FT(.5);
		CPPUNIT_ASSERT(all(tMat2.leq(ftMat2)));
		CPPUNIT_ASSERT(all(tMat2.leq(FT(1.0))));
	}
	{
		auto tMat1 = rnd<tMat>(M,N,.5,1.0); tMat1[0]=TT(.5);
		auto ctMat1 = rnd<cMat>(M,N,.0,.5); ctMat1[0]=CT(.5);
		CPPUNIT_ASSERT(all(tMat1.geq(ctMat1)));
		CPPUNIT_ASSERT(all(tMat1.geq(CT(.5))));
		
		auto tMat2 = rnd<tMat>(M,N,.0,.5); tMat2[0]=TT(.5);
		auto ctMat2 = rnd<cMat>(M,N,.5,1.0); ctMat2[0]=CT(.5);
		CPPUNIT_ASSERT(all(tMat2.leq(ctMat2)));
		CPPUNIT_ASSERT(all(tMat2.leq(CT(1.0))));
	}
#ifndef NTOLERANT__
	{
		const auto tMat1 = rnd<tMat>(M,N);
		CPPUNIT_ASSERT((all((tMat1-mtol()*2.0).leq(tMat1))));
		CPPUNIT_ASSERT((all((tMat1-mtol()*.5).leq(tMat1))));
		CPPUNIT_ASSERT((all((tMat1+mtol()*2.0).geq(tMat1))));
		CPPUNIT_ASSERT((all((tMat1+mtol()*.5).geq(tMat1))));
	}
#endif

}


template<class TT, class FT, class CT>
CppUnit::Test* test_tMat_comparison<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tMat_comparison>(
		"test_size_function", &test_tMat_comparison<TT,FT,CT>::test_msize_function));
	suite->addTest(new CppUnit::TestCaller<test_tMat_comparison>(
		"test_operator_equal_unequal_", &test_tMat_comparison<TT,FT,CT>::test_operator_equal_unequal_));
	suite->addTest(new CppUnit::TestCaller<test_tMat_comparison>(
		"test_eq_neq_leq_", &test_tMat_comparison<TT,FT,CT>::test_eq_neq_leq_));
	
	return suite;
}
