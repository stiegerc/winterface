// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_matrix_arithmetic.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


template<class TT, class FT, class CT>
void test_tMat_matrix_arithmetic<TT,FT,CT>::test_prod() {
	//const size_t M = genRndST();
	//const size_t N = genRndST();
	//const size_t K = genRndST();
	const size_t M = 2;
	const size_t N = 3;
	const size_t K = 4;
	const auto tMat1 = rnd<tMat>(M,N);

	// with fMat
	{
		const auto ftMat = rnd<fMat>(N,K);
		
		const auto sol = tMat1.prod(ftMat);
		CPPUNIT_ASSERT_EQUAL(M,sol.M());
		CPPUNIT_ASSERT_EQUAL(K,sol.N());
		if (tMat1.cpx()) CPPUNIT_ASSERT(sol.cpx());
		else CPPUNIT_ASSERT(!sol.cpx());

		tMat ck(sol.M(),sol.N());
		auto i=ck.begin();
		for (auto ki=ftMat.cBegin(),ke=ftMat.cEnd(); ki!=ke; ++ki)
			for (auto mi=tMat1.rBegin(),me=tMat1.rEnd(); mi!=me; ++mi,++i)
				*i = mi->prod(*ki);
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}
	
	// with cMat
	{
		const auto ctMat = rnd<cMat>(N,K);
		
		const auto sol = tMat1.prod(ctMat);
		CPPUNIT_ASSERT_EQUAL(M,sol.M());
		CPPUNIT_ASSERT_EQUAL(K,sol.N());
		CPPUNIT_ASSERT(sol.cpx());

		cMat ck(sol.M(),sol.N());
		auto i=ck.begin();
		for (auto ki=ctMat.cBegin(),ke=ctMat.cEnd(); ki!=ke; ++ki)
			for (auto mi=tMat1.rBegin(),me=tMat1.rEnd(); mi!=me; ++mi,++i)
				*i = mi->prod(*ki);
		
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}

	// with fRow
	{
		const auto ftMat = rnd<fMat>(1,K);
		
		const size_t n = genRndST(0,N-1);

		const auto sol = tMat1.cGet(n,1).prod(ftMat.rFront());
		const auto ck = tMat1.cGet(n,1).prod(ftMat.rGet(0,1));

		CPPUNIT_ASSERT(size(ck)==size(sol));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}
	
	// with cRow
	{
		const auto ctMat = rnd<cMat>(1,K);
		
		const size_t n = genRndST(0,N-1);

		const auto sol = tMat1.cGet(n,1).prod(ctMat.rFront());
		const auto ck = tMat1.cGet(n,1).prod(ctMat.rGet(0,1));

		CPPUNIT_ASSERT(size(ck)==size(sol));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}

	// with fCol
	{
		const auto ftMat = rnd<fMat>(N,K);

		const size_t n = genRndST(0,K-1);

		const auto sol = tMat1.prod(*(ftMat.ccBegin()+n));
		const auto ck = tMat1.prod(ftMat.cGet(n,1));
		
		CPPUNIT_ASSERT(size(ck)==size(sol));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}
	
	// with fCol
	{
		const auto ctMat = rnd<cMat>(N,K);

		const size_t n = genRndST(0,K-1);

		const auto sol = tMat1.prod(*(ctMat.ccBegin()+n));
		const auto ck = tMat1.prod(ctMat.cGet(n,1));
		
		CPPUNIT_ASSERT(size(ck)==size(sol));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}
}

template<class TT, class FT, class CT>
void test_tMat_matrix_arithmetic<TT,FT,CT>::test_leftDivide() {

	const size_t M = genRndST(5,10);
	const size_t N = genRndST(5,10);
	const auto tMat1 = rnd_sqb<tMat>(M);
	
	// leftDivide with fMat
	{
		const auto rhs = rnd<fMat>(M,N,RE__(-1.0),RE__(1.0));
		const auto sol = tMat1.leftDivide(rhs);
		
		auto tMat2 = tMat1; tMat2.inv();
		const auto ck = tMat2.prod(rhs);

		CPPUNIT_ASSERT(size(ck)==size(sol));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}
	
	// leftDivide with cMat
	{
		const auto rhs = rnd<cMat>(M,N,RE__(-1.0),RE__(1.0));
		const auto sol = tMat1.leftDivide(rhs);
		
		auto tMat2 = tMat1; tMat2.inv();
		const auto ck = tMat2.prod(rhs);

		CPPUNIT_ASSERT(size(ck)==size(sol));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}

	// leftDivide with fRow
	{
		const auto tMat1 = rnd<tMat>(1,1);
		const auto rhs = rnd<fMat>(1,N,RE__(-1.0),RE__(1.0));
		
		const auto i = rhs.crBegin();
		
		const auto sol = tMat1.leftDivide(*i);
		auto tMat2 = tMat1; tMat2.inv();
		const auto ck = tMat2.prod(*i);

		CPPUNIT_ASSERT(size(ck)==size(sol));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}
	
	// leftDivide with cRow
	{
		const auto tMat1 = rnd<tMat>(1,1);
		const auto rhs = rnd<cMat>(1,N,RE__(-1.0),RE__(1.0));
		
		const auto i = rhs.crBegin();
		
		const auto sol = tMat1.leftDivide(*i);
		auto tMat2 = tMat1; tMat2.inv();
		const auto ck = tMat2.prod(*i);

		CPPUNIT_ASSERT(size(ck)==size(sol));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}
	
	// leftDivide with fCol
	{
		const auto rhs = rnd<fMat>(M,1,RE__(-1.0),RE__(1.0));
		const auto i = rhs.ccBegin();
		
		const auto sol = tMat1.leftDivide(*i);
		auto tMat2 = tMat1; tMat2.inv();
		const auto ck = tMat2.prod(*i);

		CPPUNIT_ASSERT(size(ck)==size(sol));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}
	
	// leftDivide with cCol
	{
		const auto rhs = rnd<cMat>(M,1,RE__(-1.0),RE__(1.0));
		const auto i = rhs.ccBegin();
		
		const auto sol = tMat1.leftDivide(*i);
		auto tMat2 = tMat1; tMat2.inv();
		const auto ck = tMat2.prod(*i);

		CPPUNIT_ASSERT(size(ck)==size(sol));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}
}

template<class TT, class FT, class CT>
void test_tMat_matrix_arithmetic<TT,FT,CT>::test_leftDivideEq() {

	const size_t M = genRndST(5,10);
	const size_t N = genRndST(5,10);
	const auto tMat1 = rnd_sqb<tMat>(M);
	
	// leftDivideEq with fMat
	{
		auto rhs = rnd<fMat>(M,N,RE__(-1.0),RE__(1.0));
		
		auto tMat2 = tMat1; tMat2.inv();
		const fMat ck = tMat2.prod(rhs);
		
		const auto sol = tMat1.leftDivideEq(rhs);
		
		CPPUNIT_ASSERT(size(ck)==size(sol));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}
	
	// leftDivideEq with cMat
	{
		auto rhs = rnd<cMat>(M,N,RE__(-1.0),RE__(1.0));
		
		auto tMat2 = tMat1; tMat2.inv();
		const cMat ck = tMat2.prod(rhs);
		
		const auto sol = tMat1.leftDivideEq(rhs);
		

		CPPUNIT_ASSERT(size(ck)==size(sol));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],sol[i],delta);
	}

	// leftDivideEq with fRow
	{
		const auto tMat1 = rnd<tMat>(1,1);
		auto rhs = rnd<fMat>(1,N,RE__(-1.0),RE__(1.0));
		
		auto ck = rhs;
		tMat1.leftDivideEq(ck);
		tMat1.leftDivideEq(*rhs.rBegin());
		
		CPPUNIT_ASSERT(size(ck)==size(rhs));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],rhs[i],delta);
	}
	
	// leftDivideEq with cRow
	{
		const auto tMat1 = rnd<tMat>(1,1);
		auto rhs = rnd<cMat>(1,N,RE__(-1.0),RE__(1.0));
		
		auto ck = rhs;
		tMat1.leftDivideEq(ck);
		tMat1.leftDivideEq(*rhs.rBegin());
		
		CPPUNIT_ASSERT(size(ck)==size(rhs));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],rhs[i],delta);
	}
	// leftDivideEq with fCol
	{
		auto rhs = rnd<fMat>(M,1,RE__(-1.0),RE__(1.0));
		
		auto ck = rhs;
		tMat1.leftDivideEq(ck);
		tMat1.leftDivideEq(*rhs.cBegin());
		
		CPPUNIT_ASSERT(size(ck)==size(rhs));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],rhs[i],delta);
	}
	
	// leftDivideEq with cCol
	{
		auto rhs = rnd<cMat>(M,1,RE__(-1.0),RE__(1.0));
		
		auto ck = rhs;
		tMat1.leftDivideEq(ck);
		tMat1.leftDivideEq(*rhs.cBegin());
		
		CPPUNIT_ASSERT(size(ck)==size(rhs));
		for (size_t i=0; i!=ck.L(); ++i)
			CPPUNIT_ASSERT_DELTA(ck[i],rhs[i],delta);
	}
}



template<class TT, class FT, class CT>
CppUnit::Test* test_tMat_matrix_arithmetic<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tMat_matrix_arithmetic>(
		"test_prod", &test_tMat_matrix_arithmetic<TT,FT,CT>::test_prod));
	suite->addTest(new CppUnit::TestCaller<test_tMat_matrix_arithmetic>(
		"test_leftDivide", &test_tMat_matrix_arithmetic<TT,FT,CT>::test_leftDivide));
	suite->addTest(new CppUnit::TestCaller<test_tMat_matrix_arithmetic>(
		"test_leftDivideEq", &test_tMat_matrix_arithmetic<TT,FT,CT>::test_leftDivideEq));
	
	return suite;
}
