// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_r_tVecItr_all.h"
#include "testTools.h"
#include "lm_tMat.h"
#include <iostream>


using namespace lm__::test;

template<class TT, class FT, class CT, class VT>
void test_r_tVecItr_all<TT,FT,CT,VT>::test_ctor_assign() {

	const size_t M = genRndST();
	const size_t N = genRndST();
	auto tMat1 = rnd<tMat>(M,N);
	const ptrdiff_t i = genRndST(0,N-1);

	// default
	r_tVecItr tItr1;
	CPPUNIT_ASSERT(!tItr1.ptr());
	CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr1.i());

	// one argument
	r_tVecItr tItr2(&tMat1);
	CPPUNIT_ASSERT_EQUAL(&tMat1,tItr2.ptr());
	CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr2.i());

	// two arguments
	r_tVecItr tItr3(&tMat1,i);
	CPPUNIT_ASSERT_EQUAL(&tMat1,tItr3.ptr());
	CPPUNIT_ASSERT_EQUAL(i,tItr3.i());

	// copy
	r_tVecItr tItr4(tItr3);
	CPPUNIT_ASSERT_EQUAL(tItr3.ptr(),tItr4.ptr());
	CPPUNIT_ASSERT_EQUAL(tItr3.i(),tItr4.i());

	// assignment
	tItr2 = tItr3;
	CPPUNIT_ASSERT_EQUAL(tItr3.ptr(),tItr2.ptr());
	CPPUNIT_ASSERT_EQUAL(tItr3.i(),tItr2.i());

	// from r_tVecItr
	{
		r_tVecItr tItr1(&tMat1,i);
		r_tVecItr tItr2(tItr1);
		CPPUNIT_ASSERT_EQUAL(tItr1.ptr(),tItr2.ptr());
		CPPUNIT_ASSERT_EQUAL(tItr1.i(),tItr2.i());
	
		tItr2=tItr1;
		CPPUNIT_ASSERT_EQUAL(tItr1.ptr(),tItr2.ptr());
		CPPUNIT_ASSERT_EQUAL(tItr1.i(),tItr2.i());
	}
	
	// from tVecItr
	{
		tVecItr tItr1(&tMat1,i);
		r_tVecItr tItr2(tItr1);
		CPPUNIT_ASSERT_EQUAL(tItr1.ptr(),tItr2.ptr());
		CPPUNIT_ASSERT_EQUAL(tItr1.i()-1,tItr2.i());
	
		tItr2=tItr1;
		CPPUNIT_ASSERT_EQUAL(tItr1.ptr(),tItr2.ptr());
		CPPUNIT_ASSERT_EQUAL(tItr1.i()-1,tItr2.i());
	}
}

template<class TT, class FT, class CT, class VT>
void test_r_tVecItr_all<TT,FT,CT,VT>::test_swap() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	auto tMat1 = rnd<tMat>(M,N);
	const ptrdiff_t i1 = genRndST(0,N-1);
	
	auto tMat2 = rnd<tMat>(M,N);
	const ptrdiff_t i2 = genRndST(0,N-1);

	r_tVecItr tItr1(&tMat1,i1);
	r_tVecItr tItr2(&tMat2,i2);

	swap(tItr1,tItr2);

	CPPUNIT_ASSERT_EQUAL(&tMat1,tItr2.ptr());
	CPPUNIT_ASSERT_EQUAL(i1,tItr2.i());
	CPPUNIT_ASSERT_EQUAL(&tMat2,tItr1.ptr());
	CPPUNIT_ASSERT_EQUAL(i2,tItr1.i());
}

template<class TT, class FT, class CT, class VT>
void test_r_tVecItr_all<TT,FT,CT,VT>::test_comparison() {

	const size_t M = genRndST();
	const size_t N = genRndST();
	
	auto tMat1 = rnd<tMat>(M,N);
	const ptrdiff_t i1 = genRndST(0,N/2);
	const ptrdiff_t i2 = genRndST(N/2+1,N-1);

	r_tVecItr tItr1(&tMat1,i1);
	
	// to cr_tVecItr
	{
		cr_tVecItr tItr2(&tMat1,i2);

		CPPUNIT_ASSERT(tItr1==tItr1);
		CPPUNIT_ASSERT(tItr1!=tItr2);
		CPPUNIT_ASSERT(tItr1>tItr2);
		CPPUNIT_ASSERT(tItr2<tItr1);
		CPPUNIT_ASSERT(tItr1>=tItr2);
		CPPUNIT_ASSERT(tItr1>=tItr1);
		CPPUNIT_ASSERT(tItr2<=tItr1);
		CPPUNIT_ASSERT(tItr2<=tItr2);
	}

	// to r_tVecItr
	{
		r_tVecItr tItr2(const_cast<tMat*>(&tMat1),i2);

		CPPUNIT_ASSERT(tItr1==tItr1);
		CPPUNIT_ASSERT(tItr1!=tItr2);
		CPPUNIT_ASSERT(tItr1>tItr2);
		CPPUNIT_ASSERT(tItr2<tItr1);
		CPPUNIT_ASSERT(tItr1>=tItr2);
		CPPUNIT_ASSERT(tItr1>=tItr1);
		CPPUNIT_ASSERT(tItr2<=tItr1);
		CPPUNIT_ASSERT(tItr2<=tItr2);
	}
}

template<class TT, class FT, class CT, class VT>
void test_r_tVecItr_all<TT,FT,CT,VT>::test_arithmetic_difference() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	auto tMat1 = rnd<tMat>(M,N);

	r_tVecItr tItr1(&tMat1,0);
	
	// increment, decrement
	{
		r_tVecItr tItr2(&tMat1,0);

		auto tItr3 = tItr2++;
		CPPUNIT_ASSERT(tItr1==tItr1);
		CPPUNIT_ASSERT(tItr2>tItr1);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(1),tItr2-tItr3);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(-1),tItr3-tItr2);
		CPPUNIT_ASSERT_EQUAL(std::abs(tItr3-tItr1),distance(tItr1,tItr3));

		auto tItr4 = ++tItr2;
		CPPUNIT_ASSERT(tItr2==tItr4);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr2-tItr4);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr4-tItr2);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(2),tItr2-tItr3);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(-2),tItr3-tItr2);
		CPPUNIT_ASSERT_EQUAL(std::abs(tItr4-tItr1),distance(tItr1,tItr4));

		auto tItr5 = tItr2--;
		CPPUNIT_ASSERT(tItr5>tItr2);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(1),tItr5-tItr2);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(-1),tItr2-tItr5);
		CPPUNIT_ASSERT_EQUAL(std::abs(tItr5-tItr1),distance(tItr1,tItr5));

		auto tItr6 = --tItr2;
		CPPUNIT_ASSERT(tItr2==tItr6);
		CPPUNIT_ASSERT(tItr1==tItr2);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr2-tItr6);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr6-tItr2);
		CPPUNIT_ASSERT_EQUAL(std::abs(tItr6-tItr1),distance(tItr1,tItr6));
	}

	// += -= + -
	{
		r_tVecItr tItr2(&tMat1,0);

		const ptrdiff_t l = genRndST(0,N-1);

		auto tItr3 = tItr2+l;
		CPPUNIT_ASSERT_EQUAL(l,tItr3-tItr2);
		CPPUNIT_ASSERT_EQUAL(-l,tItr2-tItr3);
		CPPUNIT_ASSERT_EQUAL(std::abs(tItr3-tItr1),distance(tItr1,tItr3));

		auto tItr4 = tItr3-l;
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr4-tItr2);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr2-tItr4);
		CPPUNIT_ASSERT_EQUAL(std::abs(tItr4-tItr1),distance(tItr1,tItr4));

		tItr2+=l;
		CPPUNIT_ASSERT_EQUAL(l,tItr2-tItr1);
		CPPUNIT_ASSERT_EQUAL(-l,tItr1-tItr2);
		
		tItr2-=l;
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr1-tItr2);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr2-tItr1);
	}
}


template<class TT, class FT, class CT, class VT>
CppUnit::Test* test_r_tVecItr_all<TT,FT,CT,VT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_r_tVecItr_all>(
		"test_ctor_assign", &test_r_tVecItr_all<TT,FT,CT,VT>::test_ctor_assign));
	suite->addTest(new CppUnit::TestCaller<test_r_tVecItr_all>(
		"test_swap", &test_r_tVecItr_all<TT,FT,CT,VT>::test_swap));
	suite->addTest(new CppUnit::TestCaller<test_r_tVecItr_all>(
		"test_comparison", &test_r_tVecItr_all<TT,FT,CT,VT>::test_comparison));
	suite->addTest(new CppUnit::TestCaller<test_r_tVecItr_all>(
		"test_dereference", &test_r_tVecItr_all<TT,FT,CT,VT>::test_dereference));
	suite->addTest(new CppUnit::TestCaller<test_r_tVecItr_all>(
		"test_arithmetic_difference", &test_r_tVecItr_all<TT,FT,CT,VT>::test_arithmetic_difference));
	
	return suite;
}
