// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_c_tItr_all.h"
#include "testTools.h"
#include <iostream>


using namespace lm__::test;

template<class TT>
void test_c_tItr_all<TT>::test_ctor_assign() {
	
	// lambda to generate const random range
	auto gen = [](const size_t L) -> const TT* {
		TT* res = new TT [L];
		rnd(res,L);
		return res;
	};

	const size_t L = genRndST();
	const TT* rg = gen(L);
	const size_t incr = genRndST();
	
	// default
	c_tItr tItr1;
	CPPUNIT_ASSERT(!tItr1.data());
	CPPUNIT_ASSERT_EQUAL(size_t(1),tItr1.incr());

	// one argument
	c_tItr tItr2(rg);
	CPPUNIT_ASSERT_EQUAL(rg,tItr2.data());
	CPPUNIT_ASSERT_EQUAL(size_t(1),tItr2.incr());

	// two arguments
	c_tItr tItr3(rg,incr);
	CPPUNIT_ASSERT_EQUAL(rg,tItr3.data());
	CPPUNIT_ASSERT_EQUAL(incr,tItr3.incr());

	// copy
	c_tItr tItr4(tItr3);
	CPPUNIT_ASSERT_EQUAL(tItr3.data(),tItr4.data());
	CPPUNIT_ASSERT_EQUAL(tItr3.incr(),tItr4.incr());

	// assignment
	tItr2 = tItr3;
	CPPUNIT_ASSERT_EQUAL(tItr3.data(),tItr2.data());
	CPPUNIT_ASSERT_EQUAL(tItr3.incr(),tItr2.incr());

	// from itr
	{
		c_tItr tItr1(const_cast<TT*>(rg),incr);
		c_tItr tItr2(tItr1);
		CPPUNIT_ASSERT_EQUAL(const_cast<const TT*>(tItr1.data()),tItr2.data());
		CPPUNIT_ASSERT_EQUAL(tItr1.incr(),tItr2.incr());
		
		tItr2=tItr1;
		CPPUNIT_ASSERT_EQUAL(const_cast<const TT*>(tItr1.data()),tItr2.data());
		CPPUNIT_ASSERT_EQUAL(tItr1.incr(),tItr2.incr());
	}
	
	delete[] rg;
}

template<class TT>
void test_c_tItr_all<TT>::test_swap() {
	
	// lambda to generate const random range
	auto gen = [](const size_t L) -> const TT* {
		TT* res = new TT [L];
		rnd(res,L);
		return res;
	};

	const size_t L1 = genRndST();
	const TT* rg1 = gen(L1);
	const size_t incr1 = genRndST();
	
	const size_t L2 = genRndST();
	const TT* rg2 = gen(L2);
	const size_t incr2 = genRndST();

	c_tItr tItr1(rg1,incr1);
	c_tItr tItr2(rg2,incr2);

	swap(tItr1,tItr2);

	CPPUNIT_ASSERT_EQUAL(rg1,tItr2.data());
	CPPUNIT_ASSERT_EQUAL(incr1,tItr2.incr());
	CPPUNIT_ASSERT_EQUAL(rg2,tItr1.data());
	CPPUNIT_ASSERT_EQUAL(incr2,tItr1.incr());

	delete[] rg1;
	delete[] rg2;
}

template<class TT>
void test_c_tItr_all<TT>::test_comparison() {

	// lambda to generate const random range
	auto gen = [](const size_t L) -> const TT* {
		TT* res = new TT [L];
		rnd(res,L);
		return res;
	};

	const size_t L = genRndST();
	const TT* rg = gen(L);
	const size_t incr = genRndST();

	c_tItr tItr1(rg,incr);
	
	// to c_tItr
	{
		c_tItr tItr2(rg+L/2,incr);

		CPPUNIT_ASSERT(tItr1==tItr1);
		CPPUNIT_ASSERT(tItr1!=tItr2);
		CPPUNIT_ASSERT(tItr1<tItr2);
		CPPUNIT_ASSERT(tItr2>tItr1);
		CPPUNIT_ASSERT(tItr1<=tItr2);
		CPPUNIT_ASSERT(tItr1<=tItr1);
		CPPUNIT_ASSERT(tItr2>=tItr1);
		CPPUNIT_ASSERT(tItr2>=tItr2);
	}
	
	// to itr
	{
		tItr tItr2(const_cast<TT*>(rg)+L/2,incr);

		CPPUNIT_ASSERT(tItr1==tItr1);
		CPPUNIT_ASSERT(tItr1!=tItr2);
		CPPUNIT_ASSERT(tItr1<tItr2);
		CPPUNIT_ASSERT(tItr2>tItr1);
		CPPUNIT_ASSERT(tItr1<=tItr2);
		CPPUNIT_ASSERT(tItr1<=tItr1);
		CPPUNIT_ASSERT(tItr2>=tItr1);
		CPPUNIT_ASSERT(tItr2>=tItr2);
	}

	delete[] rg;
}

template<class TT>
void test_c_tItr_all<TT>::test_dereference() {

	// lambda to generate const random range
	auto gen = [](const size_t L) -> const TT* {
		TT* res = new TT [L];
		rnd(res,L);
		return res;
	};

	const size_t L = genRndST();
	const TT* rg = gen(L);
	const size_t incr = genRndST();

	c_tItr tItr1(rg,incr);

	CPPUNIT_ASSERT_EQUAL(*rg,*tItr1);
	CPPUNIT_ASSERT_EQUAL(rg,&(*tItr1));

	const ptrdiff_t l = genRndST(0,L-1);
	CPPUNIT_ASSERT_EQUAL(rg[l],tItr1[l]);

	delete[] rg;
}

template<class TT>
void test_c_tItr_all<TT>::test_arithmetic_difference() {

	// lambda to generate const random range
	auto gen = [](const size_t L) -> const TT* {
		TT* res = new TT [L];
		rnd(res,L);
		return res;
	};

	const size_t L = genRndST();
	const TT* rg = gen(L);
	const size_t incr = genRndST();

	c_tItr tItr1(rg,incr);
	
	// increment, decrement
	{
		c_tItr tItr2(rg,incr);

		auto tItr3 = tItr2++;
		CPPUNIT_ASSERT_EQUAL(tItr1.data(),tItr3.data());
		CPPUNIT_ASSERT_EQUAL(tItr2.incr(),tItr3.incr());
		CPPUNIT_ASSERT(tItr2>tItr1);
		CPPUNIT_ASSERT_EQUAL(rg+incr,tItr2.data());
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(1),tItr2-tItr3);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(-1),tItr3-tItr2);

		auto tItr4 = ++tItr2;
		CPPUNIT_ASSERT_EQUAL(tItr2.data(),tItr4.data());
		CPPUNIT_ASSERT_EQUAL(tItr2.incr(),tItr4.incr());
		CPPUNIT_ASSERT(tItr2==tItr4);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr2-tItr4);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr4-tItr2);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(2),tItr2-tItr3);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(-2),tItr3-tItr2);

		auto tItr5 = tItr2--;
		CPPUNIT_ASSERT_EQUAL(tItr4.data(),tItr5.data());
		CPPUNIT_ASSERT_EQUAL(tItr4.incr(),tItr5.incr());
		CPPUNIT_ASSERT(tItr5>tItr2);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(1),tItr5-tItr2);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(-1),tItr2-tItr5);

		auto tItr6 = --tItr2;
		CPPUNIT_ASSERT_EQUAL(tItr2.data(),tItr6.data());
		CPPUNIT_ASSERT_EQUAL(tItr2.incr(),tItr6.incr());
		CPPUNIT_ASSERT(tItr2==tItr6);
		CPPUNIT_ASSERT(tItr1==tItr2);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr2-tItr6);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr6-tItr2);
	}

	// += -= + -
	{
		c_tItr tItr2(rg,incr);

		const ptrdiff_t l = genRndST(0,L-1);

		auto tItr3 = tItr2+l;
		CPPUNIT_ASSERT_EQUAL(tItr2.data()+l*tItr2.incr(),tItr3.data());
		CPPUNIT_ASSERT_EQUAL(tItr2.incr(),tItr3.incr());
		CPPUNIT_ASSERT_EQUAL(l,tItr3-tItr2);
		CPPUNIT_ASSERT_EQUAL(-l,tItr2-tItr3);

		auto tItr4 = tItr3-l;
		CPPUNIT_ASSERT_EQUAL(tItr2.data(),tItr4.data());
		CPPUNIT_ASSERT_EQUAL(tItr2.incr(),tItr4.incr());
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr4-tItr2);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr2-tItr4);

		tItr2+=l;
		CPPUNIT_ASSERT_EQUAL(tItr3.data(),tItr2.data());
		CPPUNIT_ASSERT_EQUAL(tItr3.incr(),tItr2.incr());
		CPPUNIT_ASSERT_EQUAL(l,tItr2-tItr1);
		CPPUNIT_ASSERT_EQUAL(-l,tItr1-tItr2);
		
		tItr2-=l;
		CPPUNIT_ASSERT_EQUAL(tItr1.data(),tItr2.data());
		CPPUNIT_ASSERT_EQUAL(tItr2.incr(),tItr2.incr());
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr1-tItr2);
		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),tItr2-tItr1);
	}

	delete[] rg;
}


template<class TT>
CppUnit::Test* test_c_tItr_all<TT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_c_tItr_all>(
		"test_ctor_assign", &test_c_tItr_all<TT>::test_ctor_assign));
	suite->addTest(new CppUnit::TestCaller<test_c_tItr_all>(
		"test_swap", &test_c_tItr_all<TT>::test_swap));
	suite->addTest(new CppUnit::TestCaller<test_c_tItr_all>(
		"test_comparison", &test_c_tItr_all<TT>::test_comparison));
	suite->addTest(new CppUnit::TestCaller<test_c_tItr_all>(
		"test_dereference", &test_c_tItr_all<TT>::test_dereference));
	suite->addTest(new CppUnit::TestCaller<test_c_tItr_all>(
		"test_arithmetic_difference", &test_c_tItr_all<TT>::test_arithmetic_difference));
	
	return suite;
}
