// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_cpxItr_all.h"
#include "testTools.h"


using namespace lm__::test;

template<ptrdiff_t s>
void test_cpxItr_all<s>::test_c_ctor_assign() {
	
	// lambda to generate const random range
	auto gen = [](const size_t L) -> const CPX__* {
		CPX__* res = new CPX__ [L];
		rnd(res,L);
		return res;
	};

	const size_t L = genRndST();
	const CPX__* rg = gen(L);
	const size_t incr = genRndST();

	// default
	c_cpxItr tItr1;
	CPPUNIT_ASSERT(!tItr1.data());
	CPPUNIT_ASSERT_EQUAL(size_t(2),tItr1.incr());

	// one argument
	c_cpxItr tItr2(rg);
	CPPUNIT_ASSERT_EQUAL(reinterpret_cast<const RE__*>(rg)+s,tItr2.data());
	CPPUNIT_ASSERT_EQUAL(size_t(2),tItr2.incr());
	
	// two argument
	c_cpxItr tItr3(rg,incr);
	CPPUNIT_ASSERT_EQUAL(reinterpret_cast<const RE__*>(rg)+s,tItr3.data());
	CPPUNIT_ASSERT_EQUAL(2*incr,tItr3.incr());

	// copy
	c_cpxItr tItr4(tItr3);
	CPPUNIT_ASSERT_EQUAL(tItr3.data(),tItr4.data());
	CPPUNIT_ASSERT_EQUAL(tItr3.incr(),tItr4.incr());

	// assignment
	tItr2 = tItr3;
	CPPUNIT_ASSERT_EQUAL(tItr3.data(),tItr2.data());
	CPPUNIT_ASSERT_EQUAL(tItr3.incr(),tItr2.incr());

	// from cpxItr
	{
		cpxItr tItr1(const_cast<CPX__*>(rg),incr);
		c_cpxItr tItr2(tItr1);
		CPPUNIT_ASSERT_EQUAL(const_cast<const RE__*>(tItr1.data()),tItr2.data());
		CPPUNIT_ASSERT_EQUAL(tItr1.incr(),tItr2.incr());
	}

	delete[] rg;
}

template<ptrdiff_t s>
void test_cpxItr_all<s>::test_swap_distance() {
	
	// lambda to generate random range
	auto gen = [](const size_t L) -> CPX__* {
		CPX__* res = new CPX__ [L];
		rnd(res,L);
		return res;
	};

	const size_t L = genRndST();
	CPX__* rg = gen(L);
	const size_t incr = genRndST();

	cpxItr tItr1(rg,incr);
	cpxItr tItr2(rg+genRndST(0,L-1),incr);
	
	CPPUNIT_ASSERT_EQUAL(std::abs(tItr2-tItr1),distance(tItr1,tItr2));

	auto tItr3=tItr1;
	auto tItr4=tItr2;
	swap(tItr1,tItr2);
	CPPUNIT_ASSERT_EQUAL(tItr4,tItr1);
	CPPUNIT_ASSERT_EQUAL(tItr3,tItr2);

	delete[] rg;
}
	
template<ptrdiff_t s>
void test_cpxItr_all<s>::test_ctor_assign() {
	
	// lambda to generate random range
	auto gen = [](const size_t L) -> CPX__* {
		CPX__* res = new CPX__ [L];
		rnd(res,L);
		return res;
	};

	const size_t L = genRndST();
	CPX__* rg = gen(L);
	const size_t incr = genRndST();

	// default
	cpxItr tItr1;
	CPPUNIT_ASSERT(!tItr1.data());
	CPPUNIT_ASSERT_EQUAL(size_t(2),tItr1.incr());

	// one argument
	cpxItr tItr2(rg);
	CPPUNIT_ASSERT_EQUAL(reinterpret_cast<RE__*>(rg)+s,tItr2.data());
	CPPUNIT_ASSERT_EQUAL(size_t(2),tItr2.incr());
	
	// two argument
	cpxItr tItr3(rg,incr);
	CPPUNIT_ASSERT_EQUAL(reinterpret_cast<RE__*>(rg)+s,tItr3.data());
	CPPUNIT_ASSERT_EQUAL(2*incr,tItr3.incr());

	// copy
	cpxItr tItr4(tItr3);
	CPPUNIT_ASSERT_EQUAL(tItr3.data(),tItr4.data());
	CPPUNIT_ASSERT_EQUAL(tItr3.incr(),tItr4.incr());

	// assignment
	tItr2 = tItr3;
	CPPUNIT_ASSERT_EQUAL(tItr3.data(),tItr2.data());
	CPPUNIT_ASSERT_EQUAL(tItr3.incr(),tItr2.incr());

	delete[] rg;
}


template<ptrdiff_t s>
CppUnit::Test* test_cpxItr_all<s>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_cpxItr_all>(
		"test_c_ctor_assign", &test_cpxItr_all<s>::test_c_ctor_assign));
	suite->addTest(new CppUnit::TestCaller<test_cpxItr_all>(
		"test_ctor_assign", &test_cpxItr_all<s>::test_ctor_assign));
	suite->addTest(new CppUnit::TestCaller<test_cpxItr_all>(
		"test_swap_distance", &test_cpxItr_all<s>::test_swap_distance));
	suite->addTest(new CppUnit::TestCaller<test_cpxItr_all>(
		"test_all_dereference", &test_cpxItr_all<s>::test_all_dereference));
	
	return suite;
}
