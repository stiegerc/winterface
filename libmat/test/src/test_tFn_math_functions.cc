// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tFn_math_functions.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


template<class TT, class FT, class CT>
void test_tFn_math_functions<TT,FT,CT>::test_abs() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	const auto tMat1 = rnd<tMat>(M,N);	
	
	{
		const auto tMat2 = abs(tMat1);
		for (size_t i=0; i!=tMat2.L(); ++i)
			CPPUNIT_ASSERT_DELTA(std::abs(tMat1[i]),tMat2[i],delta);
	}
	{
		size_t i=0;
		for (auto ri=tMat1.crBegin(), re=tMat1.crEnd(); ri!=re; ++ri, ++i)
			CPPUNIT_ASSERT(abs(*ri)==abs(tMat1.rGet(i)));
	}
	{
		size_t i=0;
		for (auto ci=tMat1.ccBegin(), ce=tMat1.ccEnd(); ci!=ce; ++ci, ++i)
			CPPUNIT_ASSERT(abs(*ci)==abs(tMat1.cGet(i)));
	}
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tFn_math_functions<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tFn_math_functions>(
		"test_abs", &test_tFn_math_functions<TT,FT,CT>::test_abs));
	suite->addTest(new CppUnit::TestCaller<test_tFn_math_functions>(
		"test_ceilEq_ceil", &test_tFn_math_functions<TT,FT,CT>::test_ceilEq_ceil));
	suite->addTest(new CppUnit::TestCaller<test_tFn_math_functions>(
		"test_floorEq_floor", &test_tFn_math_functions<TT,FT,CT>::test_floorEq_floor));
	suite->addTest(new CppUnit::TestCaller<test_tFn_math_functions>(
		"test_roundEq_round", &test_tFn_math_functions<TT,FT,CT>::test_roundEq_round));
	suite->addTest(new CppUnit::TestCaller<test_tFn_math_functions>(
		"test_signEq_sign", &test_tFn_math_functions<TT,FT,CT>::test_signEq_sign));
	
	return suite;
}
