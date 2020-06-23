// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TFN_MATH_FUNCTIONS_
#define _TEST_TFN_MATH_FUNCTIONS_

#include "lm_tMat.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class test_tFn_math_functions: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;
	
	// tests
	void test_abs();
	void test_ceilEq_ceil();
	void test_floorEq_floor();
	void test_roundEq_round();
	void test_signEq_sign();

	// delta
#ifdef DOUBLE__
	static constexpr FT delta = 1e-10;
#endif
#ifdef SINGLE__
	static constexpr FT delta = 1e-5f;
#endif

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_TFN_MATH_FUNCTIONS_
