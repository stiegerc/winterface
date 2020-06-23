// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TMAT_MATRIX_ARITHMETIC_
#define _TEST_TMAT_MATRIX_ARITHMETIC_

#include "lm_tMat.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class test_tMat_matrix_arithmetic: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;
	
	// tests
	void test_leftDivide();
	void test_leftDivideEq();
	void test_prod();

	// delta
#ifdef DOUBLE__
	static constexpr FT delta = 1e-6;
#endif
#ifdef SINGLE__
	static constexpr FT delta = 1e-2f;
#endif

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_TMAT_MATRIX_ARITHMETIC_
