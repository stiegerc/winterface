// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TMAT_EL_ARITHMETIC_
#define _TEST_TMAT_EL_ARITHMETIC_

#include "lm_tMat.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class test_tMat_el_arithmetic: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;
	
	// tests
	void test_operator_minus();
	void test_operator_plus_pluseq();
	void test_operator_minus_minuseq();
	void test_operator_prod_prodeq();
	void test_operator_div_diveq();
	void test_operator_plus_tlhs();
	void test_operator_minus_tlhs();
	void test_operator_prod_tlhs();
	void test_operator_div_tlhs();

	// delta
#ifdef DOUBLE__
	static constexpr FT delta = 1e-10;
#endif
#ifdef SINGLE__
	static constexpr FT delta = 1e-3f;
#endif

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_TMAT_EL_ARITHMETIC_
