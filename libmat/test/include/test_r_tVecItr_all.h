// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_R_TVECITR_ALL_
#define _TEST_R_TVECITR_ALL_

#include "lm_tVecItr.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class lm_tMat;

template<class TT, class FT, class CT, class VT>
class test_r_tVecItr_all: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;
	typedef lm_cr_tVecItr<tMat,VT> cr_tVecItr;
	typedef lm_r_tVecItr<tMat,VT> r_tVecItr;
	typedef lm_c_tVecItr<tMat,VT> c_tVecItr;
	typedef lm_tVecItr<tMat,VT> tVecItr;

	// tests
	void test_ctor_assign();
	void test_swap();
	void test_comparison();
	void test_dereference();
	void test_arithmetic_difference();

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_R_TVECITR_ALL_
