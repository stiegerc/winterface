// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TFN_COMPARISON_
#define _TEST_TFN_COMPARISON_

#include "lm_tMat.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class test_tFn_comparison: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;
	
	// tests
	void test_all();
	void test_mall();
	void test_nall();
	void test_any();
	void test_many();
	void test_nany();
	void test_find();
	void test_nnz();

	// delta
#ifdef DOUBLE__
	static constexpr RE__ delta = 1e-10;
#endif
#ifdef SINGLE__
	static constexpr RE__ delta = 1e-5f;
#endif

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_TFN_COMPARISON_
