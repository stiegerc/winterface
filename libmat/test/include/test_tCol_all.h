// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TCOL_ALL_
#define _TEST_TCOL_ALL_

#include "lm_tCol.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class lm_tMat;

template<class TT, class FT, class CT>
class test_tCol_all: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;

	// tests
	void test_assign();
	void test_swap();
	void test_data_access();
	void test_information();
	void test_conversion();
	void test_matrix_arithmetic();
	void test_copy_algo_();

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

#endif // _TEST_TCOL_ALL_
