// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TMAT_BASIC_PROPERTIES_INFO_
#define _TEST_TMAT_BASIC_PROPERTIES_INFO_


#include "lm_tMat.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class test_tMat_basic_properties_info: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;

	// tests
	void test_row_col_cpx();
	void test_M_N_L_size_lcap_ccap_empty_incr();
	void test_square_hermitian_diag();
	void test_ob_onb();

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_TMAT_BASIC_PROPERTIES_INFO_
