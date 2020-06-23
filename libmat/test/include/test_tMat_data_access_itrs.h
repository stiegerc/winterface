// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TMAT_DATA_ACCESS_ITRS_
#define _TEST_TMAT_DATA_ACCESS_ITRS_


#include "lm_tMat.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class test_tMat_data_access_itrs: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;

	// tests
	void test_data_access();
	void test_el_itr();
	void test_diag_itr();
	void test_row_col_access();
	void test_row_itr();
	void test_col_itr();

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_TMAT_DATA_ACCESS_
