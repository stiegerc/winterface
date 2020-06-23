// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TMAT_CTOR_
#define _TEST_TMAT_CTOR_


#include "lm_tMat.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class test_tMat_ctor: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;

	// tests
	void test_mn();
	void test_m();
	void test_default();
	void test_prealloc();
	void test_from_itrs();
	void test_initl_row_col();
	void test_initl();
	void test_initl_m();
	void test_initl_mn();
	void test_file();
	void test_fstream();
	void test_vec_tt();
	void test_vec_size_t();
	void test_re_im();
	void test_copy_fArray();
	void test_copy_cArray();
	void test_copy();
	void test_move();

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_TMAT_CTOR_
