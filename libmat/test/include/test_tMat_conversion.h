// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TMAT_CONVERSION_
#define _TEST_TMAT_CONVERSION_

#include "lm_tMat.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class test_tMat_conversion: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;
	
	// tests
	void test_copy_fcopy_ccopy();
	void test_get_getl();
	void test_rGet_cGet();
	void test_rWOGet_cWOGet();
	void test_lower();
	void test_upper();

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_TMAT_CONVERSION_
