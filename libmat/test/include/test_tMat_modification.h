// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TMAT_MODIFICATION_
#define _TEST_TMAT_MODIFICATION_

#include "lm_tMat.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class test_tMat_modification: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;
	
	// tests
	void test_rShift_cShift();
	void test_rInsert();
	void test_cInsert();
	void test_set_setl();
	void test_rRm();
	void test_cRm();
	void test_dRm();
	void test_inv();
	
	// delta
#ifdef DOUBLE__
	static constexpr FT delta = 1e-10;
#endif
#ifdef SINGLE__
	static constexpr FT delta = 1e-4f;
#endif

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_TMAT_MODIFICATION_
