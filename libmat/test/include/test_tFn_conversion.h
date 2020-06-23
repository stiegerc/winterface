// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TFN_CONVERSION_
#define _TEST_TFN_CONVERSION_

#include "lm_tMat.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class test_tFn_conversion: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;
	
	// tests
	void test_msum();
	void test_nsum();
	void test_mprod();
	void test_nprod();
	void test_mmin();
	void test_nmin();
	void test_mmax();
	void test_nmax();
	void test_mmean();
	void test_nmean();
	void test_mnormsq();
	void test_nnormsq();
	void test_mnorm();
	void test_nnorm();
	void test_mdot();
	void test_ndot();
	void test_mdotu();
	void test_ndotu();
	void test_inv();
	void test_diag();
	void test_R_C();
	void test_T();
	void test_Iz_Inz();
	void test_conj();
	void test_real();
	void test_imag();

	// delta
#ifdef DOUBLE__
	static constexpr RE__ delta = 1e-10;
#endif
#ifdef SINGLE__
	static constexpr RE__ delta = 1e-3f;
#endif

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_TFN_CONVERSION_
