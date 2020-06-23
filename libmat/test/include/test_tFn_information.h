// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TFN_INFORMATION_
#define _TEST_TFN_INFORMATION_

#include "lm_tMat.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class test_tFn_information: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;
	
	// tests
	void test_isnan();
	void test_length();
	void test_rank();

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

#endif // _TEST_TFN_INFORMATION_
