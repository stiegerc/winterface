// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_TMAT_LOGICAL_
#define _TEST_TMAT_LOGICAL_

#include "lm_tMat.h"
#include <cppunit/extensions/HelperMacros.h>

template<class TT, class FT, class CT>
class test_tMat_logical: public CppUnit::TestFixture {
public:
	// types
	typedef lm_tMat<FT,FT,CT> fMat;
	typedef lm_tMat<CT,FT,CT> cMat;
	typedef lm_tMat<TT,FT,CT> tMat;
	
	// tests
	void test_logical();
	void test_permutation();
	void test_operator_not();
	void test_operator_and();
	void test_operator_or();

protected:
	static const char* test_id() noexcept;

public:
	static CppUnit::Test* suite();
};

#endif // _TEST_TMAT_LOGICAL_
