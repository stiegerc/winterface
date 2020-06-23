// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tMat_assign.h"
#include "test_tMat_assign.cc"
#include "lm_fn.h"

using namespace lm__;

template<>
void test_tMat_assign<RE__,RE__,CPX__>::test_operator_equal_cArray() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	// from cMat, same size
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		const auto ptr = tMat1.begin();
		const auto lcap = tMat1.lcap();
		const auto ccap = tMat1.ccap();

		tMat1 = tMat2;

		CPPUNIT_ASSERT_EQUAL(lcap,tMat1.lcap());
		CPPUNIT_ASSERT_EQUAL(ccap,tMat1.ccap());
		CPPUNIT_ASSERT_EQUAL(ptr,tMat1.begin());

		CPPUNIT_ASSERT(real(tMat2)==tMat1);
	}
	
	// from cMat, larger size
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<cMat>(2*M,2*N);

		const auto ptr = tMat1.begin();

		tMat1 = tMat2;

		CPPUNIT_ASSERT_EQUAL(tMat2.lcap(),tMat1.lcap());
		CPPUNIT_ASSERT_EQUAL(tMat2.ccap(),tMat1.ccap());
		CPPUNIT_ASSERT(ptr!=tMat1.begin());

		CPPUNIT_ASSERT(real(tMat2)==tMat1);
	}

	// from cRow
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		const size_t m = genRndST(0,M-1);
		auto j = tMat2.crBegin()+m;
		
		tMat1 = *j;

		CPPUNIT_ASSERT(real(*j)==tMat1);
	}
	
	// from cCol
	{
		auto tMat1 = rnd<tMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		const size_t n = genRndST(0,N-1);
		auto j = tMat2.ccBegin()+n;
		
		tMat1 = *j;

		CPPUNIT_ASSERT(real(*j)==tMat1);
	}
}


// test id
template<>
const char* test_tMat_assign<RE__,RE__,CPX__>::test_id() noexcept {
	return "test_fMat_assign";
}

// instantiation
template class test_tMat_assign<RE__,RE__,CPX__>;
