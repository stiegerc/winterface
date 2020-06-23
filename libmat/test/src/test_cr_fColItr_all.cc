// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_cr_tVecItr_all.h"
#include "test_cr_tVecItr_all.cc"


template<>
void test_cr_tVecItr_all<RE__,RE__,CPX__,lm_tCol<RE__,RE__,CPX__>>::test_dereference() {

	const size_t M = genRndST();
	const size_t N = genRndST();
	
	const auto tMat1 = rnd<fMat>(M,N);

	cr_tVecItr tItr1(&tMat1,N-1);
	for (size_t n=0; n!=tMat1.N(); ++n)
		CPPUNIT_ASSERT(tMat1.cAt(N-1-n)==tItr1[n]);
	
	for (size_t n=0; n!=tMat1.N(); ++n,++tItr1)
		CPPUNIT_ASSERT(tMat1.cAt(N-1-n)==*tItr1);
}

// test id
template<>
const char* test_cr_tVecItr_all<RE__,RE__,CPX__,lm_tCol<RE__,RE__,CPX__>>::test_id() noexcept {
	return "test_cr_fColItr_all";
}

// instantiation
template class test_cr_tVecItr_all<RE__,RE__,CPX__,lm_tCol<RE__,RE__,CPX__>>;
