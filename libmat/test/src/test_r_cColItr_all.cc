// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_r_tVecItr_all.h"
#include "test_r_tVecItr_all.cc"


template<>
void test_r_tVecItr_all<CPX__,RE__,CPX__,lm_tCol<CPX__,RE__,CPX__>>::test_dereference() {

	const size_t M = genRndST();
	const size_t N = genRndST();
	
	auto tMat1 = rnd<cMat>(M,N);

	r_tVecItr tItr1(&tMat1,N-1);
	for (size_t n=0; n!=tMat1.N(); ++n)
		CPPUNIT_ASSERT(tMat1.cAt(N-1-n)==tItr1[n]);
	
	for (size_t n=0; n!=tMat1.N(); ++n,++tItr1)
		CPPUNIT_ASSERT(tMat1.cAt(N-1-n)==*tItr1);
	
	tVecItr tItr2(&tMat1,0);
	const size_t n_ = genRndST(0,N-1);
	
	*tItr2=tItr2[n_];
	CPPUNIT_ASSERT(*tItr2==tItr2[n_]);
}

// test id
template<>
const char* test_r_tVecItr_all<CPX__,RE__,CPX__,lm_tCol<CPX__,RE__,CPX__>>::test_id() noexcept {
	return "test_r_cColItr_all";
}

// instantiation
template class test_r_tVecItr_all<CPX__,RE__,CPX__,lm_tCol<CPX__,RE__,CPX__>>;
