// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_r_tVecItr_all.h"
#include "test_r_tVecItr_all.cc"


template<>
void test_r_tVecItr_all<CPX__,RE__,CPX__,lm_tRow<CPX__,RE__,CPX__>>::test_dereference() {

	const size_t M = genRndST();
	const size_t N = genRndST();
	
	auto tMat1 = rnd<cMat>(M,N);

	r_tVecItr tItr1(&tMat1,M-1);
	for (size_t m=0; m!=tMat1.M(); ++m)
		CPPUNIT_ASSERT(tMat1.rAt(M-1-m)==tItr1[m]);
	
	for (size_t m=0; m!=tMat1.M(); ++m,++tItr1)
		CPPUNIT_ASSERT(tMat1.rAt(M-1-m)==*tItr1);
	
	r_tVecItr tItr2(&tMat1,M-1);
	const size_t m_ = genRndST(0,M-1);
	
	*tItr2=tItr2[m_];
	CPPUNIT_ASSERT(*tItr2==tItr2[m_]);
}

// test id
template<>
const char* test_r_tVecItr_all<CPX__,RE__,CPX__,lm_tRow<CPX__,RE__,CPX__>>::test_id() noexcept {
	return "test_r_cRowItr_all";
}

// instantiation
template class test_r_tVecItr_all<CPX__,RE__,CPX__,lm_tRow<CPX__,RE__,CPX__>>;
