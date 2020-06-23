// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tMat_assign.h"
#include "test_tMat_assign.cc"


// test id
template<>
const char* test_tMat_assign<CPX__,RE__,CPX__>::test_id() noexcept {
	return "test_cMat_assign";
}

// instantiation
template class test_tMat_assign<CPX__,RE__,CPX__>;
