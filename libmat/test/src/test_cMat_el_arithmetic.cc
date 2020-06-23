// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tMat_el_arithmetic.h"
#include "test_tMat_el_arithmetic.cc"


// test id
template<>
const char* test_tMat_el_arithmetic<CPX__,RE__,CPX__>::test_id() noexcept {
	return "test_cMat_el_arithmetic";
}

// instantiation
template class test_tMat_el_arithmetic<CPX__,RE__,CPX__>;
