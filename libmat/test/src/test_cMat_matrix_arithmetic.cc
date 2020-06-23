// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tMat_matrix_arithmetic.h"
#include "test_tMat_matrix_arithmetic.cc"


// test id
template<>
const char* test_tMat_matrix_arithmetic<CPX__,RE__,CPX__>::test_id() noexcept {
	return "test_cMat_matrix_arithmetic";
}

// instantiation
template class test_tMat_matrix_arithmetic<CPX__,RE__,CPX__>;
