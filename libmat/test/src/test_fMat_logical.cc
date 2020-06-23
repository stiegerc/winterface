// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tMat_logical.h"
#include "test_tMat_logical.cc"


// test id
template<>
const char* test_tMat_logical<RE__,RE__,CPX__>::test_id() noexcept {
	return "test_fMat_logical";
}

// instantiation
template class test_tMat_logical<RE__,RE__,CPX__>;
