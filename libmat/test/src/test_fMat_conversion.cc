// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tMat_conversion.h"
#include "test_tMat_conversion.cc"


// test id
template<>
const char* test_tMat_conversion<RE__,RE__,CPX__>::test_id() noexcept {
	return "test_fMat_conversion";
}

// instantiation
template class test_tMat_conversion<RE__,RE__,CPX__>;
