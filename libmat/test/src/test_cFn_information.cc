// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tFn_information.h"
#include "test_tFn_information.cc"


// test id
template<>
const char* test_tFn_information<CPX__,RE__,CPX__>::test_id() noexcept {
	return "test_cFn_information";
}

// instantiation
template class test_tFn_information<CPX__,RE__,CPX__>;
