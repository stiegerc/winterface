// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tFn_orth.h"
#include "test_tFn_orth.cc"


// test id
template<>
const char* test_tFn_orth<RE__,RE__,CPX__>::test_id() noexcept {
	return "test_fFn_orth";
}

// instantiation
template class test_tFn_orth<RE__,RE__,CPX__>;
