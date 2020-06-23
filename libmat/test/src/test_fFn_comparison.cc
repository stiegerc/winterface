// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tFn_comparison.h"
#include "test_tFn_comparison.cc"


// test id
template<>
const char* test_tFn_comparison<RE__,RE__,CPX__>::test_id() noexcept {
	return "test_fFn_comparison";
}

// instantiation
template class test_tFn_comparison<RE__,RE__,CPX__>;
