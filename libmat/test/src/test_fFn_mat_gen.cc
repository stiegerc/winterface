// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tFn_mat_gen.h"
#include "test_tFn_mat_gen.cc"


// test id
template<>
const char* test_tFn_mat_gen<RE__,RE__,CPX__>::test_id() noexcept {
	return "test_fFn_mat_gen";
}

// instantiation
template class test_tFn_mat_gen<RE__,RE__,CPX__>;
