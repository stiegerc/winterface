// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tFn_functionals.h"
#include "test_tFn_functionals.cc"


// test id
template<>
const char* test_tFn_functionals<RE__,RE__,CPX__>::test_id() noexcept {
	return "test_fFn_functionals";
}

// instantiation
template class test_tFn_functionals<RE__,RE__,CPX__>;
