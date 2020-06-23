// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_c_tItr_all.h"
#include "test_c_tItr_all.cc"

// test id
template<>
const char* test_c_tItr_all<CPX__>::test_id() noexcept {
	return "test_c_cItr_all";
}

// instantiation
template class test_c_tItr_all<CPX__>;
