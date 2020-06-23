// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tItr_all.h"
#include "test_tItr_all.cc"

// test id
template<>
const char* test_tItr_all<RE__>::test_id() noexcept {
	return "test_fItr_all";
}

// instantiation
template class test_tItr_all<RE__>;
