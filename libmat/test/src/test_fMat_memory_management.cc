// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_memory_management.h"
#include "test_tMat_memory_management.cc"

// test id
template<>
const char* test_tMat_memory_management<RE__,RE__,CPX__>::test_id() noexcept {
	return "test_fMat_memory_management";
}

// instantiation
template class test_tMat_memory_management<RE__,RE__,CPX__>;
