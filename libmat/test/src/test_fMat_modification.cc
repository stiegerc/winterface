// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tMat_modification.h"
#include "test_tMat_modification.cc"


// test id
template<>
const char* test_tMat_modification<RE__,RE__,CPX__>::test_id() noexcept {
	return "test_fMat_modification";
}

// instantiation
template class test_tMat_modification<RE__,RE__,CPX__>;
