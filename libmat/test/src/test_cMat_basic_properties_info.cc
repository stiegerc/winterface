// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tMat_basic_properties_info.h"
#include "test_tMat_basic_properties_info.cc"


// test id
template<>
const char* test_tMat_basic_properties_info<CPX__,RE__,CPX__>::test_id() noexcept {
	return "test_cMat_basic_properties_info";
}

// instantiation
template class test_tMat_basic_properties_info<CPX__,RE__,CPX__>;
