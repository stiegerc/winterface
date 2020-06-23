// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tMat_data_access_itrs.h"
#include "test_tMat_data_access_itrs.cc"


// test id
template<>
const char* test_tMat_data_access_itrs<RE__,RE__,CPX__>::test_id() noexcept {
	return "test_fMat_data_access_itrs";
}

// instantiation
template class test_tMat_data_access_itrs<RE__,RE__,CPX__>;
