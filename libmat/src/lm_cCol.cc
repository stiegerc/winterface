// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "lm_tCol.h"
#include "lm_tCol.cc"
#include "lm_types.h"
#include <cstring>


using namespace lm__;

// assignment
template<>
cCol& cCol::operator=(const cArray& rhs) noexcept {
	assert(this->L()==rhs.L());
	memcpy(this->data(),rhs.data(),M()*sizeof(CPX__));
	return *this;
}


// instantiation
template class lm_tCol<CPX__,RE__,CPX__>;
