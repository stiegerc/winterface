// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "lm_tCol.h"
#include "lm_tCol.cc"
#include "lm_types.h"
#include <cstring>


using namespace lm__;

// assignment
template<>
fCol& fCol::operator=(const fArray& rhs) noexcept {
	assert(this->L()==rhs.L());
	memcpy(this->data(),rhs.data(),M()*sizeof(RE__));
	return *this;
}


// instantiation
template class lm_tCol<RE__,RE__,CPX__>;
