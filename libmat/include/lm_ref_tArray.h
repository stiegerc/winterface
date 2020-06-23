// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup libmat
 * @{
 */

#ifndef _LM_REF_ARRAY_
#define _LM_REF_ARRAY_

#include "lm_tArray.h"
#include "lm_ops.h"
#include <string>
#include <ostream>
#include <cassert>
#include <algorithm>


/**
 * This class holds self-referential parts that cannot be included in lm_tArray.
 * lm_tMat, lm_tRow, lm_tCol all inherit directly from this with the return type
 * as itself. \n
 * Template arguments are: \n
 * - TT: the underlying type of this, either real (if this is an fMat) or complex (if this is a cMat) \n
 * - FT: the real type \n
 * - CT: the complex type \n
 * - RT: the return type \n
 */
template<class TT, class FT, class CT, class RT>
class lm_ref_tArray: public lm_tArray<TT,FT,CT> {
public:
	/** @name types
	 */
	typedef lm_tArray<FT,FT,CT> fArray;	//!< a real array
	typedef lm_tArray<CT,FT,CT> cArray;	//!< a complex array


	/** @name constructors
	 */
	//! destructor
	virtual ~lm_ref_tArray() noexcept {}


	/** @name logical
	 */
	//! elementwise AND of this with a real number
	inline RT& operator&=(const FT rhs) noexcept {
		std::for_each(this->begin(),this->end(),[rhs](TT& i){i=lm__::ops::and_s(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! elementwise AND of this with a complex number
	inline RT& operator&=(const CT& rhs) noexcept {
		std::for_each(this->begin(),this->end(),[&rhs](TT& i){i=lm__::ops::and_s(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! elementwise AND of this with a real array
	inline RT& operator&=(const fArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.cbegin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){i=lm__::ops::and_s(i,*j++);});
		return *static_cast<RT*>(this);
	}
	//! elementwise AND of this with a complex array
	inline RT& operator&=(const cArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.cbegin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){i=lm__::ops::and_s(i,*j++);});
		return *static_cast<RT*>(this);
	}
	//! elementwise OR of this with a real number
	inline RT& operator|=(const FT rhs) noexcept {
		std::for_each(this->begin(),this->end(),[rhs](TT& i){i=lm__::ops::or_s(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! elementwise OR of this with a complex number
	inline RT& operator|=(const CT& rhs) noexcept {
		std::for_each(this->begin(),this->end(),[&rhs](TT& i){i=lm__::ops::or_s(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! elementwise OR of this with a real array
	inline RT& operator|=(const fArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.cbegin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){i=lm__::ops::or_s(i,*j++);});
		return *static_cast<RT*>(this);
	}
	//! elementwise OR of this with a complex array
	inline RT& operator|=(const cArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.cbegin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){i=lm__::ops::or_s(i,*j++);});
		return *static_cast<RT*>(this);
	}


	/** @name elementwise arithmetic
	 */
	//! this + a real number
	inline RT& operator+=(const FT rhs) noexcept {
		std::for_each(this->begin(),this->end(),[rhs](TT& i){lm__::ops::plusEq(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! this + a complex number
	inline RT& operator+=(const CT& rhs) noexcept {
		std::for_each(this->begin(),this->end(),[rhs](TT& i){lm__::ops::plusEq(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! elementwise this + a real array
	inline RT& operator+=(const fArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::plusEq(i,*j++);});
		return *static_cast<RT*>(this);
	}
	//! elementwise this + a complex array
	inline RT& operator+=(const cArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::plusEq(i,*j++);});
		return *static_cast<RT*>(this);
	}
	//! this - a real number
	inline RT& operator-=(const FT rhs) noexcept {
		std::for_each(this->begin(),this->end(),[rhs](TT& i){lm__::ops::minusEq(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! this - a complex number
	inline RT& operator-=(const CT& rhs) noexcept {
		std::for_each(this->begin(),this->end(),[rhs](TT& i){lm__::ops::minusEq(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! elementwise this - a real array
	inline RT& operator-=(const fArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::minusEq(i,*j++);});
		return *static_cast<RT*>(this);
	}
	//! elementwise this - a complex array
	inline RT& operator-=(const cArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::minusEq(i,*j++);});
		return *static_cast<RT*>(this);
	}
	//! this * a real number
	inline RT& operator*=(const FT rhs) noexcept {
		std::for_each(this->begin(),this->end(),[rhs](TT& i){lm__::ops::prodEq(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! this * a complex number
	inline RT& operator*=(const CT& rhs) noexcept {
		std::for_each(this->begin(),this->end(),[rhs](TT& i){lm__::ops::prodEq(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! elementwise this * a real array
	inline RT& operator*=(const fArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::prodEq(i,*j++);});
		return *static_cast<RT*>(this);
	}
	//! elementwise this * a complex array
	inline RT& operator*=(const cArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::prodEq(i,*j++);});
		return *static_cast<RT*>(this);
	}
	//! this / a real number
	inline RT& operator/=(const FT rhs) noexcept {
		std::for_each(this->begin(),this->end(),[rhs](TT& i){lm__::ops::divEq(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! this / a complex number
	inline RT& operator/=(const CT& rhs) noexcept {
		std::for_each(this->begin(),this->end(),[rhs](TT& i){lm__::ops::divEq(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! elementwise this / a real array
	inline RT& operator/=(const fArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::divEq(i,*j++);});
		return *static_cast<RT*>(this);
	}
	//! elementwise this / a complex array
	inline RT& operator/=(const cArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::divEq(i,*j++);});
		return *static_cast<RT*>(this);
	}
	//! this % a real number
	inline RT& operator%=(const FT rhs) noexcept {
		std::for_each(this->begin(),this->end(),[rhs](TT& i){lm__::ops::modEq(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! this % a complex number
	inline RT& operator%=(const CT& rhs) noexcept {
		std::for_each(this->begin(),this->end(),[rhs](TT& i){lm__::ops::modEq(i,rhs);});
		return *static_cast<RT*>(this);
	}
	//! elementwise this % a real array
	inline RT& operator%=(const fArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::modEq(i,*j++);});
		return *static_cast<RT*>(this);
	}
	//! elementwise this % a complex array
	inline RT& operator%=(const cArray& rhs) noexcept {
		assert(msize(*this)==msize(rhs));
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::modEq(i,*j++);});
		return *static_cast<RT*>(this);
	}
};


#endif // _LM_REF_ARRAY_

/** @}
 */
