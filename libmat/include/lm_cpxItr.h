// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup libmat
 * @{
 */

#ifndef _LM_CPXITR_
#define _LM_CPXITR_

#include "lm_defs.h"
#include "lm_tItr.h"
#include <cassert>

template<ptrdiff_t s>
class lm_cpxItr;


/**
 * Custom increment const_iterator over real or imag part of complex range. \n
 * Template arguments are: \n
 * -s: offset, 0 for real parts, 1 for complex parts
 */
template<ptrdiff_t s=0>
class lm_c_cpxItr: public lm_c_tItr_b<RE__,lm_c_cpxItr<s>,lm_cpxItr<s>> {
public:
	/** @name types
	 */
	typedef lm_c_cpxItr<s> c_cpxItr;			//!< re/im only const_iterator
	typedef lm_cpxItr<s> cpxItr;				//!< re/im only iterator
	typedef lm_c_tItr_b<RE__,c_cpxItr,cpxItr> c_fItr_b;	//!< real const_iterator
	typedef lm_c_tItr<CPX__> c_cItr;			//!< complex const_iterator
	typedef lm_tItr<CPX__> cItr;				//!< complex iterator


	/** @name constructors
	 */
	//! default constructor
	inline explicit lm_c_cpxItr() noexcept: c_fItr_b(nullptr,2) {}
	/** constructor from const pointer and increment
	 * @param data const pointer to complex range
	 * @param incr custom increment
	 */
	inline explicit lm_c_cpxItr(const CPX__* data, const size_t incr=1) noexcept:
		c_fItr_b(reinterpret_cast<const RE__*>(data)+s,2*incr) {}
	//! constructor from re/im only const_iterator
	inline lm_c_cpxItr(const c_cpxItr& inp) noexcept:
		c_fItr_b(inp.data(),inp.incr()) {}
	//! constructor from re/im iterator
	inline lm_c_cpxItr(const cpxItr& inp) noexcept:
		c_fItr_b(inp.data(),inp.incr()) {}
	//! constructor from complex const_iterator over complex range
	inline lm_c_cpxItr(const c_cItr& inp) noexcept:
		c_cpxItr(inp.data(),inp.incr()) {}
	//! constructor from complex iterator over complex range
	inline lm_c_cpxItr(const cItr& inp) noexcept:
		c_cpxItr(inp.data(),inp.incr()) {}


	/** @name assigment
	 */
	//! swap function
	inline friend void swap(c_cpxItr& lhs, c_cpxItr& rhs) noexcept {
		using std::swap;
		swap(lhs.data_,rhs.data_);
		swap(lhs.incr_,rhs.incr_);
	}
	//! assignment operator from re/im only const_iterator
	inline c_cpxItr& operator=(const c_cpxItr& rhs) noexcept {
		this->data_ = rhs.data_;
		this->incr_ = rhs.incr_;
		return *this;
	}


	/** @name distance
	 */
	//! distance between two re/im only const_iterators
	inline friend ptrdiff_t distance(const c_cpxItr& lhs, c_cpxItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
	//! distance between re/im only iterator and re/im only const_iterator
	inline friend ptrdiff_t distance(const c_cpxItr& lhs, cpxItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
};


/**
 * Custom increment iterator over real or imag part of complex range. \n
 * Template arguments are: \n
 * -s: offset, 0 for real parts, 1 for complex parts
 */
template<ptrdiff_t s=0>
class lm_cpxItr: public lm_tItr_b<RE__,lm_cpxItr<s>,lm_c_cpxItr<s>> {
public:
	/** @name types
	 */
	typedef lm_c_cpxItr<s> c_cpxItr;			//!< re/im only const_iterator
	typedef lm_cpxItr<s> cpxItr;				//!< re/im only iterator
	typedef lm_tItr_b<RE__,cpxItr,c_cpxItr> fItr_b;		//!< real iterator
	typedef lm_tItr<CPX__> cItr;				//!< complex iterator


	/** @name constructors
	 */
	//! default constructor
	inline explicit lm_cpxItr() noexcept: fItr_b(nullptr,2) {}
	/** constructor from pointer and increment
	 * @param data pointer to complex range
	 * @param incr custom increment
	 */
	inline explicit lm_cpxItr(CPX__* data, const size_t incr=1) noexcept:
		fItr_b(reinterpret_cast<RE__*>(data)+s,2*incr) {}
	//! constructor from re/im only iterator
	inline lm_cpxItr(const cpxItr& inp) noexcept:
		fItr_b(inp.data(),inp.incr()) {}
	//! constructor from complex iterator over complex range
	inline lm_cpxItr(const cItr& inp) noexcept:
		cpxItr(inp.data(),inp.incr()) {}


	/** @name assignment
	 */
	//! swap functions
	inline friend void swap(cpxItr& lhs, cpxItr& rhs) noexcept {
		using std::swap;
		swap(lhs.data_,rhs.data_);
		swap(lhs.incr_,rhs.incr_);
	}
	//! assignment operator from re/im only iterator
	inline cpxItr& operator=(const cpxItr& rhs) noexcept {
		this->data_ = rhs.data_;
		this->incr_ = rhs.incr_;
		return *this;
	}
	

	/** @name distance
	 */
	//! distance between two re/im only iterators
	inline friend ptrdiff_t distance(const cpxItr& lhs, cpxItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
	//! distance between re/im only iterator and re/im only const_iterator
	inline friend ptrdiff_t distance(const cpxItr& lhs, c_cpxItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
};

#endif // _LM_CPXITR_

/** @}
 */
