// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup libmat
 * @{
 */

#ifndef _LM_TITR_
#define _LM_TITR_

#include <cassert>
#include <ostream>
#include <iterator>
#include <stddef.h>
#include <iostream>
#include <utility>


template<class TT>
class lm_tItr;


/**
 * Custom increment const_iterator base class.
 * This is used to traverse a matrix in different ways, e.g. for increment 1 along columns,
 * for increment M along rows, for increment M+1 along diagonals. \n
 * Template arguments are: \n
 * - TT: underlying data type, i.e. real or complex
 * - RT: return type, i.e. the iterator inheriting from this
 * - RT_: other return type, i.e.  the const_iterator going with the iterator type (or vice versa)
 */
template<class TT, class RT, class RT_>
class lm_c_tItr_b: public std::iterator<std::random_access_iterator_tag,	// iterator_category
				        TT,					// value_type
				        ptrdiff_t,				// difference_type
				        const TT*,				// pointer
				        const TT&> {				// reference
public:
	/** @name constructors
	 */
	/** constructor from const pointer and increment
	 * @param data data pointer
	 * @param incr custom increment
	 */
	inline explicit lm_c_tItr_b(const TT* data=nullptr, const size_t incr=1) noexcept:
		data_(data), incr_(incr) { assert(incr_); }
	

	/** @name information
	 */
	//! raw data const pointer
	inline const TT* data() const noexcept { return data_; }
	//! custom increment
	inline size_t incr() const noexcept { return incr_; }	
	

	/** @name const dereference
	 */
	//! const dereference operator
	inline const TT& operator*() const noexcept { return *data(); }
	//! const linear indexing operator
	inline const TT& operator[](const ptrdiff_t i) const noexcept { return *(data()+i); }
	//! const arrow operator
	inline const TT* operator->() const noexcept { return data(); }
	

	/** @name increment, decrement
	 */
	//! post-increment by custom increment
	inline RT& operator++() noexcept { this->data_+=this->incr(); return *static_cast<RT*>(this); }
	//! pre-increment by custom increment
	inline RT operator++(int) noexcept { RT res(*static_cast<RT*>(this)); ++(*this); return res; }
	//! post-decrement by custom increment
	inline RT& operator--() noexcept { this->data_-=this->incr(); return *static_cast<RT*>(this); }
	//! post-increment by custom increment
	inline RT operator--(int) noexcept { RT res(*static_cast<RT*>(this)); --(*this); return res; }
	

	/** @name arithmetic operators
	 */
	//! operator+=
	inline RT& operator+=(const ptrdiff_t rhs) noexcept {
		this->data_+=rhs*this->incr(); return *static_cast<RT*>(this);
	}
	//! operator+
	inline RT operator+(const ptrdiff_t rhs) const noexcept { return RT(*static_cast<const RT*>(this))+=rhs; }
	//! operator-=
	inline RT& operator-=(const ptrdiff_t rhs) noexcept {
		this->data_-=rhs*this->incr(); return *static_cast<RT*>(this);
	}
	//! operator+
	inline RT operator-(const ptrdiff_t rhs) const noexcept { return RT(*static_cast<const RT*>(this))-=rhs; }
	

	/** @name arithmetic as rhs
	 */
	//! operator+
	inline friend RT operator+(const ptrdiff_t lhs, const RT& rhs) noexcept {
		return rhs+lhs;
	}


	/** @name fundamental comparison operators
	 */
	//! comparison to iterator
	inline bool operator==(const RT& rhs) const noexcept { return data()==rhs.data(); }
	//! comparison to other iterator
	inline bool operator==(const RT_& rhs) const noexcept { return data()==rhs.data(); }
	//! comparison to iterator
	inline bool operator<(const RT& rhs) const noexcept { return data()<rhs.data(); }
	//! comparison to other iterator
	inline bool operator<(const RT_& rhs) const noexcept { return data()<rhs.data(); }

	
	/** @name derived comparison operators
	 */
	//! comparison to iterator
	inline bool operator!=(const RT& rhs) const noexcept { return !operator==(rhs); }
	//! comparison to other iterator
	inline bool operator!=(const RT_& rhs) const noexcept { return !operator==(rhs); }
	//! comparison to iterator
	inline bool operator>(const RT& rhs) const noexcept { return rhs<*static_cast<const RT*>(this); }
	//! comparison to other iterator
	inline bool operator>(const RT_& rhs) const noexcept { return rhs<*static_cast<const RT*>(this); }
	//! comparison to iterator
	inline bool operator<=(const RT& rhs) const noexcept { return !operator>(rhs); }
	//! comparison to other iterator
	inline bool operator<=(const RT_& rhs) const noexcept { return !operator>(rhs); }
	//! comparison to iterator
	inline bool operator>=(const RT& rhs) const noexcept { return !operator<(rhs); }
	//! comparison to other iterator
	inline bool operator>=(const RT_& rhs) const noexcept { return !operator<(rhs); }


	/** @name difference
	 */
	//! iterator difference
	inline friend ptrdiff_t operator-(const RT& lhs, const RT& rhs) noexcept {
		assert(lhs.incr()==rhs.incr());
		return (lhs.data()-rhs.data())/ptrdiff_t(lhs.incr());
	}
	//! other iterator difference
	inline friend ptrdiff_t operator-(const RT& lhs, const RT_& rhs) noexcept {
		assert(lhs.incr()==rhs.incr());
		return (lhs.data()-rhs.data())/ptrdiff_t(lhs.incr());
	}


	/** @name streaming
	 */
	//! streaming operator
	inline friend std::ostream& operator<<(std::ostream& os, const lm_c_tItr_b& i) noexcept {
		return (os << "i(" << i.data() << "," << i.incr() << ")");
	}


protected:
	/** @name member variables
	 */
	const TT* data_;	//!< const data pointer
	size_t incr_;		//!< custom increment
};


/**
 * Custom increment const_iterator class. \n
 * Template arguments are: \n
 * - TT: the underlying data type, i.e. real or complex
 */
template<class TT>
class lm_c_tItr: public lm_c_tItr_b<TT,lm_c_tItr<TT>,lm_tItr<TT>> {
public:
	/** @name types
	 */
	typedef lm_c_tItr<TT> c_tItr;	//!< const_iterator
	typedef lm_tItr<TT> tItr;	//!< iterator


	/** @name constructors
	 */
	using lm_c_tItr_b<TT,c_tItr,tItr>::lm_c_tItr_b;
	//! default constructor
	lm_c_tItr()=default;
	//! constructor from const_iterator
	inline lm_c_tItr(const c_tItr& inp) noexcept: c_tItr(inp.data(),inp.incr()) {}
	//! constructor from iterator
	inline lm_c_tItr(const tItr& inp) noexcept: c_tItr(inp.data(),inp.incr()) {}
	

	/** @name assignment
	 */
	//! swap function
	inline friend void swap(c_tItr& lhs, c_tItr& rhs) noexcept {
		using std::swap;
		swap(lhs.data_,rhs.data_);
		swap(lhs.incr_,rhs.incr_);
	}
	//! assignment operator from const_iterator
	inline c_tItr& operator=(const c_tItr& rhs) noexcept {
		this->data_ = rhs.data_;
		this->incr_ = rhs.incr_;
		return *this;
	}

	/** @name distance
	 */
	//! distance between two const_iterators
	inline friend ptrdiff_t distance(const c_tItr& lhs, const c_tItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
	//! distance between const_iterator and iterator
	inline friend ptrdiff_t distance(const c_tItr& lhs, const tItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
};



/**
 * Custom increment iterator base class.
 * This is used to traverse a matrix in different ways, e.g. for increment 1 along columns,
 * for increment M along rows, for increment M+1 along diagonals. \n
 * Template arguments are: \n
 * - TT: underlying data type, i.e. real or complex
 * - RT: return type, i.e. the iterator inheriting from this
 * - RT_: other return type, i.e.  the const_iterator going with the iterator type (or vice versa)
 */
template<class TT, class RT, class RT_>
class lm_tItr_b: public std::iterator<std::random_access_iterator_tag,		// iterator_category
				      TT,					// value_type
				      ptrdiff_t,				// difference_type
				      TT*,					// pointer
				      TT&> {					// reference
public:
	/** @name constructors
	 */
	/** constructor from pointer and increment
	 * @param data data pointer
	 * @param incr custom increment
	 */
	inline explicit lm_tItr_b(TT* data=nullptr, const size_t incr=1) noexcept:
		data_(data), incr_(incr) { assert(incr_); }
	
	/** @name information
	 */
	//! raw data pointer
	inline TT* data() const noexcept { return data_; }
	//! custom increment
	inline size_t incr() const noexcept { return incr_; }	
	

	/** @name dereference
	 */
	//! dereference operator
	inline TT& operator*() const noexcept { return *data(); }
	//! linear indexing operator
	inline TT& operator[](const ptrdiff_t i) const noexcept { return *(data()+i); }
	//! arrow operator
	inline TT* operator->() const noexcept { return data(); }
	

	/** @name increment, decrement
	 */
	//! post-increment by custom increment
	inline RT& operator++() noexcept { this->data_+=this->incr(); return *static_cast<RT*>(this); }
	//! pre-increment by custom increment
	inline RT operator++(int) noexcept { RT res(*static_cast<RT*>(this)); ++(*this); return res; }
	//! post-decrement by custom increment
	inline RT& operator--() noexcept { this->data_-=this->incr(); return *static_cast<RT*>(this); }
	//! pre-decrement by custom increment
	inline RT operator--(int) noexcept { RT res(*static_cast<RT*>(this)); --(*this); return res; }


	/** @name arithmetic operators
	 */
	//! operator+=
	inline RT& operator+=(const ptrdiff_t rhs) noexcept {
		this->data_+=rhs*this->incr(); return *static_cast<RT*>(this);
	}
	//! operator+
	inline RT operator+(const ptrdiff_t rhs) const noexcept { return RT(*static_cast<const RT*>(this))+=rhs; }
	//! operator-=
	inline RT& operator-=(const ptrdiff_t rhs) noexcept {
		this->data_-=rhs*this->incr(); return *static_cast<RT*>(this);
	}
	//! operator-
	inline RT operator-(const ptrdiff_t rhs) const noexcept { return RT(*static_cast<const RT*>(this))-=rhs; }


	/** @name arithmetic as rhs
	 */
	//! opperator+
	inline friend RT operator+(const ptrdiff_t lhs, const RT& rhs) noexcept {
		return rhs+lhs;
	}


	/** @name fundamental comparison operators
	 */
	//! comparison to iterator
	inline bool operator==(const RT& rhs) const noexcept { return data()==rhs.data(); }
	//! comparison to const_iterator
	inline bool operator==(const RT_& rhs) const noexcept { return data()==rhs.data(); }
	//! comparison to iterator
	inline bool operator<(const RT& rhs) const noexcept { return data()<rhs.data(); }
	//! comparison to const_iterator
	inline bool operator<(const RT_& rhs) const noexcept { return data()<rhs.data(); }


	/** derived comparison operators
	 */
	//! comparison to iterator
	inline bool operator!=(const RT& rhs) const noexcept { return !operator==(rhs); }
	//! comparison to const_iterator
	inline bool operator!=(const RT_& rhs) const noexcept { return !operator==(rhs); }
	//! comparison to iterator
	inline bool operator>(const RT& rhs) const noexcept { return rhs<*static_cast<const RT*>(this); }
	//! comparison to const_iterator
	inline bool operator>(const RT_& rhs) const noexcept { return rhs<*static_cast<const RT*>(this); }
	//! comparison to iterator
	inline bool operator<=(const RT& rhs) const noexcept { return !operator>(rhs); }
	//! comparison to const_iterator
	inline bool operator<=(const RT_& rhs) const noexcept { return !operator>(rhs); }
	//! comparison to iterator
	inline bool operator>=(const RT& rhs) const noexcept { return !operator<(rhs); }
	//! comparison to const_iterator
	inline bool operator>=(const RT_& rhs) const noexcept { return !operator<(rhs); }


	/** @name difference
	 */
	//! iterator difference
	inline friend ptrdiff_t operator-(const RT& lhs, const RT& rhs) noexcept {
		assert(lhs.incr()==rhs.incr());
		return (lhs.data()-rhs.data())/ptrdiff_t(lhs.incr());
	}
	//! other iterator difference
	inline friend ptrdiff_t operator-(const RT& lhs, const RT_& rhs) noexcept {
		assert(lhs.incr()==rhs.incr());
		return (lhs.data()-rhs.data())/ptrdiff_t(lhs.incr());
	}


	/** @name streaming
	 */
	//! streaming operator
	inline friend std::ostream& operator<<(std::ostream& os, const lm_tItr_b& i) noexcept {
		return (os << "i(" << i.data() << "," << i.incr() << ")");
	}


protected:
	/** @name member variables
	 */
	TT* data_;		//!< data pointer
	size_t incr_;		//!< custom increment
};


/**
 * Custom increment iterator class. \n
 * Template arguments are: \n
 * - TT: the underlying data type, i.e. real or complex
 */
template<class TT>
class lm_tItr: public lm_tItr_b<TT,lm_tItr<TT>,lm_c_tItr<TT>> {
public:
	/** @name types
	 */
	typedef ::lm_tItr<TT> tItr;		//!< iterator
	typedef ::lm_c_tItr<TT> c_tItr;		//!< const_iterator


	/** @name constructors
	 */
	using lm_tItr_b<TT,tItr,c_tItr>::lm_tItr_b;
	//! default constructor
	lm_tItr()=default;
	//! constructor from iterator
	inline lm_tItr(const tItr& inp) noexcept: tItr(inp.data(),inp.incr()) {}
	

	/** @name assignment
	 */
	//! swap function
	inline friend void swap(tItr& lhs, tItr& rhs) noexcept {
		using std::swap;
		swap(lhs.data_,rhs.data_);
		swap(lhs.incr_,rhs.incr_);
	}
	//! assignment operator from iterator
	inline tItr& operator=(const tItr& rhs) noexcept {
		this->data_ = rhs.data_;
		this->incr_ = rhs.incr_;
		return *this;
	}
	

	/** @name distance
	 */
	//! distance between two iterators
	inline friend ptrdiff_t distance(const tItr& lhs, const tItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
	//! distance between iterator and const_iterator
	inline friend ptrdiff_t distance(const tItr& lhs, const c_tItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
};

#endif // _LM_TITR_

/** @}
 */
