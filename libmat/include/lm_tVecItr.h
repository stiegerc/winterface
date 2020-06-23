// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup libmat
 * @{
 */

#ifndef _LM_TVECITR_
#define _LM_TVECITR_

#include <iterator>
#include <cstddef>
#include <cassert>
#include <iostream>

template<class MT, class VT>
class lm_tVecItr;
template<class MT, class VT>
class lm_cr_tVecItr;
template<class MT, class VT>
class lm_r_tVecItr;


/**
 * Base class for vector iterators. I.e. row and column iterators.
 * These are 'stashing iterators', i.e. iterators that return a reference to a member element,
 * or in other words: they themselves contain what they are pointing to.
 * Therefore std::reverse_iterator does not work with these and manual implementations of reverse
 * iterators are provided. \n
 * \n
 * Also the following condition of forward iterators does not hold
 *    If a and b are both dereferenceable,
 *    then a == b if and only if *a and *b are bound to the same object.
 * The reason is again the reference of operator* to a member item.
 * But since these items are a tVec which have no substance in themselves,
 * but are wrappers for a specific row or column in a matrix.
 * So *a and *b are bound to different wrappers that point to
 * the same object and thus the condition can be regarded as fulfilled. \n
 * \n
 * In fact tVecs are themselves almost iterators but we need certain
 * operators such as operator++ to act on the section of the matrix they point
 * to as if they were a matrix themselves and at the same time we need
 * it to work in the way an iterator would. The classes in this file thus
 * provide the iterator interface for the matrix interfaces in lm_tCol and lm_tRow.
 * The only way to make this work is for the vector iterator classes
 * to hold the vector it points to as a member and thus we end up returning
 * references to a member item. In all other ways these can be considered standard iterators. \n
 * Template arguments are: \n
 * - MT: 'host' matrix type
 * - VT: vector type, i.e. lm_tRow or lm_tCol
 * - RT: return type, i.e. the iterator inheriting from this
 * - RT_: other return type, i.e. the const_iterator going with the iterator type (or vice versa)
 */
template<class MT, class VT, class RT, class RT_>
class lm_tVecItr_b {
public:
	/** @name constructors
	 */
	/** constructor from matrix
	 * @param ptr pointer to the host matrix
	 * @param i row or column index in the matrix
	 */
	inline explicit lm_tVecItr_b(const MT* ptr=nullptr,
		const ptrdiff_t i=0) noexcept: vec_(ptr,i) {}

	/** @name information
	 */
	//! row or column index in the 'host' matrix
	inline ptrdiff_t i() const noexcept { return vec_.i(); }
	//! the number of rows or columns in the 'host' matrix
	inline ptrdiff_t iL() const noexcept { return vec_.iL(); }
	//! const pointer to the 'host' matrix
	inline const MT* ptr() const noexcept { return vec_.ptr(); }


	/** @name casts
	 */
	//! casts to signed integer, i.e. the row or column index in the 'host' matrix
	inline explicit operator ptrdiff_t() const { return i(); }
	//! casts to integer, i.e. the row or column index in the 'host' matrix
	inline explicit operator size_t() const { return size_t(i()); }


	/** @name increment, decrement
	 */
	//! iterator post-increment
	inline RT& operator++() noexcept {
		++this->vec_.i_; return *static_cast<RT*>(this);
	}
	//! iterator pre-increment
	inline RT operator++(int) noexcept {
		RT res(*static_cast<RT*>(this)); ++(*this); return res;
	}
	//! iterator post-decrement
	inline RT& operator--() noexcept {
		--this->vec_.i_; return *static_cast<RT*>(this);
	}
	//! iterator pre-decrement
	inline RT operator--(int) noexcept {
		RT res(*static_cast<RT*>(this)); --(*this); return res;
	}


	/** @name arithmetic operators
	 */
	//! iterator +=
	inline RT& operator+=(const ptrdiff_t rhs) noexcept {
		this->vec_.i_+=rhs;
		return *static_cast<RT*>(this);
	}
	//! iterator +
	inline RT operator+(const ptrdiff_t rhs) const noexcept {
		return RT(*static_cast<const RT*>(this))+=rhs;
	}
	//! iterator -=
	inline RT& operator-=(const ptrdiff_t rhs) noexcept {
		this->vec_.i_-=rhs;
		return *static_cast<RT*>(this);
	}
	//! iterator -
	inline RT operator-(const ptrdiff_t rhs) const noexcept {
		return RT(*static_cast<const RT*>(this))-=rhs;
	}
	

	/** @name arithmetic as rhs
	 */
	//! iterator +
	inline friend RT operator+(const ptrdiff_t lhs, const RT& rhs) noexcept {
		return rhs+lhs;
	}


	/** fundamental comparison operators
	 */
	//! comparison to iterator
	inline bool operator==(const RT& rhs) const noexcept {
		return i()==rhs.i() && ptr()==rhs.ptr();
	}
	//! comparison to other iterator
	inline bool operator==(const RT_& rhs) const noexcept {
		return i()==rhs.i() && ptr()==rhs.ptr();
	}
	//! comparison to iterator
	inline bool operator<(const RT& rhs) const noexcept {
		assert(ptr()==rhs.ptr()); return i()<rhs.i();
	}
	//! comparison to other iterator
	inline bool operator<(const RT_& rhs) const noexcept {
		assert(ptr()==rhs.ptr()); return i()<rhs.i();
	}


	/** derived comparison operators
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
		assert(lhs.ptr()==rhs.ptr()); return lhs.i()-rhs.i();
	}
	//! other iterator difference
	inline friend ptrdiff_t operator-(const RT& lhs, const RT_& rhs) noexcept {
		assert(lhs.ptr()==rhs.ptr()); return lhs.i()-rhs.i();
	}


	/** @name streaming
	 */
	//! streaming operator
	friend std::ostream& operator<<(std::ostream& os, const lm_tVecItr_b& i) noexcept {
		return (os << "vi(" << i.ptr() << "," << i.i() << ")");
	}

protected:
	/** @name member variables
	 */
	mutable VT vec_;	//!< vector, i.e. what is dereferenced
};



/**
 * Vector const_iterator class. Can be constructed from other vector iterators or const_iterators.
 * Const Correctness is ensured by returning const rows or columns when dereferencing thereby making
 * changes in the 'host' matrix impossible.
 */
template<class MT, class VT>
class lm_c_tVecItr final: public lm_tVecItr_b<MT,VT,lm_c_tVecItr<MT,VT>,lm_tVecItr<MT,VT>>,
		   	  public std::iterator<std::random_access_iterator_tag,	// iterator_category
	       				       const VT,			// value_type
					       ptrdiff_t,			// difference_type
					       const VT*,			// pointer
					       const VT&> {			// reference
public:
	/** @name types
	 */
	typedef lm_c_tVecItr<MT,VT> c_tVecItr;		//!< vector const_iterator
	typedef lm_tVecItr<MT,VT> tVecItr;		//!< vector iterator
	typedef lm_cr_tVecItr<MT,VT> cr_tVecItr;	//!< reverse vector const_iterator
	typedef lm_r_tVecItr<MT,VT> r_tVecItr;		//!< reverse vector iterator


	/** @name constructors
	 */
	using lm_tVecItr_b<MT,VT,c_tVecItr,tVecItr>::lm_tVecItr_b;
	//! default constructor
	lm_c_tVecItr()=default;
	//! constructor from vector const_iterator
	lm_c_tVecItr(const c_tVecItr& inp) noexcept: c_tVecItr(inp.ptr(),inp.i()) {}
	//! constructor from vector iterator
	lm_c_tVecItr(const tVecItr& inp) noexcept: c_tVecItr(inp.ptr(),inp.i()) {}
	//! constructor from reverse vector const_iterator
	lm_c_tVecItr(const cr_tVecItr& inp) noexcept: c_tVecItr(inp.ptr(),inp.i()+1) {}
	//! constructor from reverse vector iterator
	lm_c_tVecItr(const r_tVecItr& inp) noexcept: c_tVecItr(inp.ptr(),inp.i()+1) {}


	/** @name assignment
	 */
	//! assignment operator from vector const_iterator
	inline c_tVecItr& operator=(const c_tVecItr& rhs) noexcept {
		this->vec_.ptr_ = const_cast<MT*>(rhs.ptr()); this->vec_.i_ = rhs.i();
		return *this;
	}
	//! assignment operator from vector iterator
	inline c_tVecItr& operator=(const tVecItr& rhs) noexcept {
		this->vec_.ptr_ = const_cast<MT*>(rhs.ptr()); this->vec_.i_ = rhs.i();
		return *this;
	}
	//! assignment operator from reverse vector const_iterator
	inline c_tVecItr& operator=(const cr_tVecItr& rhs) noexcept {
		this->vec_.ptr_ = const_cast<MT*>(rhs.ptr()); this->vec_.i_ = rhs.i()+1;
		return *this;
	}
	//! assignment operator from reverse vector iterator
	inline c_tVecItr& operator=(const r_tVecItr& rhs) noexcept {
		this->vec_.ptr_ = const_cast<MT*>(rhs.ptr()); this->vec_.i_ = rhs.i()+1;
		return *this;
	}
	//! swap function
	inline friend void swap(c_tVecItr& lhs, c_tVecItr& rhs) noexcept {
		swap_(lhs.vec_,rhs.vec_);
	}

	
	/** @name information
	 */
	//! const pointer to 'host' matrix
	inline const MT* ptr() const noexcept { return this->vec_.ptr(); }


	/** @name const dereference
	 */
	//! const dereference operator returning const vector
	inline const VT& operator*() const noexcept {
		assert(this->i()>=0);
		assert(this->i()<this->vec_.iL());
		return this->vec_;
	}
	//! const linear indexing operator returning const vector
	inline const VT operator[](const ptrdiff_t i) const noexcept {
		assert(this->i()+i>=0);
		assert(this->i()+i<ptrdiff_t(this->vec_.iL()));
		return *(*this+i);
	}
	//! const arrow operator returning const vector
	inline const VT* operator->() const noexcept {
		return &(operator*());
	}


	/** @name distance
	 */
	//! distance between two vector const_iterators
	inline friend ptrdiff_t distance(const c_tVecItr& lhs, const c_tVecItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
	//! distance between a vector const_iterator and a vector iterator
	inline friend ptrdiff_t distance(const c_tVecItr& lhs, const tVecItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
};



/**
 * Vector iterator class. Can be constructed from other (reverse) vector iterators
 */
template<class MT, class VT>
class lm_tVecItr final: public lm_tVecItr_b<MT,VT,lm_tVecItr<MT,VT>,lm_c_tVecItr<MT,VT>>,
			public std::iterator<std::random_access_iterator_tag,	// iterator_category
	       				     VT,				// value_type
					     ptrdiff_t,				// difference_type
					     VT*,				// pointer
					     VT&> {				// reference
public:
	/** @name types
	 */
	typedef lm_c_tVecItr<MT,VT> c_tVecItr;		//!< vector const_iterator
	typedef lm_tVecItr<MT,VT> tVecItr;		//!< vector iterator
	typedef lm_r_tVecItr<MT,VT> r_tVecItr;		//!< reverse vector iterator


	/** @name constructors
	 */
	/** constructor from matrix
	 * @param ptr pointer to the host matrix
	 * @param i row or column index in the matrix
	 */
	inline explicit lm_tVecItr(MT* ptr=nullptr, const ptrdiff_t i=0) noexcept:
		lm_tVecItr_b<MT,VT,tVecItr,c_tVecItr>(ptr,i) {}
	//! constructor from vector iterator
	lm_tVecItr(const tVecItr& inp) noexcept: tVecItr(inp.ptr(),inp.i()) {}
	//! constructor from vector const_iterator
	lm_tVecItr(const r_tVecItr& inp) noexcept: tVecItr(inp.ptr(),inp.i()+1) {}

	
	/** @name assignment
	 */
	//! assignment operator from vector iterator
	inline tVecItr& operator=(const tVecItr& rhs) noexcept {
		this->vec_.ptr_ = const_cast<MT*>(rhs.ptr()); this->vec_.i_ = rhs.i();
		return *this;
	}
	//! assignment operator from reverse vector iterator
	inline tVecItr& operator=(const r_tVecItr& rhs) noexcept {
		this->vec_.ptr_ = const_cast<MT*>(rhs.ptr()); this->vec_.i_ = rhs.i()+1;
		return *this;
	}
	//! swap function
	inline friend void swap(tVecItr& lhs, tVecItr& rhs) noexcept {
		swap_(lhs.vec_,rhs.vec_);
	}
	

	/** @name information
	 */
	//! pointer to 'host' matrix
	inline MT* ptr() const noexcept { return this->vec_.ptr(); }


	/** @name dereference
	 */
	//! dereference operator returning vector
	inline VT& operator*() const noexcept {
		assert(this->i()>=0);
		assert(this->i()<this->vec_.iL());
		return this->vec_;
	}
	//! linear indexing operator returning vector
	inline VT operator[](const ptrdiff_t i) const noexcept {
		assert(this->i()+i>=0);
		assert(this->i()+i<ptrdiff_t(this->vec_.iL()));
		return *(*this+i);
	}
	//! arrow operator returning vector
	inline VT* operator->() const noexcept {
		return &(operator*());
	}


	/** distance
	 */
	//! distance between two vector iterators
	inline friend ptrdiff_t distance(const tVecItr& lhs, const tVecItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
	//! distance between a vector iterator and a vector const_iterator
	inline friend ptrdiff_t distance(const tVecItr& lhs, const c_tVecItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
};



/**
 * Base class for reverse vector iterators. I.e. reverse row and column iterators.
 * These must be implemented manually because vector iterators are 'stashing' iterators.
 * Template arguments are: \n
 * - MT: 'host' matrix type
 * - VT: vector type, i.e. lm_tRow or lm_tCol
 * - RT: return type, i.e. the reverse iterator inheriting from this
 * - RT_: other return type, i.e. the reverse const_iterator going with the iterator type (or vice versa)
 */
template<class MT, class VT, class RT, class RT_>
class lm_r_tVecItr_b {
public:
	/** constructors
	 */
	/** constructor from matrix
	 * @param ptr pointer to the host matrix
	 * @param i row or column index in the matrix
	 */
	inline explicit lm_r_tVecItr_b(const MT* ptr=nullptr,
		const ptrdiff_t i=0) noexcept: vec_(ptr,i) {}


	/** @name information
	 */
	//! row or column index in the 'host' matrix
	inline ptrdiff_t i() const noexcept { return vec_.i(); }
	//! the number of rows or columns in the 'host' matrix
	inline ptrdiff_t iL() const noexcept { return vec_.iL(); }
	//! const pointer to the 'host' matrix
	inline const MT* ptr() const noexcept { return vec_.ptr(); }


	/** @name casts
	 */
	//! casts to signed integer, i.e. the row or column index in the 'host' matrix
	inline explicit operator ptrdiff_t() const { return i(); }
	//! casts to integer, i.e. the row or column index in the 'host' matrix
	inline explicit operator size_t() const { return size_t(i()); }


	/** @name base iterator
	 */
	//! return the base forward iterator
	RT base() const noexcept { return RT(ptr(),i()+1); }


	/** @name increment, decrement
	 */
	//! iterator post-increment
	inline RT& operator++() noexcept {
		--this->vec_.i_; return *static_cast<RT*>(this);
	}
	//! iterator pre-increment
	inline RT operator++(int) noexcept {
		RT res(*static_cast<RT*>(this)); ++(*this); return res;
	}
	//! iterator post-decrement
	inline RT& operator--() noexcept {
		++this->vec_.i_; return *static_cast<RT*>(this);
	}
	//! iterator pre-decrement
	inline RT operator--(int) noexcept {
		RT res(*static_cast<RT*>(this)); --(*this); return res;
	}
	

	/** @name arithmetic operators
	 */
	//! iterator +=
	inline RT& operator+=(const ptrdiff_t rhs) noexcept {
		this->vec_.i_-=rhs;
		return *static_cast<RT*>(this);
	}
	//! iterator +
	inline RT operator+(const ptrdiff_t rhs) const noexcept {
		return RT(*static_cast<const RT*>(this))+=rhs;
	}
	//! iterator -=
	inline RT& operator-=(const ptrdiff_t rhs) noexcept {
		this->vec_.i_+=rhs;
		return *static_cast<RT*>(this);
	}
	//! iterator -
	inline RT operator-(const ptrdiff_t rhs) const noexcept {
		return RT(*static_cast<const RT*>(this))-=rhs;
	}


	/** @name arithmetic as rhs
	 */
	//! iterator +
	inline friend RT operator+(const ptrdiff_t lhs, const RT& rhs) noexcept {
		return rhs+lhs;
	}


	/** @name fundamental comparison operators
	 */
	//! comparison to iterator
	inline bool operator==(const RT& rhs) const noexcept {
		return i()==rhs.i() && ptr()==rhs.ptr();
	}
	//! comparison to other iterator
	inline bool operator==(const RT_& rhs) const noexcept {
		return i()==rhs.i() && ptr()==rhs.ptr();
	}
	//! comparison to iterator
	inline bool operator<(const RT& rhs) const noexcept {
		assert(ptr()==rhs.ptr()); return i()>rhs.i();
	}
	//! comparison to other iterator
	inline bool operator<(const RT_& rhs) const noexcept {
		assert(ptr()==rhs.ptr()); return i()>rhs.i();
	}


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
		assert(lhs.ptr()==rhs.ptr()); return rhs.i()-lhs.i();
	}
	//! other iterator difference
	inline friend ptrdiff_t operator-(const RT& lhs, const RT_& rhs) noexcept {
		assert(lhs.ptr()==rhs.ptr()); return rhs.i()-lhs.i();
	}


	/** @name streaming
	 */
	//! streaming operator
	friend std::ostream& operator<<(std::ostream& os, const lm_r_tVecItr_b& i) noexcept {
		return (os << "rvi(" << i.ptr() << "," << i.i() << ")");
	}

protected:
	/** @name member variables
	 */
	mutable VT vec_;	//!< vector, i.e. what is dereferenced
};


/**
 * Reverse vector const_iterator class. Can be constructed from other (reverse) vector iterators or const_iterators.
 * Const Correctness is ensured by returning const rows or columns when dereferencing thereby making
 * changes in the 'host' matrix impossible.
 */
template<class MT, class VT>
class lm_cr_tVecItr final: public lm_r_tVecItr_b<MT,VT,lm_cr_tVecItr<MT,VT>,lm_r_tVecItr<MT,VT>>,
			   public std::iterator<std::random_access_iterator_tag,// iterator_category
	       					const VT,			// value_type
						ptrdiff_t,			// difference_type
						const VT*,			// pointer
						const VT&> {			// reference
public:
	/** @name types
	 */
	typedef lm_cr_tVecItr<MT,VT> cr_tVecItr;	//!< reverse vector const_iterator
	typedef lm_r_tVecItr<MT,VT> r_tVecItr;		//!< reverse vector iterator
	typedef lm_c_tVecItr<MT,VT> c_tVecItr;		//!< vector const_iterator
	typedef lm_tVecItr<MT,VT> tVecItr;		//!< vector iterator


	/** @name constructor
	 */
	using lm_r_tVecItr_b<MT,VT,cr_tVecItr,r_tVecItr>::lm_r_tVecItr_b;
	//! default constructor
	lm_cr_tVecItr()=default;
	//! constructor from reverse vector const_iterator
	lm_cr_tVecItr(const cr_tVecItr& inp) noexcept: cr_tVecItr(inp.ptr(),inp.i()) {}
	//! constructor from reverse vector iterator
	lm_cr_tVecItr(const r_tVecItr& inp) noexcept: cr_tVecItr(inp.ptr(),inp.i()) {}
	//! constructor from vector const_iterator
	lm_cr_tVecItr(const c_tVecItr& inp) noexcept: cr_tVecItr(inp.ptr(),inp.i()-1) {}
	//! constructor from vector iterator
	lm_cr_tVecItr(const tVecItr& inp) noexcept: cr_tVecItr(inp.ptr(),inp.i()-1) {}


	/** @name assignment
	 */
	//! assignment operator from reverse vector const_iterator
	inline cr_tVecItr& operator=(const cr_tVecItr& rhs) noexcept {
		this->vec_.ptr_ = const_cast<MT*>(rhs.ptr()); this->vec_.i_ = rhs.i();
		return *this;
	}
	//! assignment operator from reverse vector iterator
	inline cr_tVecItr& operator=(const r_tVecItr& rhs) noexcept {
		this->vec_.ptr_ = const_cast<MT*>(rhs.ptr()); this->vec_.i_ = rhs.i();
		return *this;
	}
	//! assignment operator from vector const_iterator
	inline cr_tVecItr& operator=(const c_tVecItr& rhs) noexcept {
		this->vec_.ptr_ = const_cast<MT*>(rhs.ptr()); this->vec_.i_ = rhs.i()-1;
		return *this;
	}
	//! assignment operator from vector iterator
	inline cr_tVecItr& operator=(const tVecItr& rhs) noexcept {
		this->vec_.ptr_ = const_cast<MT*>(rhs.ptr()); this->vec_.i_ = rhs.i()-1;
		return *this;
	}
	//! swap function
	inline friend void swap(cr_tVecItr& lhs, cr_tVecItr& rhs) noexcept {
		swap_(lhs.vec_,rhs.vec_);
	}


	/** @name information
	 */
	//! const pointer to 'host' matrix
	inline const MT* ptr() const noexcept { return this->vec_.ptr(); }


	/** @name const dereference
	 */
	//! const dereference operator returning const vector
	inline const VT& operator*() const noexcept {
		assert(this->i()>=0);
		assert(this->i()<this->vec_.iL());
		return this->vec_;
	}
	//! const linear indexing operator returning const vector
	inline const VT operator[](const ptrdiff_t i) const noexcept {
		assert(this->i()-i>=0);
		assert(this->i()-i<ptrdiff_t(this->vec_.iL()));
		return *(*this+i);
	}
	//! const arrow operator returning const vector
	inline const VT* operator->() const noexcept {
		return &(operator*());
	}


	/** @name distance
	 */
	//! distance between two reverse vector const_itertators
	inline friend ptrdiff_t distance(const cr_tVecItr& lhs, const cr_tVecItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
	//! distance between a reverse vector const_itertators and a reverse vector iterator
	inline friend ptrdiff_t distance(const cr_tVecItr& lhs, const r_tVecItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
};


/**
 * Reverse vector iterator class. Can be constructed from other (reverse) vector iterators
 */
template<class MT, class VT>
class lm_r_tVecItr final: public lm_r_tVecItr_b<MT,VT,lm_r_tVecItr<MT,VT>,lm_cr_tVecItr<MT,VT>>,
			  public std::iterator<std::random_access_iterator_tag,	// iterator_category
	       				       VT,				// value_type
					       ptrdiff_t,			// difference_type
					       VT*,				// pointer
					       VT&> {				// reference
public:
	/** @name types
	 */
	typedef lm_cr_tVecItr<MT,VT> cr_tVecItr;	//!< reverse vector const_iterator
	typedef lm_r_tVecItr<MT,VT> r_tVecItr;		//!< reverse vector iterator
	typedef lm_c_tVecItr<MT,VT> c_tVecItr;		//!< vector const_iterator
	typedef lm_tVecItr<MT,VT> tVecItr;		//!< vector iterator


	/** @name constructor
	 */
	/** constructor from matrix
	 * @param ptr pointer to the host matrix
	 * @param i row or column index in the matrix
	 */
	inline explicit lm_r_tVecItr(MT* ptr=nullptr, const ptrdiff_t i=0) noexcept:
		lm_r_tVecItr_b<MT,VT,r_tVecItr,cr_tVecItr>(ptr,i) {}
	//! constructor from reverse vector iterator
	lm_r_tVecItr(const r_tVecItr& inp) noexcept: r_tVecItr(inp.ptr(),inp.i()) {}
	//! constructor from vector iterator
	lm_r_tVecItr(const tVecItr& inp) noexcept: r_tVecItr(inp.ptr(),inp.i()-1) {}


	/** @name assignment
	 */
	//! assignment operator from reverse vector iterator
	inline r_tVecItr& operator=(const r_tVecItr& rhs) noexcept {
		this->vec_.ptr_ = const_cast<MT*>(rhs.ptr()); this->vec_.i_ = rhs.i();
		return *this;
	}
	//! assignment operator from vector iterator
	inline r_tVecItr& operator=(const tVecItr& rhs) noexcept {
		this->vec_.ptr_ = const_cast<MT*>(rhs.ptr()); this->vec_.i_ = rhs.i()-1;
		return *this;
	}
	//! swap function
	inline friend void swap(r_tVecItr& lhs, r_tVecItr& rhs) noexcept {
		swap_(lhs.vec_,rhs.vec_);
	}
	

	/** @name information
	 */
	//! pointer to 'host' matrix
	inline MT* ptr() const noexcept { return this->vec_.ptr(); }


	/** @name dereference
	 */
	//! dereference operator returning vector
	inline VT& operator*() const noexcept {
		assert(this->i()>=0);
		assert(this->i()<this->vec_.iL());
		return this->vec_;
	}
	//! linear indexing operator returning vector
	inline VT operator[](const ptrdiff_t i) const noexcept {
		assert(this->i()-i>=0);
		assert(this->i()-i<ptrdiff_t(this->vec_.iL()));
		return *(*this+i);
	}
	//! arrow operator returning vector
	inline VT* operator->() const noexcept {
		return &(operator*());
	}
	

	/** @name distance
	 */
	//! distance between two reverse vector iterators
	inline friend ptrdiff_t distance(const r_tVecItr& lhs, const r_tVecItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
	//! distance between a reverse vector iterator and a reverse vector const_iterator
	inline friend ptrdiff_t distance(const r_tVecItr& lhs, const cr_tVecItr& rhs) noexcept {
		return std::abs(rhs-lhs);
	}
};

#endif // _LM_TVECITR_

/** @}
 */
