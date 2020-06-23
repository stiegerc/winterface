// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup libmat
 * @{
 */

#ifndef _LM_TCOL_
#define _LM_TCOL_

#include "lm_ref_tArray.h"
#include <algorithm>

template<class TT, class FT, class CT>
class lm_tMat;
template<class TT, class FT, class CT>
class lm_tRow;
template<class MT, class VT, class RT, class RT_>
class lm_tVecItr_b;
template<class MT, class VT>
class lm_c_tVecItr;
template<class MT, class VT>
class lm_tVecItr;
template<class MT, class VT, class RT, class RT_>
class lm_r_tVecItr_b;
template<class MT, class VT>
class lm_cr_tVecItr;
template<class MT, class VT>
class lm_r_tVecItr;


/**
 * The column class as hosted by lm_tMat.
 * The api is equivalent to matrices where applicable due to the common inheritance from lm_tArray.
 * This is a class you are not supposed to actively use since it can only be created by a matrix
 * because its only constructor is private. It is only ever used by the functions cAt, cFront, cBack
 * and indirectly by the column iterator. The only purpose is as a window into a single column of the hosting matrix.
 * The reason this class works with a pointer to the 'host' matrix rather than a raw pointer to the data
 * directly, is that this way iterator invalidation can be minimized.
 * Template arguments are: \n
 * - TT: the underlying type of this, either real (if this is a fCol) or complex (if this is a cCol) \n
 * - FT: the real type \n
 * - CT: the complex type\n
 */
template<class TT, class FT, class CT>
class lm_tCol final: public lm_ref_tArray<TT,FT,CT,lm_tCol<TT,FT,CT>> {
public:
	/** @name friend classes
	 */
	friend class lm_tMat<TT,FT,CT>;
	friend class lm_tVecItr_b<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>,
	       lm_c_tVecItr<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>>,lm_tVecItr<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>>>;
	friend class lm_c_tVecItr<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>>;
	friend class lm_tVecItr_b<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>,
	       lm_tVecItr<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>>,lm_c_tVecItr<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>>>;
	friend class lm_tVecItr<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>>;
	friend class lm_r_tVecItr_b<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>,
	       lm_cr_tVecItr<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>>,lm_r_tVecItr<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>>>;
	friend class lm_cr_tVecItr<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>>;
	friend class lm_r_tVecItr_b<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>,
	       lm_r_tVecItr<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>>,lm_cr_tVecItr<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>>>;
	friend class lm_r_tVecItr<lm_tMat<TT,FT,CT>,lm_tCol<TT,FT,CT>>;

	/** @name types
	 */
	typedef lm_tArray<FT,FT,CT> fArray;		//!< a real array
	typedef lm_tArray<CT,FT,CT> cArray;		//!< a complex array
	typedef lm_tRow<TT,FT,CT> tRow;			//!< a row, either real or complex
	typedef lm_tRow<FT,FT,CT> fRow;			//!< a real row
	typedef lm_tRow<CT,FT,CT> cRow;			//!< a complex row
	typedef lm_tCol<TT,FT,CT> tCol;			//!< a column, either real or complex
	typedef lm_tCol<FT,FT,CT> fCol;			//!< a real column
	typedef lm_tCol<CT,FT,CT> cCol;			//!< a complex column
	typedef lm_tMat<TT,FT,CT> tMat;			//!< a matrix, either real or complex
	typedef lm_tMat<FT,FT,CT> fMat;			//!< a real matrix
	typedef lm_tMat<CT,FT,CT> cMat;			//!< a complex matrix
	typedef lm__::lm_size lm_size;			//!< a size struct holding M and N


	/** @name copy constructor and destructor
	 */
	//! copy constructor
	lm_tCol(const lm_tCol& inp) noexcept;
	//! destructor
	inline ~lm_tCol() noexcept { if (i()==PTRDIFF_MAX) delete ptr_; }
	

private:
	/** @name constructor from matrix
	 */
	//! hidden constructor from matrix
	inline explicit lm_tCol(const tMat* ptr=nullptr, const ptrdiff_t i=0) noexcept:
		ptr_(const_cast<tMat*>(ptr)), i_(i) {}

public:
	/** @name assignment
	 */
	//! assignment from a real number
	tCol& operator=(const FT rhs) noexcept;
	//! assignment from a complex number
	tCol& operator=(const CT& rhs) noexcept;
	//! assignment from a real row
	tCol& operator=(const fRow& rhs) noexcept;
	//! assignment from a complex row
	tCol& operator=(const cRow& rhs) noexcept;
	//! assignment from a real array
	tCol& operator=(const fArray& rhs) noexcept;
	//! assignment from a complex array
	tCol& operator=(const cArray& rhs) noexcept;
	//! assignment from a column
	tCol& operator=(const tCol& rhs) noexcept;
	//! swap function, swapping matrix contents
	inline friend void swap(tCol&& lhs, tCol&& rhs) noexcept {
		assert(lhs.M()==rhs.M());
		for (auto il=lhs.begin(),ir=rhs.begin(),e=lhs.end(); il!=e; ++il,++ir)
			std::iter_swap(il,ir);
	}
	//! swap function, swapping hosts
	inline friend void swap_(tCol& lhs, tCol& rhs) noexcept {
		using std::swap;
		swap(lhs.i_,rhs.i_);
		swap(lhs.ptr_,rhs.ptr_);
	}

	/** @name data access
	 */
	//! const pointer to underlying data
	inline const TT* data() const noexcept {
		return this->ptr()->data() + (this->i()==PTRDIFF_MAX ? 0: this->i()*M());
	}
	//! pointer to underlying data
	inline TT* data() noexcept {
		return const_cast<TT*>(static_cast<const lm_tCol*>(this)->data());
	}
	//! const pointer to the 'host' matrix
	const tMat* ptr() const noexcept { return this->ptr_; }
	//! pointer to the 'host' matrix
	tMat* ptr() noexcept { return const_cast<tMat*>(static_cast<const lm_tCol*>(this)->ptr()); }


	/** @name conversion
	 */
	//! returns this as a copy (matrix)
	tMat copy() const noexcept;
	//! returns this as a real copy (matrix)
	fMat fcopy() const noexcept;
	//! returns this as a complex copy (matrix)
	cMat ccopy() const noexcept;


	/** @name information
	 */
	//! returns the dimensions of this
	inline lm_size S() const noexcept { return {M(),N()}; }
	//! returns the number of rows
	inline size_t M() const noexcept { return this->ptr()->M(); }
	//! returns the number of columns
	inline size_t N() const noexcept { return 1; }
	//! returns the increment in the data between successive elements
	inline size_t incr() const noexcept { return 1; }
	//! returns the column position in the 'host' matrix
	ptrdiff_t i() const noexcept { return i_; }
	//! returns the number of columns in the 'host' matrix
	inline ptrdiff_t iL() const noexcept { return this->ptr()->N(); }
	

	/** @name casts
	 */
	//! casts to signed integer, i.e. the column index in the 'host' matrix
	inline explicit operator ptrdiff_t() const { return i(); }
	//! casts to integer, i.e. the column index in the 'host' matrix
	inline explicit operator size_t() const { return size_t(i()); }

	/** @name matrix arithmetic
	 */
	//! matrix product with a real array
	tMat prod(const fArray& inp) const noexcept;
	//! matrix product with a complex array
	cMat prod(const cArray& inp) const noexcept;


	/** @name printing
	 */
	//! write this to binary file
	void writeToFile(const std::string& fileName, const bool noheader=false) const;


private:
	/** @name member variables
	 */
	tMat* ptr_;		//!< pointer to the 'host' matrix
	ptrdiff_t i_;		//!< column index
};

#endif // _LM_TCOL_

/** @}
 */
