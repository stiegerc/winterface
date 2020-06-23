// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup libmat
 * @{
 */

#ifndef _LM_TARRAY_
#define _LM_TARRAY_

#include "lm_defs.h"
#include "lm_tItr.h"
#include "lm_ops.h"
#include <string>
#include <ostream>
#include <cassert>
#include <iterator>
#include <algorithm>


template<class TT, class FT, class CT>
class lm_tMat;


/**
 * The parent class for matrices, rows and columns.
 * Implements all features common to matrices, rows and columns, thus enabling polymorphism between these types.
 * Whenever treatment as an array is sufficient, this class should be used over its children. However it is important
 * to realize that the main class is the matrix, since rows and columns cannot exist independently and are hosted by
 * a 'host' matrix where they represent rows or columns of this matrix. However standalone rows or columns are still
 * possible in the form of an appropriately shaped matrix. The number of rows is referred to as M, the number of
 * columns as N, the linear size, i.e. the total number of elements, as L. Since the underlying data is arranged
 * COLUMN-MAJOR (due to calls to FORTRAN routines), working over columns should be preferred to over rows wherever
 * possible. Template arguments are: \n
 * - TT: the underlying type of this, either real (if this is an fArray) or complex (if this is a cArray) \n
 * - FT: the real type \n
 * - CT: the complex type \n
 */
template<class TT, class FT, class CT>
class lm_tArray {
public:
	/** @name types
	 */
	typedef lm_tItr<TT> tItr;			//!< iterator using custom increment
	typedef lm_c_tItr<TT> c_tItr;			//!< const_iterator using custom increment
	typedef std::reverse_iterator<tItr> r_tItr;	//!< reverse iterator using custom increment
	typedef std::reverse_iterator<c_tItr> cr_tItr;	//!< reverse const_iterator using custom increment
	typedef lm_tArray<TT,FT,CT> tArray;		//!< an array, either real or complex
	typedef lm_tArray<FT,FT,CT> fArray;		//!< a real array
	typedef lm_tArray<CT,FT,CT> cArray;		//!< a complex array
	typedef lm_tMat<TT,FT,CT> tMat;			//!< a matrix, either real or complex
	typedef lm_tMat<FT,FT,CT> fMat;			//!< a real matrix
	typedef lm_tMat<CT,FT,CT> cMat;			//!< a complex matrix
	typedef lm__::lm_size lm_size;			//!< a size struct holding M and N


	/** @name initialization
	 */
	//! virtual destructor
	virtual ~lm_tArray() noexcept {}

	/** @name data access
	 */
	/** linear indexed data access
	 * @param i linear index
	 */
	const TT& operator[](const size_t i) const noexcept {
		assert(i<L());
		return data()[i*incr()];
	}
	/** linear indexed data access
	 * @param i linear index
	 */
	inline TT& operator[](const size_t i) noexcept {
		return const_cast<TT&>(static_cast<const tArray*>(this)->operator[](i));
	}
	/** indexed data access
	 * @param m row index
	 * @param n column index
	 */
	inline const TT& operator()(const size_t m, const size_t n) const noexcept {
		assert(m<M()); assert(n<N());
		return (*this)[n*M()+m];
	}
	/** indexed data access
	 * @param m row index
	 * @param n column index
	 */
	inline TT& operator()(const size_t m, const size_t n) noexcept {
		return const_cast<TT&>(static_cast<const tArray*>(this)->operator()(m,n));
	}

	//! pointer to underlying data
	virtual TT* data() noexcept=0;
	//! const pointer to underlying data
	virtual const TT* data() const noexcept=0;
	//! pointer to the 'host' matrix
	virtual tMat* ptr() noexcept=0;
	//! const pointer to the 'host' matrix
	virtual const tMat* ptr() const noexcept=0;
	
	//! direct access to the first element
	inline TT& front() noexcept { return (*this)[0]; }
	//! const direct access to the first element
	inline const TT& front() const noexcept { return (*this)[0]; }
	//! direct access to the last element
	inline TT& back() noexcept { return (*this)[this->L()-1]; }
	//! const direct access to the last element
	inline const TT& back() const noexcept { return (*this)[this->L()-1]; }


	/** @name iterators
	 */
	//! iterator to the first element
	inline tItr begin() noexcept { return tItr(data(),incr()); }
	//! past the end iterator
	inline tItr end() noexcept { return begin()+L(); }
	//! const_iterator to the first element
	inline c_tItr begin() const noexcept { return c_tItr(data(),incr()); }
	//! past the end const_iterator
	inline c_tItr end() const noexcept { return begin()+L(); }
	//! const_iterator to the first element
	inline c_tItr cbegin() const noexcept { return begin(); }
	//! past the end const_iterator
	inline c_tItr cend() const noexcept { return end(); }
	//! reverse iterator to the last element
	inline r_tItr rbegin() noexcept { return r_tItr(end()); }
	//! past the front reverse iterator
	inline r_tItr rend() noexcept { return r_tItr(begin()); }
	//! reverse const_iterator to the last element
	inline cr_tItr rbegin() const noexcept { return cr_tItr(end()); }
	//! past the front reverse const_iterator
	inline cr_tItr rend() const noexcept { return cr_tItr(begin()); }
	//! reverse const_iterator to the last element
	inline cr_tItr crbegin() const noexcept { return rbegin(); }
	//! past the front reverse const_iterator
	inline cr_tItr crend() const noexcept { return rend(); }

	/** diagonal iterator to the first element
	 * @param os diagonal offset, m ? offset along rows: offset along columns
	 * @param m switch for offset aling rows(true), columns(false)
	 */
	inline tItr dbegin(const size_t os=0, const bool m=true) noexcept {
		assert((empty() && !os) || (m && os<M()) || (!m && os<N()));
		return tItr(data()+(m?os:os*M()),M()+1);
	}
	/** past the end diagonal iterator
	 * @param os diagonal offset, m ? offset along rows: offset along columns
	 * @param m switch for offset aling rows(true), columns(false)
	 */
	inline tItr dend(const size_t os=0, const bool m=true) noexcept {
		return dbegin(os,m) + (m ? std::min(M()-os,N()): std::min(M(),N()-os));
	}
	/** diagonal const_iterator to the first element
	 * @param os diagonal offset, m ? offset along rows: offset along columns
	 * @param m switch for offset aling rows(true), columns(false)
	 */
	inline c_tItr dbegin(const size_t os=0, const bool m=true) const noexcept {
		assert((empty() && !os) || (m && os<M()) || (!m && os<N()));
		return c_tItr(data()+(m?os:os*M()),M()+1);
	}
	/** past the end diagonal const_iterator
	 * @param os diagonal offset, m ? offset along rows: offset along columns
	 * @param m switch for offset aling rows(true), columns(false)
	 */
	inline c_tItr dend(const size_t os=0, const bool m=true) const noexcept {
		return dbegin(os,m) + (m ? std::min(M()-os,N()): std::min(M(),N()-os));
	}
	/** diagonal const_iterator to the first element
	 * @param os diagonal offset, m ? offset along rows: offset along columns
	 * @param m switch for offset aling rows(true), columns(false)
	 */
	inline c_tItr cdbegin(const size_t os=0, const bool m=true) const noexcept {
		return dbegin(os,m);
	}
	/** past the end diagonal const_iterator
	 * @param os diagonal offset, m ? offset along rows: offset along columns
	 * @param m switch for offset aling rows(true), columns(false)
	 */
	inline c_tItr cdend(const size_t os=0, const bool m=true) const noexcept {
		return dend(os,m);
	}
	/** reverse diagonal iterator to the last element
	 * @param os diagonal offset, m ? offset along rows: offset along columns
	 * @param m switch for offset aling rows(true), columns(false)
	 */
	inline r_tItr rdbegin(const size_t os=0, const bool m=true) noexcept {
		return r_tItr(dend(os,m));
	}
	/** past the front reverse diagonal iterator
	 * @param os diagonal offset, m ? offset along rows: offset along columns
	 * @param m switch for offset aling rows(true), columns(false)
	 */
	inline r_tItr rdend(const size_t os=0, const bool m=true) noexcept {
		return r_tItr(dbegin(os,m));
	}
	/** reverse diagonal const_iterator to the last element
	 * @param os diagonal offset, m ? offset along rows: offset along columns
	 * @param m switch for offset aling rows(true), columns(false)
	 */
	inline cr_tItr rdbegin(const size_t os=0, const bool m=true) const noexcept {
		return cr_tItr(dend(os,m));
	}
	/** past the front reverse diagonal const_iterator
	 * @param os diagonal offset, m ? offset along rows: offset along columns
	 * @param m switch for offset aling rows(true), columns(false)
	 */
	inline cr_tItr rdend(const size_t os=0, const bool m=true) const noexcept {
		return cr_tItr(dbegin(os,m));
	}
	/** reverse diagonal const_iterator to the last element
	 * @param os diagonal offset, m ? offset along rows: offset along columns
	 * @param m switch for offset aling rows(true), columns(false)
	 */
	inline cr_tItr crdbegin(const size_t os=0, const bool m=true) const noexcept {
		return cr_tItr(rdbegin(os,m));
	}
	/** past the front reverse diagonal const_iterator
	 * @param os diagonal offset, m ? offset along rows: offset along columns
	 * @param m switch for offset aling rows(true), columns(false)
	 */
	inline cr_tItr crdend(const size_t os=0, const bool m=true) const noexcept {
		return cr_tItr(rdend(os,m));
	}


	/** @name basic properties
	 */
	//! check whether this is a row
	bool row() const noexcept { return (M()==1)&&N(); }
	//! check whether this is a column
	bool col() const noexcept { return (N()==1)&&M(); }
	//! check whether this is complex
	static bool cpx() noexcept;


	/** @name information
	 */
	//! returns dimensions of this
	virtual lm_size S() const noexcept=0;
	//! returns the number of rows
	virtual size_t M() const noexcept=0;
	//! returns the number of columns
	virtual size_t N() const noexcept=0;
	//! returns the total number of elements
	inline size_t L() const noexcept { return M()*N(); }
	//! returns the total number of elements
	inline size_t size() const noexcept { return L(); }
	//! returns whether this is of size 0
	inline bool empty() const noexcept { return !L(); }
	//! returns the increment in the data between successive elements
	virtual size_t incr() const noexcept=0;
	

	/** @name conversion
	 */
	//! returns this as a copy
	virtual tMat copy() const noexcept=0;
	//! returns this as a real copy
	virtual fMat fcopy() const noexcept=0;
	//! returns this as a complex copy
	virtual cMat ccopy() const noexcept=0;


	/** @name logical
	 */
	//! checks whether this is composed of only 0 and 1
	inline bool logical() const noexcept {
		return !std::all_of(cbegin(),cend(),[](const TT& i){
			return lm__::ops::nz_s(i) && lm__::ops::neq_s(i,FT(1.0));});
	}
	//! returns the NOT of this
	fMat operator~() const noexcept;
	//! returns the elementwise AND of this and a real number
	inline fMat operator&(const FT rhs) const noexcept { return copy()&=rhs; }
	//! returns the elementwise AND of this and a complex number
	inline fMat operator&(const CT& rhs) const noexcept { return copy()&=rhs; }
	//! returns the elementwise AND of this and a real array
	inline fMat operator&(const fArray& rhs) const noexcept { return copy()&=rhs; }
	//! returns the elementwise AND of this and a complex array
	inline fMat operator&(const cArray& rhs) const noexcept { return copy()&=rhs; }
	//! returns the elementwise OR of this and a real number
	inline fMat operator|(const FT rhs) const noexcept { return copy()|=rhs; }
	//! returns the elementwise OR of this and a complex number
	inline fMat operator|(const CT& rhs) const noexcept { return copy()|=rhs; }
	//! returns the elementwise OR of this and a real array
	inline fMat operator|(const fArray& rhs) const noexcept { return copy()|=rhs; }
	//! returns the elementwise OR of this and a complex array
	inline fMat operator|(const cArray& rhs) const noexcept { return copy()|=rhs; }


	/** @name comparison collapsing to bool
	 */
	//! comparison of this to a real number
	inline bool operator==(const FT rhs) const noexcept {
		return std::all_of(cbegin(),cend(),[rhs](const TT& i){return lm__::ops::eq(i,rhs);});
	}
	//! comparison of this to a complex number
	inline bool operator==(const CT& rhs) const noexcept {
		return std::all_of(cbegin(),cend(),[rhs](const TT& i){return lm__::ops::eq(i,rhs);});
	}
	//! comparison of this to a real array
	inline bool operator==(const fArray& rhs) const noexcept {
		return std::equal(cbegin(),cend(),rhs.cbegin(),rhs.cend(),
			[](const TT i, const FT j){return lm__::ops::eq(i,j);});
	}
	//! comparison of this to a complex array
	inline bool operator==(const cArray& rhs) const noexcept {
		if (M()!=rhs.M() || N()!=rhs.N()) return false;
		return std::equal(cbegin(),cend(),rhs.cbegin(),rhs.cend(),
			[](const TT i, const CT j){return lm__::ops::eq(i,j);});
	}
	//! comparison of this to a real number
	inline bool operator!=(const FT rhs) const noexcept { return !(*this==rhs); }
	//! comparison of this to a complex number
	inline bool operator!=(const CT& rhs) const noexcept { return !(*this==rhs); }
	//! comparison of this to a real array
	inline bool operator!=(const fArray& rhs) const noexcept { return !(*this==rhs); }
	//! comparison of this to a complex array
	inline bool operator!=(const cArray& rhs) const noexcept { return !(*this==rhs); }
	//! lexicographical compare to a real number
	inline bool operator<(const FT rhs) const noexcept {
		for (auto i=cbegin(), e=cend(); i!=e; ++i) {
			if (lm__::ops::lt(*i,rhs)) return true;
			if (lm__::ops::gt(*i,rhs)) return false;
		}
		return false;
	}
	//! lexicographical comparison of this to a complex number
	inline bool operator<(const CT& rhs) const noexcept {
		for (auto i=cbegin(), e=cend(); i!=e; ++i) {
			if (lm__::ops::lt(*i,rhs)) return true;
			if (lm__::ops::gt(*i,rhs)) return false;
		}
		return false;
	}
	//! lexicographical comparison of this to a real array
	inline bool operator<(const fArray& rhs) const noexcept {
		assert(L()==rhs.L());
		auto j=rhs.cbegin();
		for (auto i=cbegin(), e=cend(); i!=e; ++i,++j) {
			if (lm__::ops::lt(*i,*j)) return true;
			if (lm__::ops::gt(*i,*j)) return false;
		}
		return false;
	}
	//! lexicographical comparison of this to a complex array
	inline bool operator<(const cArray& rhs) const noexcept {
		assert(L()==rhs.L());
		auto j=rhs.cbegin();
		for (auto i=cbegin(), e=cend(); i!=e; ++i,++j) {
			if (lm__::ops::lt(*i,*j)) return true;
			if (lm__::ops::gt(*i,*j)) return false;
		}
		return false;
	}
	//! lexicographical comparison of this to a real number
	inline bool operator<=(const FT rhs) const noexcept {
		return operator<(rhs)||operator==(rhs);
	}
	//! lexicographical comparison of this to a complex number
	inline bool operator<=(const CT& rhs) const noexcept {
		return operator<(rhs)||operator==(rhs);
	}
	//! lexicographical comparison of this to a real array
	inline bool operator<=(const fArray& rhs) const noexcept {
		return operator<(rhs)||operator==(rhs);
	}
	//! lexicographical comparison of this to a complex array
	inline bool operator<=(const cArray& rhs) const noexcept {
		return operator<(rhs)||operator==(rhs);
	}
	//! lexicographical comparison of this to a real number
	inline bool operator>(const FT rhs) const noexcept {
		return !operator<=(rhs);
	}
	//! lexicographical comparison of this to a complex number
	inline bool operator>(const CT& rhs) const noexcept {
		return !operator<=(rhs);
	}
	//! lexicographical comparison of this to a real array
	inline bool operator>(const fArray& rhs) const noexcept {
		return !operator<=(rhs);
	}
	//! lexicographical comparison of this to a complex array
	inline bool operator>(const cArray& rhs) const noexcept {
		return !operator<=(rhs);
	}
	//! lexicographical comparison of this to a real number
	inline bool operator>=(const FT rhs) const noexcept {
		return !operator<(rhs);
	}
	//! lexicographical comparison of this to a complex number
	inline bool operator>=(const CT& rhs) const noexcept {
		return !operator<(rhs);
	}
	//! lexicographical comparison of this to a real array
	inline bool operator>=(const fArray& rhs) const noexcept {
		return !operator<(rhs);
	}
	//! lexicographical comparison of this to a complex array
	inline bool operator>=(const cArray& rhs) const noexcept {
		return !operator<(rhs);
	}


	/** @name comparison collapsing to logical matrix
	 */
	//! comparison of this to a real number
	fMat eq(const FT rhs) const noexcept;
	//! comparison of this to a complex number
	fMat eq(const CT& rhs) const noexcept;
	//! comparison of this to a real array
	fMat eq(const fArray& rhs) const noexcept;
	//! comparison of this to a complex array
	fMat eq(const cArray& rhs) const noexcept;
	//! comparison of this to a real number
	fMat neq(const FT rhs) const noexcept;
	//! comparison of this to a complex number
	fMat neq(const CT& rhs) const noexcept;
	//! comparison of this to a real array
	fMat neq(const fArray& rhs) const noexcept;
	//! comparison of this to a complex array
	fMat neq(const cArray& rhs) const noexcept;
	//! comparison of this to a real number
	fMat lt(const FT rhs) const noexcept;
	//! comparison of this to a complex number
	fMat lt(const CT& rhs) const noexcept;
	//! comparison of this to a real array
	fMat lt(const fArray& rhs) const noexcept;
	//! comparison of this to a complex array
	fMat lt(const cArray& rhs) const noexcept;
	//! comparison of this to a real number
	fMat leq(const FT rhs) const noexcept;
	//! comparison of this to a complex number
	fMat leq(const CT& rhs) const noexcept;
	//! comparison of this to a real array
	fMat leq(const fArray& rhs) const noexcept;
	//! comparison of this to a complex array
	fMat leq(const cArray& rhs) const noexcept;
	//! comparison of this to a real number
	fMat gt(const FT rhs) const noexcept;
	//! comparison of this to a complex number
	fMat gt(const CT& rhs) const noexcept;
	//! comparison of this to a real array
	fMat gt(const fArray& rhs) const noexcept;
	//! comparison of this to a complex array
	fMat gt(const cArray& rhs) const noexcept;
	//! comparison of this to a real number
	fMat geq(const FT rhs) const noexcept;
	//! comparison of this to a complex number
	fMat geq(const CT& rhs) const noexcept;
	//! comparison of this to a real array
	fMat geq(const fArray& rhs) const noexcept;
	//! comparison of this to a complex array
	fMat geq(const cArray& rhs) const noexcept;
	
	/** @name elementwise arithmetic
	 */
	//! unary operator minus
	tMat operator-() const noexcept;
	//! returns this + a real number
	inline tMat operator+(const FT rhs) const noexcept { return copy()+=rhs; }
	//! returns complex this + a complex number
	inline cMat operator+(const CT& rhs) const noexcept { return ccopy()+=rhs; }
	//! returns elementwise this + a real array
	inline tMat operator+(const fArray& rhs) const noexcept { return copy()+=rhs; }
	//! returns complex elementwise this + a complex array
	inline cMat operator+(const cArray& rhs) const noexcept { return ccopy()+=rhs; }
	//! returns this - a real number
	inline tMat operator-(const FT rhs) const noexcept { return copy()-=rhs; }
	//! returns complex this - a complex number
	inline cMat operator-(const CT& rhs) const noexcept { return ccopy()-=rhs; }
	//! returns elementwise this - a real array
	inline tMat operator-(const fArray& rhs) const noexcept { return copy()-=rhs; }
	//! returns complex elementwise this - a complex array
	inline cMat operator-(const cArray& rhs) const noexcept { return ccopy()-=rhs; }
	//! returns this * a real number
	inline tMat operator*(const FT rhs) const noexcept { return copy()*=rhs; }
	//! returns complex this * a complex number
	inline cMat operator*(const CT& rhs) const noexcept { return ccopy()*=rhs; }
	//! returns elementwise this * a real array
	inline tMat operator*(const fArray& rhs) const noexcept { return copy()*=rhs; }
	//! returns complex elementwise this * a complex array
	inline cMat operator*(const cArray& rhs) const noexcept { return ccopy()*=rhs; }
	//! returns this / a real number
	inline tMat operator/(const FT rhs) const noexcept { return copy()/=rhs; }
	//! returns complex this / a complex number
	inline cMat operator/(const CT& rhs) const noexcept { return ccopy()/=rhs; }
	//! returns elementwise this / a real array
	inline tMat operator/(const fArray& rhs) const noexcept { return copy()/=rhs; }
	//! returns complex elementwise this / a complex array
	inline cMat operator/(const cArray& rhs) const noexcept { return ccopy()/=rhs; }
	//! returns this % a real number
	inline tMat operator%(const FT rhs) const noexcept { return copy()%=rhs; }
	//! returns complex this % a complex number
	inline cMat operator%(const CT& rhs) const noexcept { return ccopy()%=rhs; }
	//! returns elementwise this % a real array
	inline tMat operator%(const fArray& rhs) const noexcept { return copy()%=rhs; }
	//! returns complex elementwise this % a complex array
	inline cMat operator%(const cArray& rhs) const noexcept { return ccopy()%=rhs; }

	/** @name printing
	 */
	//! print this to a std::string
	std::string print(const size_t precision=PPREC__, const size_t blanks=0) const noexcept;
	//! print this to a textfile
	void printToFile(const std::string& fileName, const size_t precision=PPREC__) const;
	//! write this to a binary file
	virtual void writeToFile(const std::string& fileName, const bool noheader) const=0;
};


//! matrix size function
template<class TT, class FT, class CT>
inline lm__::lm_size msize(const lm_tArray<TT,FT,CT>& inp) noexcept { return inp.S(); }


//! streaming operator for arrays
template<class TT, class FT, class CT>
inline std::ostream& operator<<(std::ostream& os, const lm_tArray<TT,FT,CT>& inp) noexcept { return (os<<inp.print()); }


//! returns a real number + an array
template<class TT, class FT, class CT>
inline lm_tMat<TT,FT,CT> operator+(const FT lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept { return rhs+lhs; }
//! returns complex this + a complex number
template<class TT, class FT, class CT>
inline lm_tMat<CT,FT,CT> operator+(const CT& lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept { return rhs+lhs; }
//! returns a real number - an array
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> operator-(const FT lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept;
//! returns complex this - a complex number
template<class TT, class FT, class CT>
lm_tMat<CT,FT,CT> operator-(const CT& lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept;
//! returns a real number * an array
template<class TT, class FT, class CT>
inline lm_tMat<TT,FT,CT> operator*(const FT lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept { return rhs*lhs; }
//! returns complex this * a complex number
template<class TT, class FT, class CT>
inline lm_tMat<CT,FT,CT> operator*(const CT& lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept { return rhs*lhs; }
//! returns a real number / an array
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> operator/(const FT lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept;
//! returns complex this / a complex number
template<class TT, class FT, class CT>
lm_tMat<CT,FT,CT> operator/(const CT& lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept;


#endif // _LM_TARRAY_

/** @}
 */
