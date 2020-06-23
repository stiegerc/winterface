// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup libmat
 * @{
 */

#ifndef _LM_TMAT_
#define _LM_TMAT_

#include "lm_ref_tArray.h"
#include "lm_tCol.h"
#include "lm_tRow.h"
#include "lm_tVecItr.h"
#include "lm_defs.h"
#include <string>
#include <vector>
#include <functional>
#include <ostream>
#include <iostream>
#include <cassert>
#include <cstring>


/**
 * The main matrix class.
 * This class can spawn so-called row or column iterators, which point to entire rows or columns of the matrix.
 * Each such iterator returns a row or column upon dereferencing that is hosted by its 'host' matrix. These do
 * not own their own memory but rather manipulations on them change data in the 'host' matrix. The api on rows,
 * columns and matrices is as equivalent as can be due all of them inheriting from lm_tArray. However, due to
 * the column-major arrangement of the underlying memory, this matrix class is specifically designed as a sort
 * of vector of columns. Working over columns is therefore to be prioritized. However rows, columns and matrices
 * may be intermixed in the sense that one can be assigned to any of the others provided the total number of
 * elements is the same. Additionally, the matrix may have space for more columns than it currently holds,
 * analogous to std::vector. \n
 * Template arguments are: \n
 * - TT: the underlying type of this, either real (if this is an fMat) or complex (if this is a cMat) \n
 * - FT: the real type \n
 * - CT: the complex type \n
 */
template<class TT, class FT, class CT>
class lm_tMat final: public lm_ref_tArray<TT,FT,CT,lm_tMat<TT,FT,CT>> {
public:
	/** @name types
	 */
	typedef lm_tArray<TT,FT,CT> tArray;		//!< an array, either real or complex
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
	typedef lm_tItr<TT> tItr;			//!< iterator using custom increment
	typedef lm_c_tItr<TT> c_tItr;			//!< const_iterator using custom increment
	typedef lm_tVecItr<tMat,tRow> tRowItr;		//!< a row iterator, either real or complex
	typedef lm_c_tVecItr<tMat,tRow> c_tRowItr;	//!< a row const_iterator, either real or complex
	typedef lm_tVecItr<tMat,tCol> tColItr;		//!< a column iterator, either real or complex
	typedef lm_c_tVecItr<tMat,tCol> c_tColItr;	//!< a column const_iterator, either real or complex
	typedef lm_r_tVecItr<tMat,tRow> r_tRowItr;	//!< a row reverse iterator, either real or complex
	typedef lm_cr_tVecItr<tMat,tRow> cr_tRowItr;	//!< a row reverse const_iterator, either real or complex
	typedef lm_r_tVecItr<tMat,tCol> r_tColItr;	//!< a column reverse iterator, either real or complex
	typedef lm_cr_tVecItr<tMat,tCol> cr_tColItr;	//!< a column reverse const_iterator, either real or complex
	typedef lm__::lm_size lm_size;			//!< a size struct holding M and N


	/** @name constructors
	 */
	//! default constructor generating an empty 0x0 matrix
	inline lm_tMat() noexcept: tMat(0,0) {}
	/** main constructor for a matrix of size MxN, memory is allocated but not initialized
	 * @param M the number of rows
	 * @param N the number of columns
	 */
	inline explicit lm_tMat(const size_t M, const size_t N) noexcept:
		M_(M), N_(N), C_(M*N), data_(C_ ? new TT [C_]: nullptr) {}
	//! overload for use with the size struct
	inline explicit lm_tMat(const lm_size& S) noexcept: lm_tMat(S.M,S.N) {}
	/** overload for square matrices
	 * @param M the number of rows/columns
	 */
	inline explicit lm_tMat(const size_t M) noexcept: tMat(M,M) {}
	/** constructor for passing preallocated data to the matrix, will be deallocated in the destructor
	 * @param data pointer to the preallocated data
	 * @param M the number of rows
	 * @param N the number of columns
	 */
	inline explicit lm_tMat(TT* const data, const size_t M, const size_t N) noexcept:
		M_(M), N_(N), C_(M*N), data_(data) {}
	//! constructor from a set of row const_iterators
	inline explicit lm_tMat(c_tRowItr i, const c_tRowItr& e) noexcept:
		tMat(distance(i,e),i->L()) { assert(i<=e); for (auto j=rBegin(); i!=e; ++i,++j) *j=*i; }
	//! constructor from a set of row iterators
	inline explicit lm_tMat(const tRowItr& i, const tRowItr& e) noexcept:
		tMat(c_tRowItr(i),c_tRowItr(e)) {}
	//! constructor from a set of column const_iterators
	explicit lm_tMat(const c_tColItr& i, const c_tColItr& e) noexcept:
		tMat(i->L(),distance(i,e)) { assert(i<=e); memcpy(data_,i->data(),this->L()*sizeof(TT)); }
	//! constructor from a set of column iterators
	inline explicit lm_tMat(const tColItr& i, const tColItr& e) noexcept:
		tMat(c_tColItr(i),c_tColItr(e)) {}
	//! constructor from a list of rows
	inline lm_tMat(const std::initializer_list<const tRow>& init) noexcept:
			tMat(init.size(), init.size() ? init.begin()->N(): 0) {
		std::copy(init.begin(),init.end(),rBegin());
	}
	//! constructor from a list of columns
	inline lm_tMat(const std::initializer_list<const tCol>& init) noexcept:
			tMat(init.size() ? init.begin()->M(): 0, init.size()) {
		std::copy(init.begin(),init.end(),cBegin());
	}
	//! constructor from a list of numbers, result is a column shaped matrix
	inline lm_tMat(const std::initializer_list<TT>& init) noexcept: tMat(init.size(),1) {
		if (!this->empty()) memcpy(data_,init.begin(),this->L()*sizeof(TT));
	}
	/** constructor from a list of numbers
	 * @param init the list of numbers
	 * @param M the number of rows
	 */
	inline lm_tMat(const std::initializer_list<TT>& init, const size_t M) noexcept: tMat(M,init.size()/M) {
		assert(init.size()==this->L());
		if (!this->empty()) memcpy(data_,init.begin(),this->L()*sizeof(TT));
	}
	/** constructor from a list of numbers
	 * @param init the list of numbers
	 * @param M the number of rows
	 * @param N the number of columns
	 */
	inline lm_tMat(const std::initializer_list<TT>& init, const size_t M, const size_t N) noexcept: tMat(M,N) {
		assert(init.size()==this->L());
		if (!this->empty()) memcpy(data_,init.begin(),this->L()*sizeof(TT));
	}
	/** constructor from a file stream
	 * @param file the input stream
	 * @param M the number of rows
	 * @param N the number of columns
	 */
	inline explicit lm_tMat(std::ifstream& file, const size_t M, const size_t N): tMat(M,N) { parse_(file); }
	/** constructor from a std::vector
	 * @param inp the input vector
	 * @param col switch whether to generate a row shaped matrix (false) or a column shaped one (true)
	 */
	inline explicit lm_tMat(const std::vector<TT>& inp, const bool col=true) noexcept:
			tMat(col ? lm_size{inp.size(),1}: lm_size{1,inp.size()}) {
		memcpy(data_,inp.data(),this->L()*sizeof(TT));
	}
	/** constructor from an std::vector
	 * @param inp the input vector
	 * @param col switch to generate a column(true) or row(false) shaped matrix
	 */
	template<class VT>
	inline explicit lm_tMat(const std::vector<VT>& inp, const bool col=true) noexcept:
			tMat(col ? lm_size{inp.size(),1}: lm_size{1,inp.size()}) {
		std::transform(inp.begin(),inp.end(),this->begin(),[](const size_t& i){return TT(i);});
	}
	//! constructor from real and imaginary arrays
	explicit lm_tMat(const fArray& re, const fArray& im) noexcept;
	//! constructor from a real array
	lm_tMat(const fArray& inp) noexcept;
	//! conmstructor from a complex array
	lm_tMat(const cArray& inp) noexcept;
	//! copy constructor
	inline lm_tMat(const tMat& inp) noexcept: tMat(msize(inp)) {
		memcpy(data(),inp.data(),this->L()*sizeof(TT));
	}
	//! move constructor
	inline lm_tMat(tMat&& inp) noexcept: tMat() {
		using std::swap;
		swap(*this,inp);
	}
	//! destructor
	~lm_tMat() noexcept {
		delete[] data_;
	}
	//! constructor from a filename
	explicit lm_tMat(const std::string& fileName);


	/** @name assignment
	 */
	//! swap function
	friend void swap(tMat& lhs, tMat& rhs) noexcept {
		using std::swap;
		swap(lhs.M_,rhs.M_);
		swap(lhs.N_,rhs.N_);
		swap(lhs.C_,rhs.C_);
		swap(lhs.data_,rhs.data_);
	}
	//! assignment operator from a real array
	inline tMat& operator=(const fArray& rhs) noexcept {
		resize(rhs.M(),rhs.N());
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::assign(i,*j++);});
		return *this;
	}
	//! assignment operator from a complex array
	inline tMat& operator=(const cArray& rhs) noexcept {
		resize(rhs.M(),rhs.N());
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::assign(i,*j++);});
		return *this;
	}
	//! assignment operator from a matrix of equal type
	tMat& operator=(const tMat& rhs) noexcept {
		if (this!=&rhs) {
			resize(rhs.M(),rhs.N(),false);
			memcpy(this->data(),rhs.data(),rhs.L()*sizeof(TT));
		}
		return *this;
	}
	//! move assignment operator
	tMat& operator=(tMat&& rhs) noexcept {
		swap(*this,rhs);
		return *this;
	}


	/** @name data access
	 */
	using tArray::operator[];
	//! linear indexed data access
	inline const TT& operator[](const size_t i) const noexcept {
		assert(i<this->L());
		return data()[i];
	}
	//! pointer to underlying data
	inline TT* data() noexcept { return data_; }
	//! const pointer to underlying data
	inline const TT* data() const noexcept { return data_; }
	//! const pointer to the 'host' matrix
	inline const tMat* ptr() const noexcept { return this; }
	//! pointer to the 'host' matrix
	inline tMat* ptr() noexcept { return this; }


	/** @name row, column access
	 */
	/** indexed const access to the rows of the matrix
	 * @param m the row index, 0<=m<M
	 */
	inline const tRow rAt(const size_t m) const noexcept {
		assert(m<M());
		return tRow(this,m);
	}
	/** indexed access to the rows of the matrix
	 * @param m the row index, 0<=m<M
	 */
	inline tRow rAt(const size_t m) noexcept {
		assert(m<M());
		return tRow(this,m);
	}
	//! const direct access to the first row
	inline const tRow rFront() const noexcept { return rAt(0); }
	//! direct access to the first row
	inline tRow rFront() noexcept { return rAt(0); }
	//! const direct access to the last row
	inline const tRow rBack() const noexcept { return rAt(M()-1); }
	//! direct access to the last row
	inline tRow rBack() noexcept { return rAt(M()-1); }
	/** indexed const access to the columns of the matrix
	 * @param n the column index, 0<=n<N
	 */
	inline const tCol cAt(const size_t n) const noexcept {
		assert(n<N());
		return tCol(this,n);
	}
	/** indexed access to the columns of the matrix
	 * @param n the column index, 0<=n<N
	 */
	inline tCol cAt(const size_t n) noexcept {
		assert(n<N());
		return tCol(this,n);
	}
	//! const direct access to the first column
	inline const tCol cFront() const noexcept { return cAt(0); }
	//! direct access to the first column
	inline tCol cFront() noexcept { return cAt(0); }
	//! const direct access to the last column
	inline const tCol cBack() const noexcept { return cAt(N()-1); }
	//! direct access to the last column
	inline tCol cBack() noexcept { return cAt(N()-1); }


	/** @name row iterators
	 */
	//! row iterator to the first row
	inline tRowItr rBegin() noexcept { return tRowItr(this,0); }
	//! past the end row iterator
	inline tRowItr rEnd() noexcept { return tRowItr(this,M()); }
	//! row const_iterator to the first row
	inline c_tRowItr rBegin() const noexcept { return c_tRowItr(this,0); }
	//! past the end row const_iterator
	inline c_tRowItr rEnd() const noexcept { return c_tRowItr(this,M()); }
	//! row const_iterator to the first row
	inline c_tRowItr crBegin() const noexcept { return rBegin(); }
	//! past the end row const_iterator
	inline c_tRowItr crEnd() const noexcept { return rEnd(); }
	//! reverse row iterator to the last element
	inline r_tRowItr rrBegin() noexcept { return r_tRowItr(rEnd()); }
	//! past the front reverse row iterator
	inline r_tRowItr rrEnd() noexcept { return r_tRowItr(rBegin()); }
	//! reverse row const_iterator to the last element
	inline cr_tRowItr rrBegin() const noexcept { return cr_tRowItr(rEnd()); }
	//! past the front reverse row const_iterator
	inline cr_tRowItr rrEnd() const noexcept { return cr_tRowItr(rBegin()); }
	//! reverse row const_iterator to the last element
	inline cr_tRowItr crrBegin() const noexcept { return rrBegin(); }
	//! past the front reverse row const_iterator
	inline cr_tRowItr crrEnd() const noexcept { return rrEnd(); }


	/** @name column iterators
	 */
	//! column iterator to the first column
	inline tColItr cBegin() noexcept { return tColItr(this,0); }
	//! past the end column iterator
	inline tColItr cEnd() noexcept { return tColItr(this,N()); }
	//! column const_iterator to the first column
	inline c_tColItr cBegin() const noexcept { return c_tColItr(this,0); }
	//! past the end column const_iterator
	inline c_tColItr cEnd() const noexcept { return c_tColItr(this,N()); }
	//! column const_iterator to the first column
	inline c_tColItr ccBegin() const noexcept { return cBegin(); }
	//! past the end column const_iterator
	inline c_tColItr ccEnd() const noexcept { return cEnd(); }
	//! reverse column iterator to the last element
	inline r_tColItr rcBegin() noexcept { return r_tColItr(cEnd()); }
	//! past the front column iterator
	inline r_tColItr rcEnd() noexcept { return r_tColItr(cBegin()); }
	//! reverse column const_iterator to the last element
	inline cr_tColItr rcBegin() const noexcept { return cr_tColItr(cEnd()); }
	//! past the front column const_iterator
	inline cr_tColItr rcEnd() const noexcept { return cr_tColItr(cBegin()); }
	//! reverse column const_iterator to the last element
	inline cr_tColItr crcBegin() const noexcept { return rcBegin(); }
	//! past the front column const_iterator
	inline cr_tColItr crcEnd() const noexcept { return rcEnd(); }

	/** @name memory management
	 */
	/** reserve space, may lead to allocating an copying
	 * @param nC the new number of columns
	 */
	tMat& reserve(const size_t nC) noexcept;
	//! shrink allocated memory to the actual number of columns
	tMat& shrink_to_fit() noexcept;
	//! see reserve
	inline tMat& operator<<(size_t nC) noexcept { return reserve(nC); }
	//! move internals by returning the underlying data
	TT* move() noexcept;
	

	/** @name basic properties
	 */
	//! returns whether the matrix is square
	inline bool square() const noexcept { return M()==N(); }
	//! returns whether the matrix is hermitian
	bool hermitian() const noexcept;
	//! returns whether the matrix has only diagonal entries
	bool diag() const noexcept;
	//! returns whether the columns of the matrix form an orthogonal basis
	bool ob() const noexcept;
	//! returns whether the columns of the matrix form an orthonormal basis
	bool onb() const noexcept;


	/** @name information
	 */
	//! returns the dimensions of this
	inline lm_size S() const noexcept { return {M(),N()}; }
	//! returns the number of rows
	inline size_t M() const noexcept { return M_; }
	//! returns the number of columns
	inline size_t N() const noexcept { return N_; }
	//! returns the linear capacity
	inline size_t lcap() const noexcept { return C_; }
	//! returns the columnar capacity
	inline size_t ccap() const noexcept { return !M() ? 0: lcap()/M(); }
	//! returns the columnar capacity
	inline size_t capacity() const noexcept { return ccap(); }
	//! returns the increment in the data between successive elements
	inline size_t incr() const noexcept { return 1; }


	/** @name basic modification
	 */
	//! 'rowizes' the matrix
	inline tMat& R() noexcept { N_*=M_; M_=1; return *this; }
	//! 'columnizes' the matrix
	inline tMat& C() noexcept { M_*=N_; N_=1; return *this; }
	//! transposes the matrix
	tMat& T() noexcept;
	//! clears the matrix by setting the number of columns to 0
	inline tMat& clear() noexcept { N_=0; return *this; }
	/** resizes the matrix
	 * @param nN the new number of columns
	 */
	inline tMat& resize(const size_t nN) noexcept {
		reserve(nN); N_=nN; return *this;
	}
	/** resizes the matrix
	 * @param nM the new number of rows
	 * @param nN the new number of columns
	 * @param cpy switch whether to copy the old data in case of reallocation
	 */
	inline tMat& resize(const size_t nM, const size_t nN, const bool cpy=true) noexcept {
		if (nM*nN>C_) {
			TT* ndata = new TT [nM*nN];
			if (cpy) memcpy(ndata,data_,M_*N_*sizeof(TT));
			delete[] data_; data_ = ndata; C_ = nM*nN;
		}
		M_ = nM, N_ = nN;
		return *this;
	}
	//! resized the matrix, size struct overload
	inline tMat& resize(const lm_size& S) noexcept { return resize(S.M,S.N); }
	/** push_back operator analogous to std::vector, overload for rows
	 * @param inp row to be inserted as a column in the back
	 */
	inline tMat& push_back(const tRow& inp) noexcept {
		assert(inp.L()==M());
		resize(N()+1);
		cBack()=inp;
		return *this;
	}
	/** push_back operator analogous to std::vector, overload for arrays
	 * @param inp array to be inserted as a column in the back
	 */
	inline tMat& push_back(const tArray& inp) noexcept {	
		assert(inp.M()==M());
		reserve(N()+inp.N());
		memcpy(data()+this->L(),inp.data(),inp.L()*sizeof(TT));
		N_+=inp.N();
		return *this;
	}
	//! same as the push_back operator
	inline tMat& operator<<(const tArray& inp) noexcept { return push_back(inp); }
	/** pop_back operator analogous to std::vector
	 * @param n the number of columns to be popped
	 */
	inline tMat& pop_back(const size_t n=1) noexcept {
		N_ -= (n<N() ? n: N()); 
		return *this;
	}
	

	/** @name modification
	 */
	/** shifts rows further back in the matrix
	 * @param m the row starting index
	 * @param s the number of rows to be shifted
	 */
	tMat& rShift(const size_t m, const size_t s) noexcept;
	/** shifts columns further back in the matrix
	 * @param n the column starting index
	 * @param s the number of columns to be shifted
	 */
	tMat& cShift(const size_t n, const size_t s) noexcept;
	/** inserts rows into the matrix
	 * @param m the row index where insertion takes place
	 * @param inp the array to be inserted
	 */
	tMat& rInsert(const size_t m, const tArray& inp) noexcept;
	/** inserts rows into the matrix
	 * @param m the row index where insertion takes place
	 * @param inp the matrix to be inserted
	 */
	tMat& rInsert(const size_t m, const tMat& inp) noexcept;
	/** inserts rows into the matrix
	 * @param itr the row const_iterator pointing to where insertion takes place
	 * @param inp the array to be inserted
	 */
	inline tMat& rInsert(const c_tRowItr& itr, const tArray& inp) noexcept {
		return rInsert((size_t)itr,inp);
	}
	/** inserts rows into the matrix
	 * @param itr the row iterator pointing to where insertion takes place
	 * @param inp the array to be inserted
	 */
	inline tMat& rInsert(const tRowItr& itr, const tArray& inp) noexcept {
		return rInsert((size_t)itr,inp);
	}
	/** inserts rows into the matrix
	 * @param itr the row const_iterator pointing to where insertion takes place
	 * @param inp the matrix to be inserted
	 */
	inline tMat& rInsert(const c_tRowItr& itr, const tMat& inp) noexcept {
		return rInsert((size_t)itr,inp);
	}
	/** inserts rows into the matrix
	 * @param itr the row iterator pointing to where insertion takes place
	 * @param inp the matrix to be inserted
	 */
	inline tMat& rInsert(const tRowItr& itr, const tMat& inp) noexcept {
		return rInsert((size_t)itr,inp);
	}
	/** inserts columns into the matrix
	 * @param n the column index where insertion takes place
	 * @param inp the array to be inserted
	 */
	tMat& cInsert(const size_t n, const tArray& inp) noexcept;
	/** inserts columns into the matrix
	 * @param n the column index where insertion takes place
	 * @param inp the matrix to be inserted
	 */
	tMat& cInsert(const size_t n, const tMat& inp) noexcept;
	/** inserts columns into the matrix
	 * @param itr the column const_iterator pointing to where insertion takes place
	 * @param inp the array to be inserted
	 */
	inline tMat& cInsert(const c_tColItr& itr, const tArray& inp) noexcept {
		return cInsert((size_t)itr,inp);
	}
	/** inserts columns into the matrix
	 * @param itr the column iterator pointing to where insertion takes place
	 * @param inp the array to be inserted
	 */
	inline tMat& cInsert(const tColItr& itr, const tArray& inp) noexcept {
		return cInsert((size_t)itr,inp);
	}
	/** inserts columns into the matrix
	 * @param itr the column const_iterator pointing to where insertion takes place
	 * @param inp the matrix to be inserted
	 */
	inline tMat& cInsert(const c_tColItr& itr, const tMat& inp) noexcept {
		return cInsert((size_t)itr,inp);
	}
	/** inserts columns into the matrix
	 * @param itr the column iterator pointing to where insertion takes place
	 * @param inp the matrix to be inserted
	 */
	inline tMat& cInsert(const tColItr& itr, const tMat& inp) noexcept {
		return cInsert((size_t)itr,inp);
	}
	/** generic set method for inserting data into the matrix
	 * @param inp the matrix holding the data to be inserted
	 * @param m a vector holding row indices
	 * @param n a vector holding column indices
	 */
	tMat& set(const tMat& inp, const std::vector<size_t>& m, const std::vector<size_t>& n) noexcept;
	/** generic set method for inserting data into the matrix
	 * @param inp the matrix holding the data to be inserted
	 * @param m the row index
	 * @param n the column index
	 */
	tMat& set(const tMat& inp, const size_t m, const size_t n) noexcept;
	/** generic set method for inserting data into the matrix using logical indexing
	 * @param inp the matrix holding the data to be inserted
	 * @param m a logical array for the row indices
	 * @param n a logical array for the column indices
	 */
	tMat& setl(const tMat& inp, const fArray& m, const fArray& n) noexcept;
	/** removes a row from the matrix
	 * @param m the row index
	 */
	tMat& rRm(const size_t m) noexcept;
	/** removes a row from the matrix
	 * @param itr row const_iterator
	 */
	inline tMat& rRm(const c_tRowItr& itr) noexcept { return rRm((size_t)itr); }
	/** removes a row from the matrix
	 * @param itr row iterator
	 */
	inline tMat& rRm(const tRowItr& itr) noexcept { return rRm((size_t)itr); }
	/** removes multiple rows from the matrix
	 * @param m a vector of row indices
	 */
	tMat& rRm(const std::vector<size_t>& m) noexcept;
	/** removes a column from the matrix
	 * @param n the column index
	 */
	tMat& cRm(const size_t n) noexcept;
	/** removes a column from the matrix
	 * @param itr column const_iterator
	 */
	inline tMat& cRm(const c_tColItr& itr) noexcept { return cRm((size_t)itr); }
	/** removes a column from the matrix
	 * @param itr column iterator
	 */
	inline tMat& cRm(const tColItr& itr) noexcept { return cRm((size_t)itr); }
	/** removes multiple columns from the matrix
	 * @param n a vector of column indices
	 */
	tMat& cRm(const std::vector<size_t>& n) noexcept;
	/** removes the diagonal from the matrix. The resulting matrix can be assembled in two ways: \n
	 * \verbatim
	   a b c      d b c            a b c       b c
	   d e f  ->  g h j     or     d e f  ->   d f
	   g h j                       g h j       g h
	   for M>N only the first version is possible
	   for N>M only the second version is possible
	   for N=M both are possible                       \endverbatim
	 * @param low switch as to how to assemble the matrix left(false) or right(true)
	 */
	tMat& dRm(const bool low=true) noexcept;
	//! inverts the matrix
	tMat& inv() noexcept;
	

	/** @name erase overload for compatibility with std::vector
	 */
	//! erase a row using a row const_iterator
	inline void erase(const c_tRowItr& itr) noexcept { rRm(itr); }
	//! erase a row using a row iterator
	inline void erase(const tRowItr& itr) noexcept { rRm(itr); }
	//! erase a column using a column const_iterator
	inline void erase(const c_tColItr& itr) noexcept { cRm(itr); }
	//! erase a column using a column itertator
	inline void erase(const tColItr& itr) noexcept { cRm(itr); }
	//! erase a column using a column index
	inline void erase(const size_t n) noexcept { cRm(n); }
	//! erase a single element for row or column shaped matrices using a const_iterator
	inline void erase(const c_tItr& itr) noexcept {
		assert(M()==1 || N()==1);
		if (M()==1) { cRm(distance(this->cbegin(),itr)); return; }
		if (N()==1) { rRm(distance(this->cbegin(),itr)); return; }
	}
	//! erase a single element for row or column shaped matrices using a iterator
	inline void erase(const tItr& itr) noexcept {
		assert(M()==1 || N()==1);
		if (M()==1) { cRm(distance(this->cbegin(),itr)); return; }
		if (N()==1) { rRm(distance(this->cbegin(),itr)); return; }
	}
	

	/** @name conversion
	 */
	//! returns this as a copy
	inline tMat copy() const noexcept { return *this; }
	//! returns this as a real copy
	fMat fcopy() const noexcept;
	//! returns this as a complex copy
	cMat ccopy() const noexcept;
	/** generic get method for extracting a submatrix
	 * @param m a vector of row indices
	 * @param n a vector of column indices
	 */
	tMat get(const std::vector<size_t>& m, const std::vector<size_t>& n) const noexcept;
	/** generic get method for extracting a submatrix
	 * @param m a row index
	 * @param n a column index
	 * @param lm the number of rows to extract
	 * @param ln the number of columns to extract
	 */
	tMat get(const size_t m, const size_t n, const size_t lm=NPOS__, const size_t ln=NPOS__) const noexcept;
	/** generic get method for extracting a submatrix using logical indexing
	 * @param m a logical array for rows
	 * @param n a logical array for columns
	 */
	tMat getl(const fArray& m, const fArray& n) const noexcept;
	/** generic method for extracting a elements as a linear column shaped matrix using logical indexing
	 * @param inp a logical matrix of the same size as this
	 */
	tMat getl(const fMat& inp) const noexcept;
	/** get method for extracting a submatrix of full rows
	 * @param m a row index
	 * @param l the number of rows
	 */
	tMat rGet(const size_t m, const size_t l=1) const noexcept;
	/** get method for extracting a submatrix of full columns
	 * @param n a row index
	 * @param l the number of columns
	 */
	tMat cGet(const size_t n, const size_t l=1) const noexcept;
	/** get method for extracting this without one specific row
	 * @param m a row index
	 */
	tMat rWOGet(const size_t m) const noexcept;
	/** get method for extracting this without one specific column
	 * @param n a column index
	 */
	tMat cWOGet(const size_t n) const noexcept;
	//! returns the elements above the diagonal in a column shaped matrix
	tMat upper() const noexcept;
	//! returns the elements below the diagonal in a column shaped matrix
	tMat lower() const noexcept;


	/** @name logical
	 */
	//! check whether this is a permutation matrix
	bool permutation() const noexcept;
	

	/** @name matrix arithmetic
	 */
	//! left division equals version, i.e. b = A \ b i.e. b = A^-1 * b for real rows as b
	inline fRow& leftDivideEq(fRow& inp) const noexcept {
		assert(this->L()==1);
		return inp/=*data();
	}
	//! left division equals version, i.e. b = A \ b i.e. b = A^-1 * b for complex rows as b
	inline cRow& leftDivideEq(cRow& inp) const noexcept {
		assert(this->L()==1);
		return inp/=*data();
	}
	//! left division equals version, i.e. b = A \ b i.e. b = A^-1 * b for real columns as b
	fCol& leftDivideEq(fCol& inp) const noexcept;
	//! left division equals version, i.e. b = A \ b i.e. b = A^-1 * b for complex columns as b
	cCol& leftDivideEq(cCol& inp) const noexcept;
	//! left division equals version, i.e. B = A \ B i.e. B = A^-1 * B for real matrices as B
	fMat& leftDivideEq(fMat& inp) const noexcept;
	//! left division equals version, i.e. B = A \ B i.e. B = A^-1 * B for complex matrices as B
	cMat& leftDivideEq(cMat& inp) const noexcept;
	//! left division, i.e. A \ b i.e. A^-1 * b for real arrays as b
	inline tMat leftDivide(const fArray& inp) const noexcept {
		assert(square());
		tMat tmp(inp);
		return leftDivideEq(tmp);
	}
	//! left division, i.e. A \ b i.e. A^-1 * b for complex arrays as b
	inline cMat leftDivide(const cArray& inp) const noexcept {
		assert(square());
		cMat tmp(inp);
		return leftDivideEq(tmp);
	}
	//! matrix product A * b for real rows as b
	inline tMat prod(const fRow& inp) const noexcept {
		assert(N()==1);
		return cFront().prod(inp);
	}
	//! matrix product A * b for complex rows as b
	inline cMat prod(const cRow& inp) const noexcept {
		assert(N()==1);
		return cFront().prod(inp);
	}
	//! matrix product A * b for real arrays as b
	tMat prod(const fArray& inp) const noexcept;
	//! matrix product A * b for complex arrays as b
	cMat prod(const cArray& inp) const noexcept;


	/** @name printing
	 */
	//! write this to binary file
	void writeToFile(const std::string& fileName, const bool noheader=false) const;

private:
	// member variables
	size_t M_;			//!< number of rows
	size_t N_;			//!< number of columns
	size_t C_;			//!< linear capacity >= M_*N_
	TT* data_;			//!< data array, size C_

private:
	// member functions
	void readText_(const std::string& fileName);
	void readBinary_(std::ifstream& file, const bool cpx);
	void parse_(std::ifstream& file);
	void readData_(std::ifstream& file, const std::function<TT(std::string)> pnfw);
};

#endif // _LM_TMAT_

/** @}
 */
