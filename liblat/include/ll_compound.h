// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_COMPOUND_
#define _LL_COMPOUND_

#include "lm_defs.h"
#include "lm_types.h"
#include "lm_tMat.h"
#include "lm_lambda.h"
#include "aux_sort.h"
#include "aux_io.h"
#include <vector>
#include <numeric>
#include <iomanip>


class test_compound_all;

namespace ll__ {

	/** Wrapper for a const matrix.
	 * The api exposing const access is public.
	 * The api exposing all other access is protected.
	 * This class is meant to be inherited from for classes who are essentially
	 * wrappers of a matrix. Some of the api is directly exposed such that we can
	 * write A.ccBegin() instead of A.mat_.ccBegin() without having to write a pass-
	 * through manually.
	 */
	template<class MT, class AT>
	class mat_cb {
	public:
		/** @name constructors
		 */
		//! default constructor
		inline mat_cb() noexcept: mat_() {}
		//! constructor from matrix
		inline mat_cb(MT mat) noexcept: mat_(std::move(mat)) {}
		//! virtual destructor
		inline virtual ~mat_cb() {}


		/** @name information
		 */
		//! check if the matrix is empty
		virtual inline bool empty() const noexcept { return mat_.empty(); }
		//! dimension of space
		inline size_t dim() const noexcept { return mat_.M(); }
		//! number of columns
		inline size_t N() const noexcept { return mat_.N(); }
		//! capacity for columns
		virtual inline size_t capacity() const noexcept { return mat_.ccap(); }


		/** @name const data access
		 */
		//! data const pointer
		auto data() const noexcept -> const decltype(MT().data()) { return mat_.data(); }
		//! direct column const access
		inline auto cAt(const size_t n) const noexcept -> const decltype(MT().cAt(0)) {
			assert(n<N());
			return mat_.cAt(n);
		}
		//! direct const access to the first column
		inline auto cFront() const noexcept -> const decltype(MT().cFront()) {
			assert(!mat_.empty());
			return mat_.cFront();
		}
		//! direct const access to the last column
		inline auto cBack() const noexcept -> const decltype(MT().cBack()) {
			assert(!mat_.empty());
			return mat_.cBack();
		}
	

		/** @name const_iterators
		 */
		//! const_iterator to the beginning of the range
		inline auto cBegin() const noexcept { return mat_.cBegin(); }
		//! const_iterator past the end of the range
		inline auto cEnd() const noexcept { return mat_.cEnd(); }
		//! const_iterator to the beginning of the range
		inline auto ccBegin() const noexcept { return cBegin(); }
		//! const_iterator past the end of the range
		inline auto ccEnd() const noexcept { return cEnd(); }


	protected:
		/** @name iterators
		 */
		//! iterator to the beginning of the range
		inline auto cBegin() noexcept { return this->mat_.cBegin(); }
		//! iterator past the end of the range
		inline auto cEnd() noexcept { return this->mat_.cEnd(); }
		

		/** @name memory management
		 */
		/** reserve space for columns
		 * @param nN new number of columns
		 */
		virtual inline void reserve(const size_t nN) noexcept { this->mat_.reserve(nN); }
		//! shrink allocated memory to the actual number of columns
		virtual inline void shrink_to_fit() noexcept { this->mat_.shrink_to_fit(); }


		/** @name basic modification
		 */
		/** resizes the matrix
		 * @param nN new number of columns
		 */
		virtual inline void resize(const size_t nN) noexcept { this->mat_.resize(nN); }
		/** push_back operator analogous to std::vector
		 * @param inp array to be inserted at the end of the range
		 */
		inline void push_back(const AT& inp) noexcept { this->mat_.push_back(inp); }
		//! pop_back operator analogous to std::vector
		virtual inline void pop_back() noexcept { this->mat_.pop_back(); }
		/** erase a column using a column index
		 * @param n column index
		 */
		virtual inline void erase(const size_t n) noexcept { this->mat_.cRm(n); }

		/** @name friends
		 */
		//! friend test class
		friend class ::test_compound_all;


	protected:
		/** @name member variables
		 */
		MT mat_;	//!< matrix container
	};


	/** Wrapper for a matrix.
	 * The api exposing const access is public.
	 * The api exposing all other access is public.
	 * This class is meant to be inherited from for classes who are essentially
	 * wrappers of a matrix. Some of the api is directly exposed such that we can
	 * write A.ccBegin() instead of A.mat_.ccBegin() without having to write a pass-
	 * through manually.
	 */
	template<class MT, class AT>
	class mat_b {
	public:
		/** @name constructors
		 */
		//! default constructor
		inline mat_b() noexcept: mat_() {}
		//! constructor from matrix
		inline mat_b(MT mat) noexcept: mat_(std::move(mat)) {}
		//! virtual destructor
		inline virtual ~mat_b() {}


		/** @name information
		 */
		//! check if the matrix is empty
		virtual inline bool empty() const noexcept { return mat_.empty(); }
		//! dimension of space
		inline size_t dim() const noexcept { return mat_.M(); }
		//! number of columns
		inline size_t N() const noexcept { return mat_.N(); }
		//! capacity for columns
		virtual inline size_t capacity() const noexcept { return mat_.ccap(); }
		

		/** @name const data access
		 */
		//! data const pointer
		auto data() const noexcept -> const decltype(MT().data()) { return mat_.data(); }
		//! direct column const access
		inline auto cAt(const size_t n) const noexcept -> const decltype(MT().cAt(0)) {
			assert(n<N());
			return mat_.cAt(n);
		}
		//! direct const access to the first column
		inline auto cFront() const noexcept -> const decltype(MT().cFront()) {
			assert(!mat_.empty());
			return mat_.cFront();
		}
		//! direct const access to the last column
		inline auto cBack() const noexcept -> const decltype(MT().cBack()) {
			assert(!mat_.empty());
			return mat_.cBack();
		}

	
		/** @name const iterators
		 */
		//! const_iterator to the beginning of the range
		inline auto cBegin() const noexcept { return mat_.cBegin(); }
		//! const_iterator past the end of the range
		inline auto cEnd() const noexcept { return mat_.cEnd(); }
		//! const_iterator to the beginning of the range
		inline auto ccBegin() const noexcept { return cBegin(); }
		//! const_iterator past the end of the range
		inline auto ccEnd() const noexcept { return cEnd(); }


	public:
		/** @name data access
		 */
		//! data pointer
		auto data() noexcept { return this->mat_.data(); }
		//! direct column access
		inline auto cAt(const size_t n) noexcept {
			assert(n<this->N());
			return this->mat_.cAt(n);
		}
		//! direct access to the first column
		inline auto cFront() noexcept {
			assert(!this->mat_.empty());
			return this->mat_.cFront();
		}
		//! direct access to the last column
		inline auto cBack() noexcept {
			assert(!this->mat_.empty());
			return this->mat_.cBack();
		}


		/** @name iterators
		 */
		//! iterator to the beginning of the range
		inline auto cBegin() noexcept { return this->mat_.cBegin(); }
		//! iterator past the end of the range
		inline auto cEnd() noexcept { return this->mat_.cEnd(); }


		/** @name memory management
		 */
		/** reserve space for columns
		 * @param nN new number of columns
		 */
		virtual inline void reserve(const size_t nN) noexcept { this->mat_.reserve(nN); }
		//! shrink allocated memory to the actual number columns
		virtual inline void shrink_to_fit() noexcept { this->mat_.shrink_to_fit(); }


		/** @name basic modification
		 */
		/** resizes the matrix
		 * @param nN new number of columns
		 */
		virtual inline void resize(const size_t nN) noexcept { this->mat_.resize(nN); }
		//! push_back operator analogous to std::vector
		inline void push_back(const AT& inp) noexcept { this->mat_.push_back(inp); }
		//! pop_back operator analogous to std::vector
		virtual inline void pop_back() noexcept { this->mat_.pop_back(); }
		/** erase a column using a column index
		 * @param n column index
		 */
		virtual inline void erase(const size_t n) noexcept { this->mat_.cRm(n); }
	

		/** @name friends
		 */
		//! friend test class
		friend class ::test_compound_all;


	protected:
		/** @name member variables
		 */
		MT mat_;	//!< matrix container
	};



	/** Wrapper for a const vector (std::vector or otherwise).
	 * The api exposing const access is public.
	 * The api exposing all other access is protected.
	 * This class is meant to be inherited from for classes who are essentially
	 * wrappers of a vector. Some of the api is directly exposed such that we can
	 * write A.cbegin() instead of A.vec_.cbegin() without having to write a pass-
	 * through manually.
	 */
	template<class CT, class CB>
	class vec_cb {
	public:
		/** @name constructor
		 */
		//! default constructor
		inline vec_cb() noexcept: vec_() {}
		//! constructor from vector
		inline vec_cb(CT vec) noexcept: vec_(std::move(vec)) {}
		//! virtual destructor
		inline virtual ~vec_cb() {}
		

		/** @name information
		 */
		//! check if the vector is empty
		virtual inline bool empty() const noexcept { return vec_.empty(); }
		//! size of the vector
		inline size_t size() const noexcept { return vec_.size(); }
		//! capacity of the vector
		virtual inline size_t capacity() const noexcept { return vec_.capacity(); }
		

		/** @name const data access
		 */
		//! data const pointer
		auto data() const noexcept { return vec_.data(); }
		//! direct const access to the first element
		inline const CB& front() const noexcept {
			assert(!vec_.empty());
			return vec_.front();
		}
		//! direct const access to the last element
		inline const CB& back() const noexcept {
			assert(!vec_.empty());
			return vec_.back();
		}
		//! linear indexed const access
		inline const CB& operator[](const size_t i) const noexcept {
			assert(i<vec_.size());
			return vec_[i];
		}
	

		/** @name const iterators
		 */
		//! const_iterator to the beginning of the range
		inline auto begin() const noexcept { return vec_.begin(); }
		//! const_iterator past the end of the range
		inline auto end() const noexcept { return vec_.end(); }
		//! const_iterator to the beginning of the range
		inline auto cbegin() const noexcept { return begin(); }
		//! const_iterator past the end of the range
		inline auto cend() const noexcept { return end(); }


	protected:
		/** @name iterators
		 */
		//! iterator to the beginning of the range
		inline auto begin() noexcept { return vec_.begin(); }
		//! iterator past the end of the range
		inline auto end() noexcept { return vec_.end(); }


		/** @name memory management
		 */
		/** reserve space
		 * @param nN new number of elements
		 */
		virtual inline void reserve(const size_t nN) noexcept { this->vec_.reserve(nN); }
		//! shrink allocated memory to the actual number of elements
		virtual inline void shrink_to_fit() noexcept { this->vec_.shrink_to_fit(); }


		/** @name basic modification
		 */
		/** resize the vector
		 * @param nN new number of elements
		 */
		virtual inline void resize(const size_t nN) noexcept { this->vec_.resize(nN); }
		//! push_back operator
		inline void push_back(CB inp) noexcept { this->vec_.push_back(std::move(inp)); }
		//! pop_back operator
		virtual inline void pop_back() noexcept { this->vec_.pop_back(); }
		/** erase an element using an index
		 * @param n index
		 */
		virtual inline void erase(const size_t n) noexcept {
			this->vec_.erase(this->vec_.cbegin()+n);
		}
		

		/** @name friend
		 */
		//! friend test class
		friend class ::test_compound_all;


	protected:
		/** @name member variables
		 */
		CT vec_;	//!< vector container
	};


	/** Wrapper for a vector (std::vector or otherwise).
	 * The api exposing const access is public.
	 * The api exposing all other access is public.
	 * This class is meant to be inherited from for classes who are essentially
	 * wrappers of a vector. Some of the api is directly exposed such that we can
	 * write A.cbegin() instead of A.vec_.cbegin() without having to write a pass-
	 * through manually.
	 */
	template<class CT, class CB>
	class vec_b {
	public:
		/** name constructors
		 */
		//! default constructor
		inline vec_b() noexcept: vec_() {}
		//! constructor from vector
		inline vec_b(CT vec) noexcept: vec_(std::move(vec)) {}
		//! virtual destructor
		inline virtual ~vec_b() {}
		

		/** @name information
		 */
		//! check if the vector is empty
		virtual inline bool empty() const noexcept { return vec_.empty(); }
		//! size of the vector
		inline size_t size() const noexcept { return vec_.size(); }
		//! capacity of the vector
		virtual inline size_t capacity() const noexcept { return vec_.capacity(); }
		

		/** @name const data access
		 */
		//! data const pointer
		auto data() const noexcept { return vec_.data(); }
		//! direct const access to the first element
		inline const CB& front() const noexcept {
			assert(!vec_.empty());
			return vec_.front();
		}
		//! direct const access to the last element
		inline const CB& back() const noexcept {
			assert(!vec_.empty());
			return vec_.back();
		}
		//! linear indexed const access
		inline const CB& operator[](const size_t i) const noexcept {
			assert(i<vec_.size());
			return vec_[i];
		}
	

		/** @name const iterators
		 */
		//! const_iterator to the beginning of the range
		inline auto begin() const noexcept { return vec_.begin(); }
		//! const_iterator past the end of the range
		inline auto end() const noexcept { return vec_.end(); }
		//! const_iterator to the beginning of the range
		inline auto cbegin() const noexcept { return begin(); }
		//! const_iterator past the end of the range
		inline auto cend() const noexcept { return end(); }


	public:
		/** @name data access
		 */
		//! direct access to the first element
		inline CB& front() noexcept {
			assert(!this->empty());
			return this->vec_.front();
		}
		//! direct access to the last element
		inline CB& back() noexcept {
			assert(!this->empty());
			return this->vec_.back();
		}
		//! linear indexed access
		inline CB& operator[](const size_t i) noexcept {
			assert(i<this->size());
			return this->vec_[i];
		}
		//! linear indexed access
		inline CB& operator[](const lm__::c_fColItr& i) noexcept {
			assert(ptrdiff_t(i)>=0 && ptrdiff_t(i)<this->size());
			return (*this)[i.i()];
		}
		

		/** @name iterators
		 */
		//! iterator to the beginning of the range
		inline auto begin() noexcept { return this->vec_.begin(); }
		//! iterator past the end of the range
		inline auto end() noexcept { return this->vec_.end(); }
		

		/** @name memory management
		 */
		/** reserve space
		 * @param nN new number of elements
		 */
		virtual inline void reserve(const size_t nN) noexcept { this->vec_.reserve(nN); }
		//! shrink allocated memory to the actual number of elements
		virtual inline void shrink_to_fit() noexcept { this->vec_.shrink_to_fit(); }


		/** @name basic modification
		 */
		/** resize the vector
		 * @param nN new number of elements
		 */
		virtual inline void resize(const size_t nN) noexcept { this->vec_.resize(nN); }
		//! push_back operator
		inline void push_back(CB inp) noexcept { this->vec_.push_back(std::move(inp)); }
		//! pop_back operator
		virtual inline void pop_back() noexcept { this->vec_.pop_back(); }
		/** erase an element using an index
		 * @param n index
		 */
		virtual inline void erase(const size_t n) noexcept {
			this->vec_.erase(this->vec_.cbegin()+n);
		}

		
		/** @name  friends
		 */
		//! friend test class
		friend class ::test_compound_all;


	protected:
		/** @name member variables
		 */
		CT vec_;	//!< vector container
	};



	/** Wrapper for a class holding a matrix and vector of the same number of columns/elements.
	 * The api exposing const access is public.
	 * The api exposing all other access is protected.
	 * Some of the api is directly exposed such that we can write A.cbegin() instead of A.vec_.cbegin()
	 * without having to write a passthrough manually. Some of the functions are such that the operation
	 * is done on both containers simultaneously. Both containers must be of the same size at all times.
	 */
	template<class MT, class AT, class CT, class CB>
	class mat_vec_cb: public mat_cb<MT,AT>, public vec_cb<CT,CB> {
	public:
		/** @name constuctors
		 */
		//! default constructor
		inline mat_vec_cb() noexcept: mat_cb<MT,AT>(), vec_cb<CT,CB>() {}
		//! constructor from matrix and vector
		inline mat_vec_cb(MT mat, CT vec) noexcept:
			mat_cb<MT,AT>(std::move(mat)), vec_cb<CT,CB>(std::move(vec)) {
			assert(this->N()==this->size());
		}
		//! virtual destructor
		inline virtual ~mat_vec_cb() {}
		

		/** @name information
		 */
		//! check whether the matrix/vector is empty
		virtual inline bool empty() const noexcept { return mat_cb<MT,AT>::empty(); }
		//! capacity of the matrix/vector
		virtual inline size_t capacity() const noexcept { return vec_cb<CT,CB>::capacity(); }
		

		/** data const access
		 */
		//! data const pointer of the matrix
		auto mdata() const noexcept { return this->mat_.data(); }
		//! data const pointer of the vector
		auto vdata() const noexcept { return this->vec_.data(); }


	protected:
		/** @name data access
		 */
		//! data pointer of the matrix
		auto mdata() noexcept { return this->mat_.data(); }
		//! data pointer of the vector
		auto vdata() noexcept { return this->vec_.data(); }
		

		/** @name memory management
		 */
		/** reserve space both containers
		 * @param nN new number of columns/elements
		 */
		virtual inline void reserve(const size_t nN) noexcept {
			mat_cb<MT,AT>::reserve(nN); vec_cb<CT,CB>::reserve(nN);
		}
		//! shrink allocated memory to the actual number of elements
		virtual inline void shrink_to_fit() noexcept {
			mat_cb<MT,AT>::shrink_to_fit(); vec_cb<CT,CB>::shrink_to_fit();
		}


		/** @name basic modification
		 */
		/** resize the matrix/vector
		 * @param nN new number of columns/elements
		 */
		virtual inline void resize(const size_t nN) noexcept {
			mat_cb<MT,AT>::resize(nN); vec_cb<CT,CB>::resize(nN);
		}
		/** push_back operator for both containers
		 * @param minp the new column to be inserted at the back back of the matrix
		 * @param vinp the new element to be inserted at the back of the vector
		 */
		inline void push_back(const AT& minp, CB vinp) noexcept {
			mat_cb<MT,AT>::push_back(minp); vec_cb<CT,CB>::push_back(std::move(vinp));
		}
		//! pop_back function applied to both containers
		virtual inline void pop_back() noexcept {
			mat_cb<MT,AT>::pop_back(); vec_cb<CT,CB>::pop_back();
		}
		/** erase function applied to both containers
		 * @param n (column) index
		 */
		virtual inline void erase(const size_t n) noexcept {
			mat_cb<MT,AT>::erase(n); vec_cb<CT,CB>::erase(n);
		}
		
		/** @name friends
		 */
		//! friend test class
		friend class ::test_compound_all;
	};
	

	/** Wrapper for a class holding a matrix and vector of the same number of columns/elements.
	 * The api exposing const access is public.
	 * The api exposing all other access is public.
	 * Some of the api is directly exposed such that we can write A.cbegin() instead of A.vec_.cbegin()
	 * without having to write a passthrough manually. Some of the functions are such that the operation
	 * is done on both containers simultaneously. Both containers must be of the same size at all times.
	 */
	template<class MT, class AT, class CT, class CB>
	class mat_vec_b: public mat_b<MT,AT>, public vec_b<CT,CB> {
	public:
		/** @name constuctors
		 */
		//! default constructor
		inline mat_vec_b() noexcept: mat_b<MT,AT>(), vec_b<CT,CB>() {}
		//! constructor from matrix and vector
		inline mat_vec_b(MT mat, CT vec) noexcept:
			mat_b<MT,AT>(std::move(mat)), vec_b<CT,CB>(std::move(vec)) {
			assert(this->N()==this->size());
		}
		//! virtual destructor
		inline virtual ~mat_vec_b() {}
		

		/** @name information
		 */
		//! check whether the matrix/vector is empty
		virtual inline bool empty() const noexcept { return mat_b<MT,AT>::empty(); }
		//! capacity of the matrix/vector
		virtual inline size_t capacity() const noexcept { return vec_b<CT,CB>::capacity(); }
		

		/** @name data const access
		 */
		//! data const pointer of the matrix
		auto mdata() const noexcept { return this->mat_.data(); }
		//! data const pointer of the vector
		auto vdata() const noexcept { return this->vec_.data(); }


	public:
		/** @name data access
		 */
		//! data pointer of the matrix
		auto mdata() noexcept { return this->mat_.data(); }
		//! data pointer of the vector
		auto vdata() noexcept { return this->vec_.data(); }


		/** @name memory management
		 */
		/** reserve space both containers
		 * @param nN new number of columns/elements
		 */
		virtual inline void reserve(const size_t nN) noexcept {
			mat_b<MT,AT>::reserve(nN); vec_b<CT,CB>::reserve(nN);
		}
		//! shrink allocated memory to the actual number of elements
		virtual inline void shrink_to_fit() noexcept {
			mat_b<MT,AT>::shrink_to_fit(); vec_b<CT,CB>::shrink_to_fit();
		}


		/** @name basic modification
		 */
		/** resize the matrix/vector
		 * @param nN new number of columns/elements
		 */
		virtual inline void resize(const size_t nN) noexcept {
			mat_b<MT,AT>::resize(nN); vec_b<CT,CB>::resize(nN);
		}
		/** push_back operator for both containers
		 * @param minp the new column to be inserted at the back back of the matrix
		 * @param vinp the new element to be inserted at the back of the vector
		 */
		inline void push_back(const AT& minp, CB vinp) noexcept {
			mat_b<MT,AT>::push_back(minp); vec_b<CT,CB>::push_back(std::move(vinp));
		}
		//! pop_back function applied to both containers
		virtual inline void pop_back() noexcept {
			mat_b<MT,AT>::pop_back(); vec_b<CT,CB>::pop_back();
		}
		/** erase function applied to both containers
		 * @param n (column) index
		 */
		virtual inline void erase(const size_t n) noexcept {
			mat_b<MT,AT>::erase(n); vec_b<CT,CB>::erase(n);
		}

		
		/** @name friend
		 */
		//! friend test class
		friend class ::test_compound_all;
	};
}

#endif // _LL_COMPOUND_

/** @}
 */
