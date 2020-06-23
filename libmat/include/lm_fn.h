// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup libmat
 * @{
 */

#ifndef _LM_FN_
#define _LM_FN_

#include "lm_tMat.h"
#include <vector>

namespace lm__ {

	/* helper structs return types
	 */
	//! matrix and position struct
	template<class TT, class FT, class CT>
	struct mnp {
		lm_tMat<TT,FT,CT> mat;		//!< matrix
		std::vector<size_t> pos;	//!< position
	};
	//! eigenvalues struct for real matrices
	struct fevr {
		fMat Er;	//!< eigenvalues real parts
		fMat Ei;	//!< eigenvalies complex parts
		fMat Vr;	//!< right eigenvectors
	};
	//! eigenvalues struct for complex matrices
	struct cevr {
		cMat E;		//!< complex eigenvalues
		cMat Vr;	//!< right eigenvectors
	};
	//! eigenvalues and left eigenvectors struct for real matrices
	struct fevl {
		fMat Er;	//!< eigenvalues real parts
		fMat Ei;	//!< eigenvalues imaginary parts
		fMat Vl;	//!< left eigenvectors
	};
	//! eigenvalues and left eigenvectors struct for complex matrices
	struct cevl {
		cMat E;		//!< eigenvalues
		cMat Vl;	//!< left eigenvectors
	};
	//! eigenvalues, left and right eigenvectors struct for real matrices
	struct fevrl {
		fMat Er;	//!< eigenvalues real parts
		fMat Ei;	//!< eigenvalues complex parts
		fMat Vr;	//!< right eigenvectors
		fMat Vl;	//!< left eigenvectors
	};
	//! eigenvalues, left and right eigenvectors struct for complex matrices
	struct cevrl {
		cMat E;		//!< eigenvalues
		cMat Vr;	//!< right eigenvectors
		cMat Vl;	//!< left eigenvectors
	};
	//! complex eigenvalue eigenvector
	struct env {
		CPX__ e;	//!< eigenvalue
		cMat v;		//!< eigenvector
	};
	//! real eigenvalues and real eigenvectors for real hermition matrices
	struct fehv {
		fMat E;		//!< real eigenvalues
		fMat V;		//!< real eigenvectors
	};
	//! real eigenvalues and complex eigenvectors for complex hermition matrices
	struct cehv {
		fMat E;		//!< real eigenvalues
		cMat V;		//!< complex eigenvectors
	};
	//! row, column index pair struct
	struct Ipair {
		size_t m;	//!< row index
		size_t n;	//!< column index
	};


	/* math functions
	 */
	//! absolute values for arrays
	template<class TT, class FT, class CT> fMat abs(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! ceil assign for arrays
	template<class TT, class FT, class CT> lm_tArray<TT,FT,CT>& ceilEq(lm_tArray<TT,FT,CT>& inp) noexcept;
	//! floor assign for arrays
	template<class TT, class FT, class CT> lm_tArray<TT,FT,CT>& floorEq(lm_tArray<TT,FT,CT>& inp) noexcept;
	//! round assign for arrays
	template<class TT, class FT, class CT> lm_tArray<TT,FT,CT>& roundEq(lm_tArray<TT,FT,CT>& inp) noexcept;
	//! sign assign for arrays
	template<class TT, class FT, class CT> lm_tArray<TT,FT,CT>& signEq(lm_tArray<TT,FT,CT>& inp) noexcept;
	//! returns ceil for arrays
	template<class TT, class FT, class CT> lm_tMat<TT,FT,CT> ceil(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! returns floor for arrays
	template<class TT, class FT, class CT> lm_tMat<TT,FT,CT> floor(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! returns round for arrays
	template<class TT, class FT, class CT> lm_tMat<TT,FT,CT> round(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! returns sign for arrays
	template<class TT, class FT, class CT> lm_tMat<TT,FT,CT> sign(const lm_tArray<TT,FT,CT>& inp) noexcept;


	/* eigenvalue computation
	 */
	//! eigenvalues from real matrix
	cMat eig(fMat inp) noexcept;
	//! eigenvalues from complex matrix
	cMat eig(cMat inp) noexcept;
	//! eigenvalues and right eigenvectors from real matrix
	fevr eigr(fMat inp) noexcept;
	//! eigenvalues and right eigenvectors from complex matrix
	cevr eigr(cMat inp) noexcept;
	//! eigenvalues and left eigenvectors from real matrix
	fevl eigl(fMat inp) noexcept;
	//! eigenvalues and left eigenvectors from complex matrix
	cevl eigl(cMat inp) noexcept;
	//! eigenvalues and left and right eigenvectors from real matrix
	fevrl eigrl(fMat inp) noexcept;
	//! eigenvalues and left and right eigenvectors from complex matrix
	cevrl eigrl(cMat inp) noexcept;
	//! eigenvalue extraction function
	env eGet(const fMat& Er, const fMat& Ei, const fMat& V, const size_t i) noexcept;
	//! eigenvalues from hermitian real matrix
	fMat eigh(fMat inp) noexcept;
	//! eigenvalues from hermitian complex matrix
	fMat eigh(cMat inp) noexcept;
	//! eigenvalues and eigenvectors from hermitian real matrix
	fehv eighv(fMat inp) noexcept;
	//! eigenvalues and eigenvectors from hermitian complex matrix
	cehv eighv(cMat inp) noexcept;


	/* matrix generation functions
	 */
	/** square matrix filled with 1
	 * @param M number of rows/columns
	 */
	template<class MT>
	MT ones(const size_t M) noexcept;
	/** matrix filled with 1
	 * @param M number of rows
	 * @param N number of columns
	 */
	template<class MT>
	MT ones(const size_t M, const size_t N) noexcept;
	/** matrix filled with 1
	 * @param S size struct
	 */
	template<class MT>
	MT ones(const lm_size& S) noexcept;
	/** square matrix filled with 0
	 * @param M number of rows/columns
	 */
	template<class MT>
	MT zeros(const size_t M) noexcept;
	/** matrix filled with 0
	 * @param M number of rows
	 * @param N number of columns
	 */
	template<class MT>
	MT zeros(const size_t M, const size_t N) noexcept;
	/** matrix filled with 0
	 * @param S size struct
	 */
	template<class MT>
	MT zeros(const lm_size& S) noexcept;
	/** square identity matrix
	 * @param M number of rows/columns
	 */
	template<class MT>
	MT eye(const size_t M) noexcept;
	/** identity matrix
	 * @param M number of rows
	 * @param N number of columns
	 */
	template<class MT>
	MT eye(const size_t M, const size_t N) noexcept;
	/** identity matrix
	 * @param S size struct
	 */
	template<class MT>
	MT eye(const lm_size& S) noexcept;
	/** identity row
	 * @param N length of row
	 * @param i index of 1
	 */
	template<class MT>
	MT rId(const size_t N, const size_t i) noexcept;
	/** identity column
	 * @param M length of row
	 * @param i index of 1
	 */
	template<class MT>
	MT cId(const size_t M, const size_t i) noexcept;
	/** square matrix filled with random numbers
	 * @param M number of rows/columns
	 */
	template<class MT>
	MT rand(const size_t M) noexcept;
	/** matrix filled with random numbers
	 * @param M number of rows
	 * @param N number of columns
	 * @param rand_min lower bound
	 * @param rand_max upper bound
	 */
	template<class MT>
	MT rand(const size_t M, const size_t N, const RE__ rand_min=0.0, const RE__ rand_max=1.0) noexcept;
	/** matrix filled with random numbers
	 * @param S size struct
	 * @param rand_min lower bound
	 * @param rand_max upper bound
	 */
	template<class MT>
	MT rand(const lm_size& S, const RE__ rand_min=0.0, const RE__ rand_max=1.0) noexcept;
	/** matrix filled with random numbers of size(l)/size(u)
	 * @param l lower bound for each element
	 * @param u upper bound for each element
	 */
	template<class MT>
	MT rand(const fMat& l, const fMat& u) noexcept;
	/** square matrix filled with random integers
	 * @param M number of rows/columns
	 */
	template<class MT>
	MT randi(const size_t M) noexcept;
	/** matrix filled with random integers
	 * @param M number of rows
	 * @param N number of columns
	 * @param rand_min lower bound
	 * @param rand_max upper bound
	 */
	template<class MT>
	MT randi(const size_t M, const size_t N, const long rand_min=0, const long rand_max=100) noexcept;
	/** matrix filled with random integers
	 * @param S size struct
	 * @param rand_min lower bound
	 * @param rand_max upper bound
	 */
	template<class MT>
	MT randi(const lm_size& S, const long rand_min=0, const long rand_max=100) noexcept;
	/** matrix filled with random integers of size(l)/size(u)
	 * @param l lower bound for each element
	 * @param u upper bound for each element
	 */
	template<class MT>
	MT randi(const fMat& l, const fMat& u) noexcept;
	//! concatenate along rows [lhs;rhs]
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> mcat(const lm_tMat<TT,FT,CT>& lhs, const lm_tMat<TT,FT,CT>& rhs) noexcept;
	//! concatenate along rows [lhs;rhs]
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> mcat(const lm_tMat<TT,FT,CT>& lhs, const lm_tRow<TT,FT,CT>& rhs) noexcept;
	//! concatenate along rows [lhs;rhs]
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> mcat(const lm_tRow<TT,FT,CT>& lhs, const lm_tMat<TT,FT,CT>& rhs) noexcept;
	//! concatenate along rows [lhs;rhs]
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> mcat(const lm_tRow<TT,FT,CT>& lhs, const lm_tRow<TT,FT,CT>& rhs) noexcept;
	//! concatenate along columns [lhs rhs]
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> ncat(const lm_tArray<TT,FT,CT>& lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept;


	/* matrix orthogonalisation
	 */
	//! Gram-Schmidt orthogonalization process
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT>& gsorth(lm_tMat<TT,FT,CT>& inp) noexcept;
	//! Gram-Schmidt orthogonalization process
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> gsorth(lm_tMat<TT,FT,CT>&& inp) noexcept;
	//! compute orthogonal complement
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> complement(lm_tMat<TT,FT,CT> inp) noexcept;


	/* functionals
	 */
	//! sum functional for arrays
	template<class TT, class FT, class CT> TT sum(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! product functional for arrays
	template<class TT, class FT, class CT> TT prod(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! min functional for arrays
	template<class TT, class FT, class CT> TT min(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! max functional for arrays
	template<class TT, class FT, class CT> TT max(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! mean functional for arrays
	template<class TT, class FT, class CT> TT mean(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! standard deviation functional for arrays
	template<class TT, class FT, class CT> FT stdd(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! standard deviation functional for arrays given a mean value
	template<class TT, class FT, class CT> FT stdd(const lm_tArray<TT,FT,CT>& inp, const TT& m) noexcept;
	//! squared Frobenius norm for real arrays
	RE__ normsq(const fArray& inp) noexcept;
	//! squared Frobenius norm for complex arrays
	RE__ normsq(const cArray& inp) noexcept;
	//! Frobenius norm for arrays
	template<class TT, class FT, class CT> FT norm(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! trace
	template<class TT, class FT, class CT> TT trace(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! determinante
	template<class TT, class FT, class CT> TT det(lm_tMat<TT,FT,CT> inp) noexcept;
	

	/* dot products
	 */
	//! dot product for real array and real array
	RE__ dot(const fArray& lhs, const fArray& rhs) noexcept;
	//! dot product for real array and complex array
	CPX__ dot(const fArray& lhs, const cArray& rhs) noexcept;
	//! dot product for complex array and real array
	CPX__ dot(const cArray& lhs, const fArray& rhs) noexcept;
	//! dot product for complex array and complex array
	CPX__ dot(const cArray& lhs, const cArray& rhs) noexcept;
	//! dot product without complex adjoint for real array and real array
	RE__ dotu(const fArray& lhs, const fArray& rhs) noexcept;
	//! dot product without complex adjoint for real array and complex array
	CPX__ dotu(const fArray& lhs, const cArray& rhs) noexcept;
	//! dot product without complex adjoint for complex array and real array
	CPX__ dotu(const cArray& lhs, const fArray& rhs) noexcept;
	//! dot product without complex adjoint for complex array and complex array
	CPX__ dotu(const cArray& lhs, const cArray& rhs) noexcept;


	/* vector sums of the form lhs = lhs + f*rhs
	 */
	//! vector sum for real row with real row, real factor
	fRow& sum(fRow& lhs, const fRow& rhs, const RE__ f) noexcept;
	//! vector sum for real row with complex row, real factor
	fRow& sum(fRow& lhs, const cRow& rhs, const RE__ f) noexcept;
	//! vector sum for real row with complex row, complex factor
	fRow& sum(fRow& lhs, const cRow& rhs, const CPX__ f) noexcept;
	//! vector sum for complex row with real row, real factor
	cRow& sum(cRow& lhs, const fRow& rhs, const RE__ f) noexcept;
	//! vector sum for complex row with real row, complex factor
	cRow& sum(cRow& lhs, const fRow& rhs, const CPX__ f) noexcept;
	//! vector sum for complex row with complex row, complex factor
	cRow& sum(cRow& lhs, const cRow& rhs, const CPX__ f) noexcept;
	//! vector sum for real column with real column, real factor
	fCol& sum(fCol& lhs, const fCol& rhs, const RE__ f) noexcept;
	//! vector sum for real column with complex column, real factor
	fCol& sum(fCol& lhs, const cCol& rhs, const RE__ f) noexcept;
	//! vector sum for real column with complex column, complex factor
	fCol& sum(fCol& lhs, const cCol& rhs, const CPX__ f) noexcept;
	//! vector sum for complex column with real column, real factor
	cCol& sum(cCol& lhs, const fCol& rhs, const RE__ f) noexcept;
	//! vector sum for complex column with real column, complex factor
	cCol& sum(cCol& lhs, const fCol& rhs, const CPX__ f) noexcept;
	//! vector sum for complex column with complex column, complex factor
	cCol& sum(cCol& lhs, const cCol& rhs, const CPX__ f) noexcept;
	//! vector sum for real row with real row, real factor
	fRow&& sum(fRow&& lhs, const fRow& rhs, const RE__ f) noexcept;
	//! vector sum for real row with complex row, real factor
	fRow&& sum(fRow&& lhs, const cRow& rhs, const RE__ f) noexcept;
	//! vector sum for real row with complex row, complex factor
	fRow&& sum(fRow&& lhs, const cRow& rhs, const CPX__ f) noexcept;
	//! vector sum for complex row with real row, real factor
	cRow&& sum(cRow&& lhs, const fRow& rhs, const RE__ f) noexcept;
	//! vector sum for complex row with real row, complex factor
	cRow&& sum(cRow&& lhs, const fRow& rhs, const CPX__ f) noexcept;
	//! vector sum for complex row with complex row, complex factor
	cRow&& sum(cRow&& lhs, const cRow& rhs, const CPX__ f) noexcept;
	//! vector sum for real column with real column, real factor
	fCol&& sum(fCol&& lhs, const fCol& rhs, const RE__ f) noexcept;
	//! vector sum for real column with complex column, real factor
	fCol&& sum(fCol&& lhs, const cCol& rhs, const RE__ f) noexcept;
	//! vector sum for real column with complex column, complex factor
	fCol&& sum(fCol&& lhs, const cCol& rhs, const CPX__ f) noexcept;
	//! vector sum for complex column with real column, real factor
	cCol&& sum(cCol&& lhs, const fCol& rhs, const RE__ f) noexcept;
	//! vector sum for complex column with real column, complex factor
	cCol&& sum(cCol&& lhs, const fCol& rhs, const CPX__ f) noexcept;
	//! vector sum for complex column with complex column, complex factor
	cCol&& sum(cCol&& lhs, const cCol& rhs, const CPX__ f) noexcept;


	/* R^3 only
	 */
	//! cross product for real array with real array
	fMat cross(const fArray& lhs, const fArray& rhs) noexcept;
	//! normalized cross product for real array with real array
	fMat crossn(const fArray& lhs, const fArray& rhs) noexcept;
	/** rotation matrix around vector
	 * @param phi rotation angle
	 * @param n normalized vector around which to rotate
	 */
	fMat getR(const RE__ phi, const fArray& n) noexcept;
	/** rotation matrix such that one vector becomes another
	 * @param v the initial vector
	 * @param w the target vector
	 */
	fMat getR(const fArray& v, const fArray& w) noexcept;


	/* conversion functions, i.e. functionals along one dimension
	 */
	//! sum of each column
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> msum(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! sum of each row
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> nsum(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! product of each row
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> mprod(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! product of each column
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> nprod(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! find minimum of each column
	template<class TT, class FT, class CT>
	mnp<TT,FT,CT> mmin(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! find minimum of each row
	template<class TT, class FT, class CT>
	mnp<TT,FT,CT> nmin(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! find maximum of each column
	template<class TT, class FT, class CT>
	mnp<TT,FT,CT> mmax(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! find maximum of each row
	template<class TT, class FT, class CT>
	mnp<TT,FT,CT> nmax(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! find mean of each column
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> mmean(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! find mean of each row
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> nmean(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! find squared norm of each column
	template<class TT, class FT, class CT>
	fMat mnormsq(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! find squared norm of each row
	template<class TT, class FT, class CT>
	fMat nnormsq(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! find norm of each column
	template<class TT, class FT, class CT>
	fMat mnorm(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! find norm of each row
	template<class TT, class FT, class CT>
	fMat nnorm(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! find dot product of the columns of a matrix with a real matrix
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> mdot(const lm_tMat<TT,FT,CT>& lhs, const fMat& rhs) noexcept;
	//! find dot product of the columns of a matrix with a complex matrix
	template<class TT, class FT, class CT>
	cMat mdot(const lm_tMat<TT,FT,CT>& lhs, const cMat& rhs) noexcept;
	//! find dot product of the rows of a matrix with a real matrix
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> ndot(const lm_tMat<TT,FT,CT>& lhs, const fMat& rhs) noexcept;
	//! find dot product of the rows of a matrix with a complex matrix
	template<class TT, class FT, class CT>
	cMat ndot(const lm_tMat<TT,FT,CT>& lhs, const cMat& rhs) noexcept;
	//! find dot product without adjoint of the columns of a matrix with a real matrix
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> mdotu(const lm_tMat<TT,FT,CT>& lhs, const fMat& rhs) noexcept;
	//! find dot product without adjoint of the columns of a matrix with a complex matrix
	template<class TT, class FT, class CT>
	cMat mdotu(const lm_tMat<TT,FT,CT>& lhs, const cMat& rhs) noexcept;
	//! find dot product without adjoint of the rows of a matrix with a real matrix
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> ndotu(const lm_tMat<TT,FT,CT>& lhs, const fMat& rhs) noexcept;
	//! find dot product without adjoint of the rows of a matrix with a complex matrix
	template<class TT, class FT, class CT>
	cMat ndotu(const lm_tMat<TT,FT,CT>& lhs, const cMat& rhs) noexcept;
	

	/* matrix conversions
	 */
	//! find the inverse of a matrix
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> inv(lm_tMat<TT,FT,CT> inp) noexcept;
	/** get the diagonal of a matrix or generate a matrix with a vector on the diagonal
	 * @param inp the matrix to extract the diagonal OR the row/col to use as the diagonal in a matrix
	 * @param os offset, m ? to the right: to the left
	 * @param m offset switch
	 */
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> diag(const lm_tArray<TT,FT,CT>& inp, const size_t os=0, const bool m=true) noexcept;
	//! get 'rowized' matrix, i.e. size MxN -> 1xM*N
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> R(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! get 'columnized' matrix, i.e. size MxN -> M*Nx1
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> C(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! get transpose
	template<class TT, class FT, class CT>
	lm_tMat<TT,FT,CT> T(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! get index pairs where entries in an array are 0
	template<class TT, class FT, class CT>
	std::vector<Ipair> Iz(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! get index pairs where entries in an array are not 0
	template<class TT, class FT, class CT>
	std::vector<Ipair> Inz(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! complex conjugate of a real array
	fMat conj(const fArray& inp) noexcept;
	//! complex conjugate of a complex array
	cMat conj(const cArray& inp) noexcept;
	//! real part of a real array
	fMat real(const fArray& inp) noexcept;
	//! real part of a complex array
	fMat real(const cArray& inp) noexcept;
	//! imaginary part of a real array
	fMat imag(const fArray& inp) noexcept;
	//! imaginary part of a complex array
	fMat imag(const cArray& inp) noexcept;


	/* information functions
	 */
	//! check for nan in an array
	template<class TT, class FT, class CT>
	fMat isnan(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! returns max(M,N)
	template<class TT, class FT, class CT>
	size_t length(const lm_tArray<TT,FT,CT>& inp) noexcept;
	/** find the rank of a real matrix
	 * @param inp real matrix
	 * @param min_sv minimum tolerable singular value
	 */
	size_t rank(fMat inp, const RE__ min_sv=MTOL__) noexcept;
	/** find the rank of a complex matrix
	 * @param inp complex matrix
	 * @param min_sv minimum tolerable singular value
	 */
	size_t rank(cMat inp, const RE__ min_sv=MTOL__) noexcept;


	/* comparison functions
	 */
	//! check if all are true for logical arrays
	template<class TT, class FT, class CT>
	bool all(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! check if all in each column are true for logical arrays
	template<class TT, class FT, class CT>
	fMat mall(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! check if all in each row are true for logical arrays
	template<class TT, class FT, class CT>
	fMat nall(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! check if any are true for logical arrays
	template<class TT, class FT, class CT>
	bool any(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! check if any in each column are true for logical arrays
	template<class TT, class FT, class CT>
	fMat many(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! check if any in each row are true for logical arrays
	template<class TT, class FT, class CT>
	fMat nany(const lm_tMat<TT,FT,CT>& inp) noexcept;
	//! find indices to non zero entries
	template<class TT, class FT, class CT>
	std::vector<size_t> find(const lm_tArray<TT,FT,CT>& inp) noexcept;
	//! find the amount of non zero entries
	template<class TT, class FT, class CT>
	size_t nnz(const lm_tArray<TT,FT,CT>& inp) noexcept;
}

#endif // _LM_FN_

/** @}
 */
