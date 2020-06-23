// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup libmat
 * @{
 */

#ifndef _LM_TYPES_
#define _LM_TYPES_

#include "lm_defs.h"
#include <cstddef>
#include <iterator>
#include <vector>

template<class TT>
class lm_c_tItr;
template<class TT>
class lm_tItr;
template<ptrdiff_t>
class lm_c_cpxItr;
template<ptrdiff_t>
class lm_cpxItr;
template<class MT, class VT>
class lm_c_tVecItr;
template<class MT, class VT>
class lm_tVecItr;
template<class MT, class VT>
class lm_cr_tVecItr;
template<class MT, class VT>
class lm_r_tVecItr;
template<class TT, class FT, class CT>
class lm_tArray;
template<class TT, class FT, class CT>
class lm_tRow;
template<class TT, class FT, class CT>
class lm_tCol;
template<class TT, class FT, class CT>
class lm_tMat;

namespace lm__ {

	/* iterator types
	 */
	typedef lm_c_tItr<RE__> c_fItr;			//!< real custom increment const_iterator
	typedef lm_c_tItr<CPX__> c_cItr;		//!< complex custom increment const_iterator
	typedef lm_tItr<RE__> fItr;			//!< real custom increment iterator
	typedef lm_tItr<CPX__> cItr;			//!< complex custom increment iterator
	typedef std::reverse_iterator<c_fItr> cr_fItr;	//!< reverse real custom increment const_iterator
	typedef std::reverse_iterator<c_cItr> cr_cItr;	//!< reverse complex custom increment const_iterator
	typedef std::reverse_iterator<fItr> r_fItr;	//!< reverse real custom increment iterator
	typedef std::reverse_iterator<cItr> r_cItr;	//!< reverse complex custom increment iterator
	typedef lm_c_cpxItr<0> c_reItr;			//!< re only custom increment const_iterator
	typedef lm_c_cpxItr<1> c_imItr;			//!< im only custom increment const_iterator
	typedef lm_cpxItr<0> reItr;			//!< re only custom increment iterator
	typedef lm_cpxItr<1> imItr;			//!< im only custom increment iterator


	/* vector iterator types
	 */
	//! real row const_iterator
	typedef lm_c_tVecItr<lm_tMat<RE__,RE__,CPX__>,lm_tRow<RE__,RE__,CPX__>> c_fRowItr;
	//! complex row const_iterator
	typedef lm_c_tVecItr<lm_tMat<CPX__,RE__,CPX__>,lm_tRow<CPX__,RE__,CPX__>> c_cRowItr;
	//! real column const_iterator
	typedef lm_c_tVecItr<lm_tMat<RE__,RE__,CPX__>,lm_tCol<RE__,RE__,CPX__>> c_fColItr;
	//! complex column const_iterator
	typedef lm_c_tVecItr<lm_tMat<CPX__,RE__,CPX__>,lm_tCol<CPX__,RE__,CPX__>> c_cColItr;
	//! real row iterator
	typedef lm_tVecItr<lm_tMat<RE__,RE__,CPX__>,lm_tRow<RE__,RE__,CPX__>> fRowItr;
	//! complex row iterator
	typedef lm_tVecItr<lm_tMat<CPX__,RE__,CPX__>,lm_tRow<CPX__,RE__,CPX__>> cRowItr;
	//! real column iterator
	typedef lm_tVecItr<lm_tMat<RE__,RE__,CPX__>,lm_tCol<RE__,RE__,CPX__>> fColItr;
	//! complex column iterator
	typedef lm_tVecItr<lm_tMat<CPX__,RE__,CPX__>,lm_tCol<CPX__,RE__,CPX__>> cColItr;
	//! reverse real row const_iterator
	typedef lm_cr_tVecItr<lm_tMat<RE__,RE__,CPX__>,lm_tRow<RE__,RE__,CPX__>> cr_fRowItr;
	//! reverse complex row const_iterator
	typedef lm_cr_tVecItr<lm_tMat<CPX__,RE__,CPX__>,lm_tRow<CPX__,RE__,CPX__>> cr_cRowItr;
	//! reverse real column const_iterator
	typedef lm_cr_tVecItr<lm_tMat<RE__,RE__,CPX__>,lm_tCol<RE__,RE__,CPX__>> cr_fColItr;
	//! reverse complex column const_iterator
	typedef lm_cr_tVecItr<lm_tMat<CPX__,RE__,CPX__>,lm_tCol<CPX__,RE__,CPX__>> cr_cColItr;
	//! reverse real row iterator
	typedef lm_r_tVecItr<lm_tMat<RE__,RE__,CPX__>,lm_tRow<RE__,RE__,CPX__>> r_fRowItr;
	//! reverse complex row iterator
	typedef lm_r_tVecItr<lm_tMat<CPX__,RE__,CPX__>,lm_tRow<CPX__,RE__,CPX__>> r_cRowItr;
	//! reverse real column iterator
	typedef lm_r_tVecItr<lm_tMat<RE__,RE__,CPX__>,lm_tCol<RE__,RE__,CPX__>> r_fColItr;
	//! reverse complex column iterator
	typedef lm_r_tVecItr<lm_tMat<CPX__,RE__,CPX__>,lm_tCol<CPX__,RE__,CPX__>> r_cColItr;


	/* array, row, column and matrix types
	 */
	typedef lm_tArray<RE__,RE__,CPX__> fArray;	//!< real array
	typedef lm_tArray<CPX__,RE__,CPX__> cArray;	//!< complex array
	typedef lm_tRow<RE__,RE__,CPX__> fRow;		//!< real row
	typedef lm_tRow<CPX__,RE__,CPX__> cRow;		//!< complex row
	typedef lm_tCol<RE__,RE__,CPX__> fCol;		//!< real column
	typedef lm_tCol<CPX__,RE__,CPX__> cCol;		//!< complex column
	typedef lm_tMat<RE__,RE__,CPX__> fMat;		//!< real matrix
	typedef lm_tMat<CPX__,RE__,CPX__> cMat;		//!< complex matrix


	//! size struct
	struct lm_size {
		size_t M;	//!< number of rows
		size_t N;	//!< number of columns
		//! comparison operator
		inline bool operator==(const lm_size& rhs) const noexcept {
			return (M==rhs.M) && (N==rhs.N);
		}
		//! comparison operator
		inline bool operator!=(const lm_size& rhs) const noexcept {
			return !(*this==rhs);
		}
	};
	//! streaming operator
	inline std::ostream& operator<<(std::ostream& os, const lm_size& inp) noexcept {
		return (os<<"("<<inp.M<<","<<inp.N<<")");
	}
}

#endif // _LM_TYPES_

/** @}
 */
