// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup libmat
 * @{
 */

#ifndef _LM_OPS_
#define _LM_OPS_

#include "lm_defs.h"
#include "lm_cpxItr.h"
#include "lm_types.h"
#include <iostream>
#include <stack>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace std {
	//! nan for complex
	inline bool isnan(const CPX__& arg) noexcept {
		return std::isnan(std::real(arg)) || std::isnan(std::imag(arg));
	}
}

namespace lm__ {
#ifndef NVARTOL__
	/**
	 * Variable tolerance implemented strough stack.
	 * Instanced in file 'lm_fn.h'. Thread safety demands multiple stacks in parallel regions.
	 */
	namespace hide {
		#ifdef _OPENMP
		extern std::vector<std::stack<RE__>> mtol_st_;		//!< stacks of stored mtol_ values
		#else
		extern std::stack<RE__> mtol_st_;			//!< stack of stored mtol_ values
		#endif
	}
	//! returns size of the (master) tolerance stack
	inline size_t mtol_depth() noexcept {
		#ifdef _OPENMP
		return hide::mtol_st_.front().size();			// size of master stack
		#else
		return hide::mtol_st_.size();				// size of stack
		#endif
	}
	//! returns the current number of tolerance stacks
	inline size_t mtol_nstacks() noexcept {
		#ifdef _OPENMP
		return hide::mtol_st_.size();				// current number of stacks
		#else
		return 1;						// only one stack
		#endif
	}
	//! returns the top of the tolerance stack
	inline RE__ mtol() noexcept {
		#ifdef _OPENMP
		return hide::mtol_st_[omp_get_thread_num()].top();	// or top of stack
		#else
		return hide::mtol_st_.top();				// or top of stack
		#endif
	}
	//! pushes tolerance onto the tolerance stack
	inline void set_mtol(const double tol) noexcept {
		#ifdef _OPENMP
		hide::mtol_st_[omp_get_thread_num()].push(tol);		// push tol to top of nth stack
		#else
		hide::mtol_st_.push(tol);				// push tol to top of the stack
		#endif
	}
	//! pops the top of the tolerance stack
	inline void reset_mtol() noexcept {
		#ifdef _OPENMP
		if (mtol_depth()>1)
			hide::mtol_st_[omp_get_thread_num()].pop();	// remove top element
		#else
		if (mtol_depth()>1)
			hide::mtol_st_.pop();				// remove top element
		#endif
	}
	/** spawns multiple tolerance stacks.
	 * @param N the number of tolerance stacks.
	 */
	inline void spawn_mtol_stacks(const size_t N) noexcept {
		#ifdef _OPENMP
		hide::mtol_st_.resize(N);
		for (auto i=hide::mtol_st_.begin()+1,
			  e=hide::mtol_st_.end(); i!=e; ++i)
			*i = hide::mtol_st_.front();
		#endif
		return;
	}
	//! collapses tolerance stacks back to master stack only
	inline void collapse_mtol_stacks() noexcept {
		#ifdef _OPENMP
		hide::mtol_st_.resize(1);
		#endif
		return;
	}

#else
	//! no variable tolerance for version without tolerance stacks
	namespace hide {
		extern RE__ tol;
	}
	//! returns static tolerance level
	inline RE__ mtol() noexcept {
		return hide::mtol_;
	}
#endif
}


// operators
namespace lm__ {

	/** contains arithmetic operators including a tolerance level. \n
	 * Template arguments are either real or complex. In comparisons the imaginary parts are ignored.
	 */
	namespace ops {
		
		/* strict logical
		 */
		//! strict cast to bool
		template<class T1> bool bool_s(const T1& inp) noexcept {
			return std::real(inp);
		}
		//! strict logical NOT
		template<class T1> bool not_s(const T1& inp) noexcept {
			return !std::real(inp);
		}
		//! strict logical AND
		template<class T1, class T2> bool and_s(const T1& lhs, const T2& rhs) noexcept {
			return std::real(lhs)&&std::real(rhs);
		}
		//! strict logical OR
		template<class T1, class T2> bool or_s(const T1& lhs, const T2& rhs) noexcept {
			return std::real(lhs)||std::real(rhs);
		}


		/* strict comparison
		 */
		//! strict comparison to 0
		template<class T1> bool z_s(const T1& inp) noexcept {
			return not_s(inp);
		}
		//! strict comparison to not 0
		template<class T1> bool nz_s(const T1& inp) noexcept {
			return bool_s(inp);
		}
		//! strict equality
		template<class T1, class T2> bool eq_s(const T1& lhs, const T2& rhs) noexcept {
			return lhs==rhs;
		}
		//! strict inequality
		template<class T1, class T2> bool neq_s(const T1& lhs, const T2& rhs) noexcept {
			return lhs!=rhs;
		}
		//! strict >=
		template<class T1, class T2> bool geq_s(const T1& lhs, const T2& rhs) noexcept {
			return std::real(lhs)>=std::real(rhs);
		}
		//! strict <=
		template<class T1, class T2> bool leq_s(const T1& lhs, const T2& rhs) noexcept {
			return std::real(lhs)<=std::real(rhs);
		}
		//! strict >
		template<class T1, class T2> bool gt_s(const T1& lhs, const T2& rhs) noexcept {
			return std::real(lhs)>std::real(rhs);
		}
		//! strict <
		template<class T1, class T2> bool lt_s(const T1& lhs, const T2& rhs) noexcept {
			return std::real(lhs)<std::real(rhs);
		}


#ifndef NTOLERANT__

		/* tolerant comparison		
		 */
		//! tolerant comparison to 0
		template<class T1> bool z(const T1& inp) noexcept {
			return std::abs(inp)<=mtol();
		}
		//! tolerant comparison to not 0
		template<class T1> bool nz(const T1& inp) noexcept {
			return std::abs(inp)>mtol();
		}
		//! tolerant equality
		template<class T1, class T2> bool eq(const T1& lhs, const T2& rhs) noexcept {
			return std::abs(lhs-rhs)<=mtol();
		}
		//! tolerant inequality
		template<class T1, class T2> bool neq(const T1& lhs, const T2& rhs) noexcept {
			return std::abs(lhs-rhs)>mtol();
		}
		//! tolerant >=
		template<class T1, class T2> bool geq(const T1& lhs, const T2& rhs) noexcept {
			return std::real(rhs)-std::real(lhs)<=mtol();
		}
		//! tolerant <=
		template<class T1, class T2> bool leq(const T1& lhs, const T2& rhs) noexcept {
			return std::real(lhs)-std::real(rhs)<=mtol();
		}
		//! tolerant >
		template<class T1, class T2> bool gt(const T1& lhs, const T2& rhs) noexcept {
			return std::real(lhs)-std::real(rhs)>mtol();
		}
		//! tolerant <
		template<class T1, class T2> bool lt(const T1& lhs, const T2& rhs) noexcept {
			return std::real(rhs)-std::real(lhs)>mtol();
		}
	

		/* tolerant logical
		 */
		//! tolerant cast to bool
		template<class T1> bool bool_(const T1& inp) noexcept {
			return nz(inp);
		}
		//! tolerant logical NOT
		template<class T1> bool not_(const T1& inp) noexcept {
			return z(inp);
		}
		//! tolerant logical AND
		template<class T1, class T2> bool and_(const T1& lhs, const T2& rhs) noexcept {
			return nz(std::real(lhs))&&nz(std::real(rhs));
		}
		//! tolerant logical OR
		template<class T1, class T2> bool or_(const T1& lhs, const T2& rhs) noexcept {
			return nz(std::real(lhs))||nz(std::real(rhs));
		}
#else

		/* intolerant passthrough comparison for version without tolerance
		 */
		//! comparison to 0
		template<class T1> bool z(const T1& inp) noexcept {
			return z_s(inp);
		}
		//! comparison to not 0
		template<class T1> bool nz(const T1& inp) noexcept {
			return nz_s(inp);
		}
		//! equality
		template<class T1, class T2> bool eq(const T1& lhs, const T2& rhs) noexcept {
			return eq_s(lhs,rhs);
		}
		//! inequality
		template<class T1, class T2> bool neq(const T1& lhs, const T2& rhs) noexcept {
			return neq_s(lhs,rhs);
		}
		//! >=
		template<class T1, class T2> bool geq(const T1& lhs, const T2& rhs) noexcept {
			return geq_s(lhs,rhs);
		}
		//! <=
		template<class T1, class T2> bool leq(const T1& lhs, const T2& rhs) noexcept {
			return leq_s(lhs,rhs);
		}
		//! >
		template<class T1, class T2> bool gt(const T1& lhs, const T2& rhs) noexcept {
			return gt_s(lhs,rhs);
		}
		//! <
		template<class T1, class T2> bool lt(const T1& lhs, const T2& rhs) noexcept {
			return lt_s(lhs,rhs);
		}
		

		/* intolerant passthrough logical for version without tolerance
		 */
		//! cast to bool
		template<class T1> bool bool_(const T1& inp) noexcept {
			return bool_s(inp);
		}
		//! logical NOT
		template<class T1> bool not_(const T1& inp) noexcept {
			return not_s(inp);
		}
		//! logical AND
		template<class T1, class T2> bool and_(const T1& lhs, const T2& rhs) noexcept {
			return and_s(lhs,rhs);
		}
		//! logical OR
		template<class T1, class T2> bool or_(const T1& lhs, const T2& rhs) noexcept {
			return or_s(lhs,rhs);
		}
#endif
		
		/* assignment operators generalizing between real and complex
		 */
		//! assign real to real
		inline RE__& assign(RE__& lhs, const RE__ rhs) noexcept {
			return lhs=rhs;
		}
		//! assign complex to real
		inline RE__& assign(RE__& lhs, const CPX__& rhs) noexcept {
			return lhs=std::real(rhs);
		}
		//! assign complex to complex
		inline CPX__& assign(CPX__& lhs, const CPX__& rhs) noexcept {
			return lhs=rhs;
		}
		//! assign real to complex
		inline CPX__& assign(CPX__& lhs, const RE__& rhs) noexcept {
			return lhs=rhs;
		}


		/* assignment equals operators
		 */
		//! add assign real to real
		inline RE__& plusEq(RE__& lhs, const RE__ rhs) noexcept {
			return lhs+=rhs;
		}
		//! subtract assign real to real
		inline RE__& minusEq(RE__& lhs, const RE__ rhs) noexcept {
			return lhs-=rhs;
		}
		//! produce assign real to real
		inline RE__& prodEq(RE__& lhs, const RE__ rhs) noexcept {
			return lhs*=rhs;
		}
		//! divide assign real to real
		inline RE__& divEq(RE__& lhs, const RE__ rhs) noexcept {
			return lhs/=rhs;
		}
		//! add assign complex to real
		inline RE__& plusEq(RE__& lhs, const CPX__& rhs) noexcept {
			return lhs+=std::real(rhs);
		}
		//! subtract assign complex to real
		inline RE__& minusEq(RE__& lhs, const CPX__& rhs) noexcept {
			return lhs-=std::real(rhs);
		}
		//! produce assign complex to real
		inline RE__& prodEq(RE__& lhs, const CPX__& rhs) noexcept {
			return lhs*=std::real(rhs);
		}
		//! divide assign complex to real
		inline RE__& divEq(RE__& lhs, const CPX__& rhs) noexcept {
			return lhs*=(std::real(rhs)/
				(std::real(rhs)*std::real(rhs)+std::imag(rhs)*std::imag(rhs)));
		}
		//! add assign complex to complex
		inline CPX__& plusEq(CPX__& lhs, const CPX__& rhs) noexcept {
			return lhs+=rhs;
		}
		//! add subtract complex to complex
		inline CPX__& minusEq(CPX__& lhs, const CPX__& rhs) noexcept {
			return lhs-=rhs;
		}
		//! add produce complex to complex
		inline CPX__& prodEq(CPX__& lhs, const CPX__& rhs) noexcept {
			return lhs*=rhs;
		}
		//! add divide complex to complex
		inline CPX__& divEq(CPX__& lhs, const CPX__& rhs) noexcept {
			return lhs/=rhs;
		}
		//! add assign real to complex
		inline CPX__& plusEq(CPX__& lhs, const RE__ rhs) noexcept {
			return lhs+=rhs;
		}
		//! subtract assign real to complex
		inline CPX__& minusEq(CPX__& lhs, const RE__ rhs) noexcept {
			return lhs-=rhs;
		}
		//! produce assign real to complex
		inline CPX__& prodEq(CPX__& lhs, const RE__ rhs) noexcept {
			return lhs*=rhs;
		}
		//! divide assign real to complex
		inline CPX__& divEq(CPX__& lhs, const RE__ rhs) noexcept {
			return lhs/=rhs;
		}


		/* ceil, floor, ...
		 */
#ifndef NTOLERANT__
		//! tolerant real ceil function
		inline RE__ ceil(const RE__ inp) noexcept {
			const RE__ tmp=std::floor(inp);	
			return (inp-tmp)<mtol() ? tmp: std::ceil(inp);
		}
		//! tolerant real floor function
		inline RE__ floor(const RE__ inp) noexcept {
			const RE__ tmp=std::ceil(inp);
			return (tmp-inp)<mtol() ? tmp: std::floor(inp);
		}
		//! tolerant real sign function
		inline RE__ sign(const RE__ inp) noexcept {
			return (inp>mtol())-(inp<-mtol());
		}
		//! tolerant complex sign function
		inline CPX__ sign(const CPX__& inp) noexcept {
			return std::abs(inp)<=mtol() ? CPX__(0.0): inp/std::abs(inp);
		}
		//! tolerant real modulus assign function
		inline RE__& modEq(RE__& lhs, const RE__ rhs) {
			if (eq(lhs,RE__(0.0)),eq(rhs,RE__(0.0))) return lhs;
			
			lhs -= lm__::ops::floor(lhs/rhs)*rhs;
			return lhs;
		}
#else
		//! real ceil function
		inline RE__ ceil(const RE__ inp) noexcept {
			return std::ceil(inp);
		}
		//! real floor function
		inline RE__ floor(const RE__ inp) noexcept {
			return std::floor(inp);
		}
		//! real sign function
		inline RE__ sign(const RE__ inp) noexcept {
			return (inp>RE__(0.0))-(inp<RE__(0.0));
		}
		//! complex sign function
		inline CPX__ sign(const CPX__& inp) noexcept {
			return std::abs(inp)<=RE__(0.0) ? CPX__(0.0): inp/std::abs(inp);
		}
		//! real modulus assign function
		inline RE__& modEq(RE__& lhs, const RE__ rhs) {
			if (eq(rhs,RE__(0.0))) return lhs;
			lhs -= std::floor(lhs/rhs)*rhs;
			return lhs;
		}
#endif
		//! (tolerant) complex ceil function
		inline CPX__ ceil(const CPX__& inp) noexcept {
			return {ceil(std::real(inp)),ceil(std::imag(inp))};
		}
		//! (tolerant) complex floor function
		inline CPX__ floor(const CPX__& inp) noexcept {
			return {floor(std::real(inp)),floor(std::imag(inp))};
		}
		//! (tolerant) complex round function
		inline RE__ round(const RE__ inp) noexcept {
			return std::round(inp);
		}
		//! (tolerant) complex round function
		inline CPX__ round(const CPX__& inp) noexcept {
			return {round(std::real(inp)),round(std::imag(inp))};
		}
		//! (tolerant) real and complex modulus assign function
		inline RE__& modEq(RE__& lhs, const CPX__& rhs) {
			return modEq(lhs,std::real(rhs));
		}
		//! (tolerant) complex and complex modulus assign function
		inline CPX__& modEq(CPX__& lhs, const CPX__& rhs) {
			modEq(*reItr(&lhs),std::real(rhs));
			modEq(*imItr(&lhs),std::imag(rhs));
			return lhs;
		}
		//! (tolerant) complex and real modulus assign function
		inline CPX__& modEq(CPX__& lhs, const RE__ rhs) {
			modEq(*reItr(&lhs),rhs);
			modEq(*imItr(&lhs),rhs);
			return lhs;
		}
	}
}

#endif // _OPS_

/** @}
 */
