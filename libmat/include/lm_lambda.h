// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup libmat
 * @{
 */

#ifndef _LM_LAMBDA_
#define _LM_LAMBDA_

#include "lm_ops.h"
#include "lm_fn.h"
#include <algorithm>

namespace lm__ {

	//! standard compare lambda for vectors
	inline auto vcmp = [](const auto& v1, const auto& v2) -> bool {
		assert(v1.size()==v2.size());
		auto j=v2.cbegin();
		for (auto i=v1.cbegin(), e=v1.cend(); i!=e; ++i,++j) {
			if (ops::lt(*i,*j)) return true;
			if (ops::gt(*i,*j)) return false;
		}
		return false;
	};
	//! compare lambda for vectors using custom majority order
	inline auto vcmpmaj = [](const auto& v1, const auto& v2,
			const std::vector<size_t>& majority) -> bool {
		assert(v1.size()==v2.size());
		assert(majority.size()==v1.size());
		
		for (auto m=majority.cbegin(),me=majority.cend(); m!=me; ++m) {
			if (ops::lt(v1[*m],v2[*m])) return true;
			if (ops::gt(v1[*m],v2[*m])) return false;
		}
		return false;
	};

	//! add column to each column of matrix lambda
	inline auto cadd = [](auto& M, const auto& c) -> decltype(M) {
		assert(M.M()==c.L());
		std::for_each(M.cBegin(),M.cEnd(),[&c](auto& i)->void{i+=c;});
		return M;
	};
	//! add row to each row of matrix lambda
	inline auto radd = [](auto& M, const auto& r) -> decltype(M) {
		assert(M.N()==r.L());
		std::for_each(M.rBegin(),M.rEnd(),[&r](auto& i)->void{i+=r;});
		return M;
	};
	
	//! produce column with each column of matrix lambda
	inline auto cprod = [](auto& M, const auto& c) -> decltype(M) {
		assert(c.N()==1 && M.M()==c.M());
		std::for_each(M.cBegin(),M.cEnd(),[&c](auto& i)->void{i*=c;});
		return M;
	};
	//! produce row with each row of matrix lambda
	inline auto rprod = [](auto& M, const auto& r) -> decltype(M) {
		assert(r.M()==1 && M.N()==r.N());
		std::for_each(M.rBegin(),M.rEnd(),[&r](auto& i)->void{i*=r;});
		return M;
	};

	//! check of rows are unique in matrix lambda
	inline auto runique = [](auto M) -> bool {
		std::sort(M.rBegin(),M.rEnd(),vcmp);
		return std::unique(M.rBegin(),M.rEnd())==M.rEnd();
	};
	//! check of columns are unique in matrix lambda
	inline auto cunique = [](auto M) -> bool {
		std::sort(M.cBegin(),M.cEnd(),vcmp);
		return std::unique(M.cBegin(),M.cEnd())==M.cEnd();
	};

	//! check if two arrays are parallel, i.e. equal up to a constant lambda
	inline auto parallel = [](const auto& vec1, const auto& vec2) -> auto {
		assert(vec1.L()==vec2.L());

		const auto chi = lm__::normsq(vec1)/lm__::dot(vec1,vec2);
		auto i2=vec2.cbegin();
		for (auto i1=vec1.cbegin(),e1=vec1.cend(); i1!=e1; ++i1, ++i2)
			if (lm__::ops::neq(*i1,chi * *i2)) return decltype(chi)(0.0);
		return chi;
	};

	//! logical matrix to index vector lambda
	inline auto logicalToI = [](const auto& inp) -> std::vector<size_t> {
		std::vector<size_t> res; res.reserve(inp.L());
		for (size_t i=0; i!=inp.L(); ++i)
			if (inp[i]) res.push_back(i);
		return res;
	};
}

#endif // _LM_LAMBDA_

/** @}
 */
