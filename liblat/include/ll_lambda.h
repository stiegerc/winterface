// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_LAMBDA_
#define _LL_LAMBDA_

#include <algorithm>
#include <numeric>
#include "lm_lambda.h"
#include "ll_types.h"

namespace ll__ {
	using namespace lm__;

	//! compare hamiltonian element to tolerance
	inline auto cmph = [](const hel& h, const double tol) -> bool {
		return std::abs(std::real(h))>=tol ||
		       std::abs(std::imag(h))>=tol;
	};


	/* volume calculation
	 */
	//! volume of a parallelepiped
	inline auto vol = [](const fMat& M) -> double {
		assert(M.square());
		return M.empty() ? 0.0: std::abs(det(M));
	};
	//! normalized volume of a parallelepiped, i.e. vol^(1/DIM)
	inline auto volnorm = [](const fMat& M) -> double {
		assert(M.square());
		return std::pow(vol(M),1.0/M.M());
	};


	/** generate range, like MATLAB s:d:l
	 */
	inline auto rg = [](size_t s, const size_t d, const size_t l) -> std::vector<size_t> {	
		if (!l) return {};
		std::vector<size_t> res(l); res.front()=s;
		std::generate(res.begin()+1,res.end(),[&s,d]{return s+=d;});
		return res;
	};


	/* rv tools
	 */
	//! find the indices where true, like find in MATLAB
	inline auto inds = [](const auto& inp) -> std::vector<size_t> {
		std::vector<size_t> res; res.reserve(inp.size());
		for (size_t i=0; i!=inp.size(); ++i)
			if (inp[i]) res.push_back(i);
		res.shrink_to_fit();
		return res;
	};
	//! find the indices where false, like find in MATLAB
	inline auto ninds = [](const auto& inp) -> std::vector<size_t> {
		std::vector<size_t> res; res.reserve(inp.size());
		for (size_t i=0; i!=inp.size(); ++i)
			if (!inp[i]) res.push_back(i);
		res.shrink_to_fit();
		return res;
	};
	//! find the number true entries
	inline auto Ntrue = [](const auto& inp) -> size_t {
		size_t res=0;
		for (const auto i: inp)
			if (i) ++res;
		return res;
	};
	//! find the number false entries
	inline auto Nfalse = [](const auto& inp) -> size_t {
		size_t res=0;
		for (const auto i: inp)
			if (!i) ++res;
		return res;
	};
}

#endif // _LL_LAMBDA_

/** @}
 */
