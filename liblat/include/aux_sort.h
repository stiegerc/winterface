// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _AUX_SORT_SEARCH_
#define _AUX_SORT_SEARCH_

#include <vector>
#include <cstddef>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <cassert>

#ifndef NPOS__
#define NPOS__ std::string::npos	//!< npos value
#endif


namespace aux {

	/** Function to find the sorted order of a range.
	 * @param b iterator to the beginning of the range
	 * @param e iterator past the end of the range
	 * @param cmp compare lambda implementing operator<
	 */
	template<class IT, class CF>
	std::vector<size_t> sorted_order(IT b, const IT e, CF cmp) noexcept {
		using std::distance;
		
		std::vector<size_t> tmp(distance(b,e)); std::iota(tmp.begin(),tmp.end(),size_t(0));
		std::stable_sort(tmp.begin(),tmp.end(),
			[&b,cmp](const size_t i1, const size_t i2) -> bool {return cmp(b[i1],b[i2]);});
		
		std::vector<size_t> res(tmp.size());
		for (size_t i=0; i!=res.size(); ++i)
			res[tmp[i]] = i;
		return res;
	}
	//! overload using standard operator<
	template<class IT>
	std::vector<size_t> sorted_order(IT b, const IT e) noexcept {
		return sorted_order(std::move(b),std::move(e),
			[](const auto& i1, const auto& i2) -> bool {return i1<i2;});
	}
	

	/** Function to reorder a range provided an order vector
	 * @param b iterator to the beginning of the range
	 * @param I the sorted order vector as produced by sorted_order
	 */
	template<class IT>
	void reorder(IT b, std::vector<size_t> I) noexcept {
		assert(!I.empty());
		
		size_t r=I.size()-1;
		for (size_t i=0; i!=I.size() && r; ++i) {
			size_t pos=i;
			while (I[pos]!=NPOS__ && r) {
				std::iter_swap(b+i,b+I[pos]);
				auto tmp=pos;
				pos=I[pos]; I[tmp]=NPOS__; --r;
			}
		}
	}
}

#endif // _AUX_SORT_SEARCH_

/** @}
 */
