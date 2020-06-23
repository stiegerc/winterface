// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "lm_tArray.h"
#include "lm_tMat.h"
#include "lm_ops.h"
#include "lm_cpxItr.h"
#include <vector>
#include <ctime>
#include <fstream>

// logical
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::operator~() const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[](const TT i){return lm__::ops::not_(i);});
	return res;
}


// comparison to logical matrix
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::eq(const FT rhs) const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[rhs](const TT& i){return lm__::ops::eq(i,rhs);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::eq(const CT& rhs) const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[&rhs](const TT& i){return lm__::ops::eq(i,rhs);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::eq(const fArray& rhs) const noexcept {
	assert(M()==rhs.M());
	assert(N()==rhs.N());
	fMat res(M(),N());
	std::transform(cbegin(),cend(),rhs.cbegin(),res.begin(),
		[](const TT& i, const FT j) -> FT {return lm__::ops::eq(i,j);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::eq(const cArray& rhs) const noexcept {
	assert(M()==rhs.M());
	assert(N()==rhs.N());
	fMat res(M(),N());
	std::transform(cbegin(),cend(),rhs.cbegin(),res.begin(),
		[](const TT& i, const CT& j) -> FT {return lm__::ops::eq(i,j);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::neq(const FT rhs) const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[rhs](const TT& i){return lm__::ops::neq(i,rhs);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::neq(const CT& rhs) const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[&rhs](const TT& i){return lm__::ops::neq(i,rhs);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::neq(const fArray& rhs) const noexcept {
	assert(M()==rhs.M());
	assert(N()==rhs.N());
	fMat res(M(),N());
	std::transform(cbegin(),cend(),rhs.cbegin(),res.begin(),
		[](const TT& i, const FT j) -> FT {return lm__::ops::neq(i,j);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::neq(const cArray& rhs) const noexcept {
	assert(M()==rhs.M());
	assert(N()==rhs.N());
	fMat res(M(),N());
	std::transform(cbegin(),cend(),rhs.cbegin(),res.begin(),
		[](const TT& i, const CT& j) -> FT {return lm__::ops::neq(i,j);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::lt(const FT rhs) const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[rhs](const TT& i){return lm__::ops::lt(i,rhs);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::lt(const CT& rhs) const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[&rhs](const TT& i){return lm__::ops::lt(i,rhs);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::lt(const fArray& rhs) const noexcept {
	assert(M()==rhs.M());
	assert(N()==rhs.N());
	fMat res(M(),N());
	std::transform(cbegin(),cend(),rhs.cbegin(),res.begin(),
		[](const TT& i, const FT j) -> FT {return lm__::ops::lt(i,j);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::lt(const cArray& rhs) const noexcept {
	assert(M()==rhs.M());
	assert(N()==rhs.N());
	fMat res(M(),N());
	std::transform(cbegin(),cend(),rhs.cbegin(),res.begin(),
		[](const TT& i, const CT& j) -> FT {return lm__::ops::lt(i,j);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::leq(const FT rhs) const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[rhs](const TT& i){return lm__::ops::leq(i,rhs);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::leq(const CT& rhs) const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[&rhs](const TT& i){return lm__::ops::leq(i,rhs);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::leq(const fArray& rhs) const noexcept {
	assert(M()==rhs.M());
	assert(N()==rhs.N());
	fMat res(M(),N());
	std::transform(cbegin(),cend(),rhs.cbegin(),res.begin(),
		[](const TT& i, const FT j) -> FT {return lm__::ops::leq(i,j);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::leq(const cArray& rhs) const noexcept {
	assert(M()==rhs.M());
	assert(N()==rhs.N());
	fMat res(M(),N());
	std::transform(cbegin(),cend(),rhs.cbegin(),res.begin(),
		[](const TT& i, const CT& j) -> FT {return lm__::ops::leq(i,j);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::gt(const FT rhs) const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[rhs](const TT& i){return lm__::ops::gt(i,rhs);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::gt(const CT& rhs) const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[&rhs](const TT& i){return lm__::ops::gt(i,rhs);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::gt(const fArray& rhs) const noexcept {
	assert(M()==rhs.M());
	assert(N()==rhs.N());
	fMat res(M(),N());
	std::transform(cbegin(),cend(),rhs.cbegin(),res.begin(),
		[](const TT& i, const FT j) -> FT {return lm__::ops::gt(i,j);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::gt(const cArray& rhs) const noexcept {
	assert(M()==rhs.M());
	assert(N()==rhs.N());
	fMat res(M(),N());
	std::transform(cbegin(),cend(),rhs.cbegin(),res.begin(),
		[](const TT& i, const CT& j) -> FT {return lm__::ops::gt(i,j);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::geq(const FT rhs) const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[rhs](const TT& i){return lm__::ops::geq(i,rhs);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::geq(const CT& rhs) const noexcept {
	fMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[&rhs](const TT& i){return lm__::ops::geq(i,rhs);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::geq(const fArray& rhs) const noexcept {
	assert(M()==rhs.M());
	assert(N()==rhs.N());
	fMat res(M(),N());
	std::transform(cbegin(),cend(),rhs.cbegin(),res.begin(),
		[](const TT& i, const FT j) -> FT {return lm__::ops::geq(i,j);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tArray<TT,FT,CT>::geq(const cArray& rhs) const noexcept {
	assert(M()==rhs.M());
	assert(N()==rhs.N());
	fMat res(M(),N());
	std::transform(cbegin(),cend(),rhs.cbegin(),res.begin(),
		[](const TT& i, const CT& j) -> FT {return lm__::ops::geq(i,j);});
	return res;
}


// elementwise arithmetic
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tArray<TT,FT,CT>::operator-() const noexcept {
	tMat res(M(),N());
	std::transform(cbegin(),cend(),res.begin(),[](const TT& i){return -i;});
	return res;
}


// printing
template<class TT, class FT, class CT>
void lm_tArray<TT,FT,CT>::printToFile(const std::string& fileName, const size_t precision) const {
	
	// open file
	std::ofstream file;
	file.open(fileName);
	if (!file.good())
		throw(std::invalid_argument("open file \'"+fileName+"\' failed"));
	
	// print data to file
	file << print(precision);
	
	file.flush();
	file.close();
}


// operators with numeric type as lhs
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> operator-(const FT lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept {
	lm_tMat<TT,FT,CT> res(msize(rhs));
	std::transform(rhs.cbegin(),rhs.cend(),res.begin(),[lhs](const TT& i){return lhs-i;});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<CT,FT,CT> operator-(const CT& lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept {
	lm_tMat<CT,FT,CT> res(msize(rhs));
	std::transform(rhs.cbegin(),rhs.cend(),res.begin(),[lhs](const TT& i){return lhs-i;});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> operator/(const FT lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept {
	lm_tMat<TT,FT,CT> res(msize(rhs));
	std::transform(rhs.cbegin(),rhs.cend(),res.begin(),[lhs](const TT& i){return lhs/i;});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<CT,FT,CT> operator/(const CT& lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept {
	lm_tMat<CT,FT,CT> res(msize(rhs));
	std::transform(rhs.cbegin(),rhs.cend(),res.begin(),[lhs](const TT& i){return lhs/i;});
	return res;
}
