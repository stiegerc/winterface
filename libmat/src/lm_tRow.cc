// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "lm_tRow.h"
#include "lm_tMat.h"
#include "lm_fn.h"
#include "lm_ops.h"
#include <fstream>

using namespace lm__;


// copy constructor
template<class TT, class FT, class CT>
lm_tRow<TT,FT,CT>::lm_tRow(const tRow& inp) noexcept: lm_tRow(new tMat(inp),PTRDIFF_MAX) {}


// assignment
template<class TT, class FT, class CT>
lm_tRow<TT,FT,CT>& lm_tRow<TT,FT,CT>::operator=(const FT rhs) noexcept {
	std::for_each(this->begin(),this->end(),[rhs](TT& i){ops::assign(i,rhs);});
	return *this;
}
template<class TT, class FT, class CT>
lm_tRow<TT,FT,CT>& lm_tRow<TT,FT,CT>::operator=(const CT& rhs) noexcept {
	std::for_each(this->begin(),this->end(),[&rhs](TT& i){ops::assign(i,rhs);});
	return *this;
}
template<class TT, class FT, class CT>
lm_tRow<TT,FT,CT>& lm_tRow<TT,FT,CT>::operator=(const lm_tArray<FT,FT,CT>& rhs) noexcept {
	assert(this->L()==rhs.L());
	auto j=rhs.begin();
	std::for_each(this->begin(),this->end(),[&j](TT& i){ops::assign(i,*j++);});
	return *this;
}
template<class TT, class FT, class CT>
lm_tRow<TT,FT,CT>& lm_tRow<TT,FT,CT>::operator=(const lm_tArray<CT,FT,CT>& rhs) noexcept {
	assert(this->L()==rhs.L());
	auto j=rhs.begin();
	std::for_each(this->begin(),this->end(),[&j](TT& i){ops::assign(i,*j++);});
	return *this;
}
template<class TT, class FT, class CT>
lm_tRow<TT,FT,CT>& lm_tRow<TT,FT,CT>::operator=(const lm_tRow<TT,FT,CT>& rhs) noexcept {
	assert(this->L()==rhs.L());
	if (this!=&rhs) {
		auto j=rhs.begin();
		std::for_each(this->begin(),this->end(),[&j](TT& i){ops::assign(i,*j++);});
	}
	return *this;
}


// conversion
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tRow<TT,FT,CT>::copy() const noexcept {
	return tMat(*this);
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tRow<TT,FT,CT>::fcopy() const noexcept {
	return fMat(*this);
}
template<class TT, class FT, class CT>
lm_tMat<CT,FT,CT> lm_tRow<TT,FT,CT>::ccopy() const noexcept {
	return cMat(*this);
}


// matrix arithmetic
template<class TT, class FT, class CT>
TT lm_tRow<TT,FT,CT>::prod(const lm_tCol<FT,FT,CT>& inp) const noexcept {
	assert(this->L()==inp.L());
	return dotu(*this,inp);
}
template<class TT, class FT, class CT>
CT lm_tRow<TT,FT,CT>::prod(const lm_tCol<CT,FT,CT>& inp) const noexcept {
	assert(this->L()==inp.L());
	return dotu(*this,inp);
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tRow<TT,FT,CT>::prod(const fMat& inp) const noexcept {
	assert(this->L()==inp.M());
	tMat res(1,inp.N());
	
	std::transform(inp.ccBegin(),inp.ccEnd(),res.begin(),
		[this](const lm_tCol<FT,FT,CT>& i){return dotu(*this,i);});
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<CT,FT,CT> lm_tRow<TT,FT,CT>::prod(const cMat& inp) const noexcept {
	assert(this->L()==inp.M());
	cMat res(1,inp.N());
	std::transform(inp.ccBegin(),inp.ccEnd(),res.begin(),
		[this](const lm_tCol<CT,FT,CT>& i){return dotu(*this,i);});
	return res;
}


// printing
template<class TT, class FT, class CT>
void lm_tRow<TT,FT,CT>::writeToFile(const std::string& fileName, const bool noheader) const {
	// open file
	std::ofstream file;
	file.open(fileName,std::ios::binary);
	if (!file.good())
		throw(std::invalid_argument("open file \'"+fileName+"\' failed"));

	// write header
	const unsigned long M = this->M(), N = this->N();
	if (!noheader) {
		const std::string hdr = this->cpx() ? "fMat": "cMat";
		file.write(hdr.c_str(),5);
		file.write((char*) &M, sizeof(unsigned long));
		file.write((char*) &N, sizeof(unsigned long));
	}

	// write matrix
	for (auto i=this->cbegin(),e=this->cend(); i!=e; ++i)
		file.write((char*) i.data(), sizeof(TT));
	
	file.flush();
	file.close();
}
