// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "lm_tCol.h"
#include "lm_tRow.h"
#include "lm_tMat.h"
#include "lm_ops.h"
#include <fstream>

using namespace lm__;


// copy constructor
template<class TT, class FT, class CT>
lm_tCol<TT,FT,CT>::lm_tCol(const tCol& inp) noexcept: lm_tCol(new tMat(inp),PTRDIFF_MAX) {}


// assignment
template<class TT, class FT, class CT>
lm_tCol<TT,FT,CT>& lm_tCol<TT,FT,CT>::operator=(const FT rhs) noexcept {
	std::for_each(this->begin(),this->end(),[rhs](TT& i){ops::assign(i,rhs);});
	return *this;
}
template<class TT, class FT, class CT>
lm_tCol<TT,FT,CT>& lm_tCol<TT,FT,CT>::operator=(const CT& rhs) noexcept {
	std::for_each(this->begin(),this->end(),[&rhs](TT& i){ops::assign(i,rhs);});
	return *this;
}
template<class TT, class FT, class CT>
lm_tCol<TT,FT,CT>& lm_tCol<TT,FT,CT>::operator=(const fRow& rhs) noexcept {
	assert(this->L()==rhs.L());
	auto j=rhs.begin();
	std::for_each(this->begin(),this->end(),[&j](TT& i){ops::assign(i,*j++);});
	return *this;
}
template<class TT, class FT, class CT>
lm_tCol<TT,FT,CT>& lm_tCol<TT,FT,CT>::operator=(const cRow& rhs) noexcept {
	assert(this->L()==rhs.L());
	auto j=rhs.begin();
	std::for_each(this->begin(),this->end(),[&j](TT& i){ops::assign(i,*j++);});
	return *this;
}
template<class TT, class FT, class CT>
lm_tCol<TT,FT,CT>& lm_tCol<TT,FT,CT>::operator=(const fArray& rhs) noexcept {
	assert(this->L()==rhs.L());
	auto j=rhs.begin();
	std::for_each(this->begin(),this->end(),[&j](TT& i){ops::assign(i,*j++);});
	return *this;
}
template<class TT, class FT, class CT>
lm_tCol<TT,FT,CT>& lm_tCol<TT,FT,CT>::operator=(const cArray& rhs) noexcept {
	assert(this->L()==rhs.L());
	auto j=rhs.begin();
	std::for_each(this->begin(),this->end(),[&j](TT& i){ops::assign(i,*j++);});
	return *this;
}
template<class TT, class FT, class CT>
lm_tCol<TT,FT,CT>& lm_tCol<TT,FT,CT>::operator=(const tCol& rhs) noexcept {
	assert(this->L()==rhs.L());
	if (this!=&rhs)
		memcpy(this->data(),rhs.data(),rhs.L()*sizeof(TT));
	return *this;
}


// conversion
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tCol<TT,FT,CT>::copy() const noexcept {
	return tMat(*this);
}
template<class TT, class FT, class CT>
lm_tMat<FT,FT,CT> lm_tCol<TT,FT,CT>::fcopy() const noexcept {
	return fMat(*this);
}
template<class TT, class FT, class CT>
lm_tMat<CT,FT,CT> lm_tCol<TT,FT,CT>::ccopy() const noexcept {
	return cMat(*this);
}


// matrix arithmetic
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tCol<TT,FT,CT>::prod(const fArray& inp) const noexcept {
	tMat res(this->M(),inp.N());
	
	auto j = res.begin();
	for (auto ir=inp.cbegin(),lr=inp.cend(); ir!=lr; ++ir)
		for (auto il=this->cbegin(), ll=this->cend(); il!=ll; ++il,++j)
			*j = (*il) * (*ir);
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<CT,FT,CT> lm_tCol<TT,FT,CT>::prod(const cArray& inp) const noexcept {
	cMat res(this->M(),inp.N());
	
	auto j = res.begin();
	for (auto ir=inp.cbegin(),lr=inp.cend(); ir!=lr; ++ir)
		for (auto il=this->cbegin(), ll=this->cend(); il!=ll; ++il,++j)
			*j = (*il) * (*ir);
	return res;
}


// printing
template<class TT, class FT, class CT>
void lm_tCol<TT,FT,CT>::writeToFile(const std::string& fileName, const bool noheader) const {
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
	file.write((char*) data(), this->L()*sizeof(TT));
	file.flush();
	file.close();
}
