// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "lm_defs.h"
#include "lm_types.h"
#include "lm_cpxItr.h"
#include "blas.h"
#include <cstring>
#include <cassert>
#include <fstream>
#include <regex>

#include "lm_tMat.h"
#include "lm_tMat.cc"


using namespace lm__;

// constructors
template<>
cMat::lm_tMat(const fArray& re, const fArray& im) noexcept: tMat(re.M(),re.N()) {
	assert(re.M()==im.M());
	assert(re.N()==im.N());
	std::copy(re.cbegin(),re.cend(),reItr(this->begin()));
	std::copy(im.cbegin(),im.cend(),imItr(this->begin()));
}
template<>
cMat::lm_tMat(const cArray& inp) noexcept: tMat(inp.M(),inp.N()) {
	if (inp.incr()==1)
		memcpy(data(),inp.data(),this->L()*sizeof(CPX__));
	else {
		auto j=inp.begin();
		std::for_each(this->begin(),this->end(),[&j](CPX__& i){ops::assign(i,*j++);});
	}
}


// assignment
template<>
cMat& cMat::operator=(const cArray& rhs) noexcept {
	if ((void*)this!=(void*)&rhs) {
		resize(rhs.M(),rhs.N());

		if (rhs.incr()==1)
			memcpy(data(),rhs.data(),this->L()*sizeof(CPX__));
		else {
			auto j=rhs.begin();
			std::for_each(this->begin(),this->end(),[&j](CPX__& i){ops::assign(i,*j++);});
		}
	}
	return *this;
}

// conversion
template<>
fMat cMat::fcopy() const noexcept {
	return fMat(*this);
}
template<>
cMat cMat::ccopy() const noexcept {
	return *this;
}


// basic modification
template<>
cMat& cMat::T() noexcept {
	if (this->empty() || M()==1 || N()==1) {
		std::swap(M_,N_);
		for (imItr i=this->begin(),e=this->end(); i!=e; ++i)
			*i = -(*i);
		return *this;
	}	
	if (square()) {
		auto ri = rBegin();
		auto ci = cBegin();
		for (size_t i=1,l=this->M(); i!=l; ++i,++ri,++ci)
			for (auto rj=ri->begin()+i,cj=ci->begin()+i,ce=ci->end(); cj!=ce; ++rj,++cj) {
				const auto tmp = std::conj(*cj);
				*cj = std::conj(*rj); *rj = tmp;
			}
		for (imItr i=this->dbegin(),e=this->dend(); i!=e; ++i)
			*i = -*i;
		return *this;
	}

	cMat tmp(N(),M());
	auto ri = tmp.begin();
	for (auto i=crBegin(),ie=crEnd(); i!=ie; ++i)
		for (auto j=i->begin(),je=i->end(); j!=je; ++j,++ri)
			*ri = std::conj(*j);
	return *this=std::move(tmp);
}


// matrix arithmetic
template<>
fCol& cMat::leftDivideEq(fCol& inp) const noexcept {
	// left matrix division (aka A\B) equals version, equivalent to B=A^-1*B
	return inp = fMat(leftDivide(inp));
}
template<>
cCol& cMat::leftDivideEq(cCol& inp) const noexcept {
	// left matrix division (aka A\B) equals version, equivalent to B=A^-1*B
	assert(square());
	assert(N()==inp.M());

	cMat tmp = *this;
	int* ipiv = new int [M()]; int info;
	c_xgesv(M(),1,tmp.data(),M(),ipiv,inp.data(),M(),&info);
	assert(!info);
	
	delete[] ipiv;
	return inp;
}
template<>
fMat& cMat::leftDivideEq(fMat& inp) const noexcept {
	// left matrix division (aka A\B) equals version, equivalent to B=A^-1*B
	return inp = fMat(leftDivide(inp));
}
template<>
cMat& cMat::leftDivideEq(cMat& inp) const noexcept {
	// left matrix division (aka A\B) equals version, equivalent to B=A^-1*B
	assert(square());
	assert(N()==inp.M());

	cMat tmp = *this;
	int* ipiv = new int [M()]; int info;
	c_xgesv(M(),inp.N(),tmp.data(),M(),ipiv,inp.data(),M(),&info);
	assert(!info);
	
	delete[] ipiv;
	return inp;
}
template<>
cMat cMat::prod(const fArray& inp) const noexcept {
	assert(!empty());
	assert(!inp.empty());
	assert(N()==inp.M());
	
	cMat res(M(),inp.N());

	const RE__* this_ptr = c_reItr(begin()).data();
	const RE__* const this_end=this_ptr+2*M();
	RE__* res_ptr = reItr(res.begin()).data();
	const size_t thisM=2*M(), resM=2*res.M();

	while (this_ptr!=this_end) {
		c_xgemv('T',inp.M(),inp.N(),1.0,inp.data(),inp.M(),this_ptr++,thisM,0.0,res_ptr++,resM);
		c_xgemv('T',inp.M(),inp.N(),1.0,inp.data(),inp.M(),this_ptr++,thisM,0.0,res_ptr++,resM);
	}
	
	return res;
}
template<>
cMat cMat::prod(const cArray& inp) const noexcept {
	assert(!empty());
	assert(!inp.empty());
	assert(N()==inp.M());
	
	cMat res(M(),inp.N());
	c_xgemm('N','N',M(),inp.N(),N(),1.0,data(),M(),inp.data(),N(),0.0,res.data(),M());
	
	return res;
}


// member functions
template<>
void cMat::readBinary_(std::ifstream& file, const bool cpx) {
	uint64_t M,N;
	file.read((char*) &M, sizeof(uint64_t));
	file.read((char*) &N, sizeof(uint64_t));
	(*this) = cMat(M,N);
	
	if (cpx) file.read((char*) data(),L()*sizeof(CPX__));
	else
		for (auto& i: *this) {
			RE__ buff;
			file.read((char*) &buff, sizeof(RE__));
			i = {buff,0.0};
		}
}
template<>
void cMat::parse_(std::ifstream& file) {
	assert(file.good());
	assert(!empty());
	
	// regex'
	const std::regex rn("(?=[-+i])");

	// lambda to parse numbers from string
	auto pnfw = [&rn](const std::string& inp) -> std::complex<RE__> {
		// tokenize number
		std::sregex_token_iterator i(inp.begin(),inp.end(),rn,-1), e;
		if (!i->length()) ++i; // inp starts with delimiter

		// find position of 'i'
		auto ipos = std::find(inp.begin(),inp.end(),'i');
		switch (std::distance(ipos,inp.end())) {
			case 0: // number has only real part
			if (std::distance(i,e)>1) throw(std::invalid_argument("parsing error, bad format (1)"));
			return std::complex<RE__>(std::stod(*i));

			case 1: // number has imaginary part
				
			switch (std::distance(i,e)) {
				case 1: return std::complex<RE__>(0.0,1.0);
				case 2: return std::complex<RE__>(0.0,(*i=="+") ? 1.0: (*i=="-") ? -1.0: std::stod(*i));
				case 3:
				{
					RE__ re = RE__(std::stod(*i)); ++i;
					RE__ im = RE__((*i=="+") ? 1.0: (*i=="-") ? -1.0: std::stod(*i));
					return std::complex<RE__>(re,im);
				}
					
				default:
				throw(std::invalid_argument("parsing error, bad format (2)"));
			}
				
			default: // 'i' in bad position
			throw(std::invalid_argument("parsing error, bad format (3)"));
		}
	};

	// read data line by line
	readData_(file,pnfw);
}


// instantiation
template class lm_tMat<CPX__,RE__,CPX__>;
