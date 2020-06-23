// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "lm_defs.h"
#include "lm_types.h"
#include "lm_ops.h"
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
fMat::lm_tMat(const fArray& inp) noexcept: tMat(inp.M(),inp.N()) {
	if (inp.incr()==1)
		memcpy(data(),inp.data(),this->L()*sizeof(RE__));
	else {
		auto j=inp.begin();
		std::for_each(this->begin(),this->end(),[&j](RE__& i){ops::assign(i,*j++);});
	}
}


// assignment
template<>
fMat& fMat::operator=(const fArray& rhs) noexcept {
	if ((void*)this!=(void*)&rhs) {
		resize(rhs.M(),rhs.N());

		if (rhs.incr()==1)
			memcpy(data(),rhs.data(),this->L()*sizeof(RE__));
		else {
			auto j=rhs.begin();
			std::for_each(this->begin(),this->end(),[&j](RE__& i){ops::assign(i,*j++);});
		}
	}
	return *this;
}

// conversion
template<>
fMat fMat::fcopy() const noexcept {
	return copy();
}
template<>
cMat fMat::ccopy() const noexcept {
	return cMat(*this);
}


// matrix arithmetic
template<>
fCol& fMat::leftDivideEq(fCol& inp) const noexcept {
	// left matrix division (aka A\B) equals version, equivalent to B=A^-1*B
	assert(square());
	assert(N()==inp.M());
	if (M()==1) return inp/=*data();
		
	fMat tmp(*this);
	int* ipiv = new int [M()]; int info;
	c_xgesv(M(),1,tmp.data(),M(),ipiv,inp.data(),M(),&info);
	assert(!info);
	delete[] ipiv;

	return inp;
}
template<>
cCol& fMat::leftDivideEq(cCol& inp) const noexcept {
	// left matrix division (aka A\B) equals version, equivalent to B=A^-1*B
	assert(square());
	assert(N()==inp.M());
	if (M()==1) return inp/=*data();
	
	fMat LU(*this);
	int* ipiv = new int [M()]; int info;
	c_xgetrf(M(),N(),LU.data(),M(),ipiv,&info);
	assert(!info);

	// lambdas to solve system
	auto Ly_Pb = [&LU,ipiv](const cItr& bi, const cItr& be) mutable -> void {
		// solve Ly=Pb
		auto j = ipiv;
		for (auto i=bi; i!=be; ++i,++j)
			std::swap(*i,*(bi+*j-1));
		
		auto c=LU.ccBegin(); auto i=bi;
		for (size_t m=1; m!=LU.M(); ++m,++c,++i) {
			auto cj=c->begin()+m;
			for (auto ij=i+1; ij!=be; ++ij,++cj)
				*ij -= (*i)*(*cj);
		}
		
	};
	auto Ux_y = [&LU](const r_cItr& bi, const r_cItr& be) -> void {
		// solve Ux=y
		auto c=LU.crcBegin(); auto i=bi;
		for (size_t m=0; m!=LU.M(); ++m,++c,++i) {
			auto cj=c->rbegin()+m;
			*i/=*cj++;
			for (auto ij=i+1; ij!=be; ++ij,++cj)
				*ij -= (*i)*(*cj);
		}
	};

	// solve for inp
	Ly_Pb(inp.begin(),inp.end());
	delete[] ipiv;
	Ux_y(inp.rbegin(),inp.rend());


	return inp;
}
template<>
fMat& fMat::leftDivideEq(fMat& inp) const noexcept {
	// left matrix division (aka A\B) equals version, equivalent to B=A^-1*B
	assert(square());
	assert(N()==inp.M());
	if (M()==1) return inp/=*data();
		
	fMat tmp(*this);
	int* ipiv = new int [M()]; int info;
	c_xgesv(M(),inp.N(),tmp.data(),M(),ipiv,inp.data(),M(),&info);
	assert(!info);
	delete[] ipiv;

	return inp;
}
template<>
cMat& fMat::leftDivideEq(cMat& inp) const noexcept {
	// left matrix division (aka A\B) equals version, equivalent to B=A^-1*B
	assert(square());
	assert(N()==inp.M());
	if (M()==1) return inp/=*data();
	
	fMat LU(*this);
	int* ipiv = new int [M()]; int info;
	c_xgetrf(M(),N(),LU.data(),M(),ipiv,&info);
	assert(!info);

	// lambdas to solve system
	auto Ly_Pb = [&LU,ipiv](const cItr& bi, const cItr& be) mutable -> void {
		// solve Ly=Pb
		auto j = ipiv;
		for (auto i=bi; i!=be; ++i,++j)
			std::swap(*i,*(bi+*j-1));
		
		auto c=LU.ccBegin(); auto i=bi;
		for (size_t m=1; m!=LU.M(); ++m,++c,++i) {
			auto cj=c->begin()+m;
			for (auto ij=i+1; ij!=be; ++ij,++cj)
				*ij -= (*i)*(*cj);
		}
		
	};
	auto Ux_y = [&LU](const r_cItr& bi, const r_cItr& be) -> void {
		// solve Ux=y
		auto c=LU.crcBegin(); auto i=bi;
		for (size_t m=0; m!=LU.M(); ++m,++c,++i) {
			auto cj=c->rbegin()+m;
			*i/=*cj++;
			for (auto ij=i+1; ij!=be; ++ij,++cj)
				*ij -= (*i)*(*cj);
		}
	};

	// solve for inp
	for (auto i=inp.cBegin(),e=inp.cEnd(); i!=e; ++i) {
		Ly_Pb(i->begin(),i->end());
	}
	delete[] ipiv;
	for (auto i=inp.cBegin(),e=inp.cEnd(); i!=e; ++i)
		Ux_y(i->rbegin(),i->rend());

	return inp;
}
template<>
fMat fMat::prod(const fArray& inp) const noexcept {
	assert(!empty());
	assert(!inp.empty());
	assert(N()==inp.M());
	
	fMat res(M(),inp.N());
	c_xgemm('N','N',M(),inp.N(),N(),1.0,data(),M(),inp.data(),N(),0.0,res.data(),M());
	
	return res;
}
template<>
cMat fMat::prod(const cArray& inp) const noexcept {
	assert(!empty());
	assert(!inp.empty());
	assert(N()==inp.M());
	
	cMat res(M(),inp.N());

	const RE__* inp_ptr = c_reItr(inp.begin()).data();
	const RE__* const inp_end = c_reItr(inp.end()).data();
	RE__* res_ptr = reItr(res.begin()).data();

	const size_t resM=2*res.M(), inpM=2*inp.M();
	while (inp_ptr<inp_end) {
		c_xgemv('N',M(),N(),1.0,data(),M(),inp_ptr,2,0.0,res_ptr,2);
		c_xgemv('N',M(),N(),1.0,data(),M(),inp_ptr+1,2,0.0,res_ptr+1,2);
		inp_ptr+=inpM; res_ptr+=resM;
	}

	return res;
}


// member functions
template<>
void fMat::readBinary_(std::ifstream& file, const bool cpx) {
	uint64_t M,N;
	file.read((char*) &M, sizeof(uint64_t));
	file.read((char*) &N, sizeof(uint64_t));
	(*this) = fMat(M,N);
	
	if (cpx) {
		RE__ buff;
		for (auto i=data(),e=data()+L(); i!=e; ++i) {
			file.read((char*) i, sizeof(RE__));
			file.read((char*) &buff, sizeof(RE__));
		}
	} else file.read((char*) data(),L()*sizeof(RE__));
}
template<>
void fMat::parse_(std::ifstream& file) {
	assert(file.good());
	assert(!empty());
	
	// regex
	const std::regex rn("(?=[-+i])");

	// lambda to parse numbers from string
	auto pnfw = [&rn](const std::string& inp) -> RE__ {
		// tokenize number
		std::sregex_token_iterator i(inp.begin(),inp.end(),rn,-1), e;
		if (!i->length()) ++i; // inp starts with delimiter

		// find position of 'i'
		auto ipos = std::find(inp.begin(),inp.end(),'i');
		switch (std::distance(ipos,inp.end())) {
			case 0: // number has only real part
			if (std::distance(i,e)>1) throw(std::invalid_argument("parsing error, bad format (1)"));
			return RE__(std::stod(*i));

			case 1: // number has imaginary part
			
			switch (std::distance(i,e)) {
				case 1: return RE__(0.0);		// number is 'i'
				case 2: return RE__(0.0);		// number is only imaginary
				case 3: return RE__(std::stod(*i));	// number has real part
			
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
template class lm_tMat<RE__,RE__,CPX__>;
