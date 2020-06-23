// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "lm_tMat.h"
#include "lm_fn.h"
#include "lm_ops.h"
#include "blas.h"
#include <cassert>
#include <regex>
#include <fstream>


using namespace lm__;

// constructor
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>::lm_tMat(const fArray& re, const fArray& im) noexcept: tMat(re) {
	assert(re.M()==im.M());
	assert(re.N()==im.N());
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>::lm_tMat(const fArray& inp) noexcept: tMat(inp.M(),inp.N()) {
	auto j=inp.begin();
	std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::assign(i,*j++);});
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>::lm_tMat(const cArray& inp) noexcept: tMat(inp.M(),inp.N()) {
	auto j=inp.begin();
	std::for_each(this->begin(),this->end(),[&j](TT& i){lm__::ops::assign(i,*j++);});
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>::lm_tMat(const std::string& fileName): tMat(0,0) {
	// open file
	std::ifstream file;
	file.open(fileName,std::ios::binary);
	if (!file.good()) throw(std::invalid_argument("file \'"+fileName+"\' not found"));

	// read header
	char buff[5];
	file.read(buff,5);

	// read binary or text depending on header
	if (strcmp(buff,"fMat") && strcmp(buff,"cMat")) {
		file.close();
		readText_(fileName);
	} else {
		uint64_t dat;
		file.read((char*) &dat,sizeof(uint64_t));
		if (dat!=2)
			throw(std::invalid_argument("data in file '"+fileName+"' is not 2 dimensional"));
		readBinary_(file,strcmp(buff,"fMat"));
		file.close();
	}
}


// memory management
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::reserve(const size_t nC) noexcept {
	if (M() && nC>ccap()) {
		TT* ndata = new TT [M()*nC];
		memcpy(ndata,data_,this->L()*sizeof(TT));
		delete[] data_;
		data_ = ndata;
		C_ = M()*nC;
	}
	return *this;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::shrink_to_fit() noexcept {
	if (this->L()<lcap()) {
		TT* ndata = new TT [this->L()];
		memcpy(ndata,data_,this->L()*sizeof(TT));
		delete[] data_;

		data_ = ndata;
		C_ = this->L();
	}
	return *this;
}
template<class TT, class FT, class CT>
TT* lm_tMat<TT,FT,CT>::move() noexcept {
	auto tmp = data_;
	data_ = nullptr;
	N_ = C_ = 0;
	return tmp;
}


// basic properties
template<class TT, class FT, class CT>
bool lm_tMat<TT,FT,CT>::hermitian() const noexcept {
	if (!square()) return false;

	for (size_t n=0; n!=N(); ++n)
		for (size_t m=0; m!=M(); ++m)
			if (lm__::ops::neq((*this)(m,n),std::conj((*this)(n,m))))
				return false;
	return true;
}
template<class TT, class FT, class CT>
bool lm_tMat<TT,FT,CT>::diag() const noexcept {
	auto i=this->begin(), e = M()<N() ? this->begin()+M()*M(): this->end();
		
	while (i<e) {
		if (lm__::ops::z(*i++)) return false;
		
		auto ib = i+M();
		while (i<ib && i<e)
			if (lm__::ops::nz(*i++)) return false;
	}
	while (i<this->end()) 
		if (lm__::ops::nz(*i++)) return false;
	
	return this->L();
}
template<class TT, class FT, class CT>
bool lm_tMat<TT,FT,CT>::ob() const noexcept {
	if (!square()) return false;

	const auto n = mnorm(*this);
	if (!all(n)) return false;

	auto ci=ccBegin();
	for (size_t i=0; i!=M(); ++i, ++ci) {
		auto cj = ci+1;
		for (size_t j=i+1; j!=M(); ++j, ++cj)
			if (ops::nz(dot(*ci,*cj)/(n[i]*n[j])))
				return false;
	}
	return true;
}
template<class TT, class FT, class CT>
bool lm_tMat<TT,FT,CT>::onb() const noexcept {
	if (!square()) return false;
	for (auto i=ccBegin(),e=ccEnd(); i!=e; ++i) {
		auto j=i;
		if (ops::neq(dot(*i,*j++),TT(1.0)))
			return false;
		while (j!=e)
			if (ops::nz(dot(*i,*j++)))
				return false;
	}
	return true;
}


// basic modification
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::T() noexcept {
	if (this->empty() || M()==1 || N()==1) {
		std::swap(M_,N_);
		return *this;
	}	
	if (square()) {
		auto ri = rBegin();
		auto ci = cBegin();
		for (size_t i=1,l=this->M(); i!=l; ++i,++ri,++ci)
			for (auto rj=ri->begin()+i,cj=ci->begin()+i,ce=ci->end(); cj!=ce; ++rj,++cj)
				std::swap(*rj,*cj);
		return *this;
	}

	tMat tmp(N(),0); tmp.reserve(M());
	for (auto i=crBegin(),e=crEnd(); i!=e; ++i)
		tmp.push_back(*i);
	return *this=std::move(tmp);
}


// modification
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::rShift(const size_t m, const size_t s) noexcept {
	assert(m<M());
	if (!s) return *this;

	// insufficient space
	if (lcap()<(M()+s)*N()) {
		tMat tmp(M()+s,N());
		auto i=data(), it=tmp.data();
		for (size_t n=0; n!=N(); ++n, i+=M(), it+=tmp.M()) {
			memcpy(it,i,m*sizeof(TT));
			memcpy(it+m+s,i+m,(M()-m)*sizeof(TT));
		}
		return *this=std::move(tmp);
	}

	// sufficient space
	// move tail
	const auto src_ = data()+M()*(N()-1)+m;
	const auto dst_ = src_+N()*s;
	if (N()*s>M()-m) memcpy(dst_,src_,(M()-m)*sizeof(TT));
	else memmove(dst_,src_,(M()-m)*sizeof(TT));
	
	// move rest
	if (N()>1) {
		auto src_ = data()+M()*(N()-2)+m;
		auto dst_ = src_+(N()-1)*s;

		for (; src_>=data() && size_t(dst_-src_)>M(); src_-=M(), dst_-=(M()+s))
			memcpy(dst_,src_,M()*sizeof(TT));
		for (; src_>=data(); src_-=M(), dst_-=(M()+s))
			memmove(dst_,src_,M()*sizeof(TT));
	}
	
	M_+=s;
	return *this;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::cShift(const size_t n, const size_t s) noexcept {
	assert(n<=N());
	
	// insufficent space
	if (ccap()<N()+s) {
		tMat tmp(M(),N()+s);
		if (n) memcpy(tmp.data(),data(),M()*n*sizeof(TT));
		memcpy(tmp.data()+M()*(n+s),data()+M()*n,M()*(N()-n)*sizeof(TT));
		return *this=std::move(tmp);
	}

	// sufficient space
	if (s<N()-n) memmove(data()+M()*(n+s),data()+M()*n,M()*(N()-n)*sizeof(TT));
	else memcpy(data()+M()*(n+s),data()+M()*n,M()*(N()-n)*sizeof(TT));

	N_+=s;
	return *this;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::rInsert(const size_t m, const tArray& inp) noexcept {
	if (inp.ptr()!=this) {
		rShift(m,1);
		rAt(m)=inp;
		return *this;
	} else return rInsert(m,tMat(inp));
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::rInsert(const size_t m, const tMat& inp) noexcept {
	if (inp.ptr()!=this) {
		rShift(m,inp.M());
		return set(inp,m,0);
	} else return rInsert(m,tMat(inp));
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::cInsert(const size_t n, const tArray& inp) noexcept {
	if (inp.ptr()!=this) {
		cShift(n,1);
		cAt(n)=inp;
		return *this;
	} else return cInsert(n,tMat(inp));
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::cInsert(const size_t n, const tMat& inp) noexcept {
	if (inp.ptr()!=this) {
		cShift(n,inp.N());
		return set(inp,0,n);
	} else return cInsert(n,tMat(inp));
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::set(const tMat& inp,
		const std::vector<size_t>& m, const std::vector<size_t>& n) noexcept {
	if (!m.size() && !n.size()) {
		assert(msize(*this)==msize(inp));
		return *this = inp;
	}
	if (!m.size()) {
		assert(inp.M()==M() && inp.N()==n.size());
		assert(std::all_of(n.begin(), n.end(), [&](size_t i){return i<N();}));
		
		auto j = inp.ccBegin();
		for (size_t i=0; i!=inp.N(); ++i,++j)
			cAt(n[i]) = *j;
		return *this;
	}
	if (!n.size()) {
		assert(inp.M()==m.size() && inp.N()==N());
		assert(std::all_of(m.begin(), m.end(), [&](size_t i){return i<M();}));
		
		auto j = inp.crBegin();
		for (size_t i=0; i!=inp.M(); ++i,++j)
			rAt(m[i]) = *j;
		return *this;
	}

	assert(inp.M()==m.size());
	assert(inp.N()==n.size());
	assert(std::all_of(m.begin(), m.end(), [&](size_t i){return i<M();}));
	assert(std::all_of(n.begin(), n.end(), [&](size_t i){return i<N();}));

	auto ii = inp.cbegin();
	for (auto j: n)
		for (auto i: m)
			(*this)(i,j) = *ii++;
	return *this;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::set(const tMat& inp, const size_t m, const size_t n) noexcept {
	assert(m+inp.M()<=M());
	assert(n+inp.N()<=N());

	auto s=inp.data();
	for (auto d=data()+n*M()+m,e=data()+(n+inp.N())*M()+m; d!=e; d+=M(),s+=inp.M())
		memcpy(d,s,inp.M()*sizeof(TT));
	return *this;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::setl(const tMat& inp, const fArray& m, const fArray& n) noexcept {
	assert(m.L()==M());
	assert(n.L()==N());
	assert(inp.M()==size_t(sum(~~m)));
	assert(inp.N()==size_t(sum(~~n)));

	auto ii = inp.cbegin();
	for (size_t j=0; j!=N(); ++j)
		if (n[j]) {
			for (size_t i=0; i!=M(); ++i)
				if (m[i]) (*this)(i,j) = *ii++;
		}
	return *this;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::rRm(const size_t m) noexcept {
	// remove row by moving the data in place, no dealloc
	assert(!this->empty());
	assert(m<M());
	
	auto i=data()+m; size_t move=0; size_t max_move=N()-1;
	while(move!=max_move) {
		memmove(i-move, i+1, (M()-1)*sizeof(TT));
		i+=M(); ++move;
	}
	memmove(i-move, i+1, (M()-m-1)*sizeof(TT));
	--M_;
	return *this;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::rRm(const std::vector<size_t>& m) noexcept {
	// remove rows by moving the data in place, no dealloc or new alloc
	if (m.empty()) return *this;
	assert(std::is_sorted(m.begin(),m.end()));
	assert(std::all_of(m.begin(),m.end()-1,[](const size_t& i){return i!=*(&i+1);}));
	assert(std::all_of(m.begin(),m.end(),[&](const size_t& i){return i<M();}));
	if (m.size()==M()) return resize(0,N());
	
	std::vector<size_t> s; s.reserve(m.size());
	std::vector<size_t> k; k.reserve(m.size());
	
	size_t cnt=0;
	for (size_t i=0; i!=m.size()-1; ++i)
		if (m[i+1]-m[i]==1) ++cnt;
		else {
			k.push_back(cnt+1);
			s.push_back(m[i+1]-m[i]+cnt);
			cnt=0;
		}
	k.push_back(cnt+1);
	s.push_back(M()-m.back()+m.front()+cnt);
	
	auto i=data()+m[0]; size_t move=0; size_t max_move=m.size()*N()-k.back(); size_t j=0;
	while (move!=max_move) {
		memmove(i-move,i+k[j],(s[j]-k[j])*sizeof(TT));
		move+=k[j]; i+=s[j]; ++j%=k.size();
	}
	memmove(i-move, i+k[j], (M()-m.back()-1)*sizeof(TT));
	M_-=m.size();

	return *this;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::cRm(const size_t n) noexcept {
	// remove col by moving the data in place, no dealloc
	assert(!this->empty());
	assert(n<N());
	
	memmove(data()+n*M(), data()+(n+1)*M(), (N()-n-1)*M()*sizeof(TT));
	--N_;
	return *this;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::cRm(const std::vector<size_t>& n) noexcept {
	// remove cols by moving the data in place, no dealloc or new alloc
	if (n.empty()) return *this;
	assert(std::is_sorted(n.begin(),n.end()));
	assert(std::all_of(n.begin(),n.end()-1,[](const size_t& i){return i!=*(&i+1);}));
	assert(std::all_of(n.begin(),n.end(),[&](const size_t& i){return i<N();}));
	if (n.size()==N()) return resize(0);

	size_t i;
	for (i=0; i!=n.size()-1; ++i)
		memmove(data()+(n[i]-i)*M(), data()+(n[i]+1)*M(), (n[i+1]-n[i]-1)*M()*sizeof(TT));
	memmove(data()+(n[i]-i)*M(), data()+(n[i]+1)*M(), (N()-n[i]-1)*M()*sizeof(TT));
	N_-=n.size();
	
	return *this;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::dRm(const bool low) noexcept {
	// remove diagonal by moving the data in place, no dealloc or new alloc
	// if this is diagonal, parts may be added from the left or from below
	// if this has more rows than cols, parts must be added from the left
	// if this has more cols than rows, parts must be added from below
	assert(!this->empty());

	if (this->L()==1) {
		M_=N_=0;
		return *this;
	}

	auto i=data()+M(); size_t l;
	if (N()>M() || (square()&&low)) {
		// add parts from below
		for (l=1; l!=M(); i+=M(), ++l)
			memcpy(i-M(),i,l*sizeof(TT));
		memmove(i-M(),i,(data()+this->L()-i)*sizeof(TT));
		--N_;
	} else {
		// add parts from the left
		for (l=M()-1; l!=M()-N(); i+=M(), --l)
			memmove(i-M(),i-l,M()*sizeof(TT));
		memmove(i-M(),i-l,l*sizeof(TT));
		--M_;
	}
	return *this;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm_tMat<TT,FT,CT>::inv() noexcept {
	assert(square());

	int* ipiv = new int [N()]; int info;
	c_xgetrf(M(),M(),data(),M(),ipiv,&info);
	assert(!info);

	TT* work = new TT [M()];
	c_xgetri(M(),data(),M(),ipiv,work,M(),&info);
	assert(!info);

	delete[] ipiv; delete[] work;
	return *this;
}


// conversion
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tMat<TT,FT,CT>::get(const std::vector<size_t>& m, const std::vector<size_t>& n) const noexcept {
	if (!m.size() && !n.size()) return *this;
	if (!m.size()) {
		assert(std::all_of(n.begin(), n.end(), [&](size_t i){return i<N();}));
		
		tMat res(M(),n.size());
		for (size_t i=0; i!=res.N(); ++i)
			memcpy(res.data()+i*M(),data()+n[i]*M(),M()*sizeof(TT));
		return res;
	}
	if (!n.size()) {
		assert(std::all_of(m.begin(), m.end(), [&](size_t i){return i<M();}));
		
		tMat res(m.size(),N());
		for (size_t i=0; i!=res.N(); ++i)
			for (size_t j=0; j!=res.M(); ++j)
				res(j,i) = (*this)(m[j],i);
		return res;
	}


	assert(std::all_of(m.begin(), m.end(), [&](size_t i){return i<M();}));
	assert(std::all_of(n.begin(), n.end(), [&](size_t i){return i<N();}));
	
	tMat res(m.size(),n.size());
	auto r = res.begin();
	for (auto i: n)
		for (auto j: m)
			*r++ = (*this)(j,i);
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tMat<TT,FT,CT>::get(const size_t m, const size_t n,
		const size_t lm, const size_t ln) const noexcept {
	assert(m<M());
	assert(n<N());
	
	tMat res(M()-m<lm ? M()-m: lm, N()-n<ln ? N()-n: ln);
	auto dptr = res.data(), eptr = res.data()+res.L(); auto sptr = data()+n*M()+m;
	while (dptr!=eptr) {
		memcpy(dptr, sptr, res.M()*sizeof(TT));
		dptr+=res.M(); sptr+=M();
	}
	
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tMat<TT,FT,CT>::getl(const fArray& m, const fArray& n) const noexcept {
	assert(m.col());
	assert(m.M()==M());
	assert(n.row());
	assert(n.N()==N());
		
	tMat res(M()-size_t(sum(~m)),N()-size_t(sum(~n)));
	
	size_t ri=0; size_t rj=0;
	for (size_t j=0; j!=N(); ++j)
		if (bool(n[j])) {
			for (size_t i=0; i!=M(); ++i)
				if (bool(m[i]))
					res(ri++,rj) = (*this)(i,j);
			ri=0; ++rj;
		}
		
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tMat<TT,FT,CT>::getl(const fMat& inp) const noexcept {
	assert(msize(*this)==msize(inp));

	size_t cnt=0;
	for (auto i: inp)
		if (i) ++cnt;
	tMat res(cnt,1);

	auto i=this->cbegin(); auto r=res.begin();
	for (auto j=inp.cbegin(),e=inp.cend(); j!=e; ++j,++i)
		if (*j) *r++ = *i;

	return res;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tMat<TT,FT,CT>::rGet(const size_t m, const size_t l) const noexcept {
	assert(!this->empty());
	assert(m+l<=M());
	
	return !l ? tMat(0,N()): get(m,0,l,N());
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tMat<TT,FT,CT>::cGet(const size_t n, const size_t l) const noexcept {
	assert(!this->empty());
	assert(n+l<=N());

	tMat res(M(),l);
	memcpy(res.data(),data()+n*M(),M()*l*sizeof(TT));
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tMat<TT,FT,CT>::rWOGet(const size_t m) const noexcept {
	assert(m<M());
	
	tMat res(M()-1,N());
	auto res_ptr=res.data(), res_end=res.data()+res.L();
	auto this_ptr=data();
	const size_t rm=M()-1-m;
	
	while (res_ptr!=res_end) {
		memcpy(res_ptr, this_ptr, m*sizeof(TT));
		res_ptr+=m; this_ptr+=m+1;
		memcpy(res_ptr, this_ptr, rm*sizeof(TT));
		res_ptr+=rm; this_ptr+=rm;
	}
	
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tMat<TT,FT,CT>::cWOGet(const size_t n) const noexcept {
	assert(n<N());
	
	tMat res(M(),N()-1);
	auto res_ptr=res.data();
	auto this_ptr=data();
	
	memcpy(res_ptr, this_ptr, n*M()*sizeof(TT));
	res_ptr+=n*M(); this_ptr+=(n+1)*M();
	memcpy(res_ptr, this_ptr, (N()-n-1)*M()*sizeof(TT));
	
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tMat<TT,FT,CT>::upper() const noexcept {
	if (this->empty()) return tMat(0,1);
	
	const size_t n = M()<N() ? M(): N();
	const size_t d = M()<N() ? N()-M(): 0;
	tMat res(n*(n-1)/2 + d*M(),1);

	size_t l=1;
	auto ditr=res.data(); auto sitr=data()+M();
	while (l!=n) {
		memcpy(ditr,sitr,l*sizeof(TT));
		ditr+=l; sitr+=M(); ++l;
	}
	if (d) memcpy(ditr,sitr,d*M()*sizeof(TT));

	return res;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm_tMat<TT,FT,CT>::lower() const noexcept {
	if (this->empty()) return tMat(0,1);

	const size_t n = M()<N() ? M(): N();
	const size_t d = M()<N() ? 0: M()-N();
	tMat res(n*(n-1)/2 + d*N(),1);

	size_t l=M()-1;
	auto ditr=res.data(); auto sitr=data()+1, e=data()+this->L();
	while (l && sitr<e) {
		memcpy(ditr,sitr,l*sizeof(TT));
		ditr+=l; sitr+=(M()+1); --l;
	}

	return res;
}


// logical
template<class TT, class FT, class CT>
bool lm_tMat<TT,FT,CT>::permutation() const noexcept {
	if (!square()) return false;

	// bools to check each 1 is in unique row
	// one extra space to check if there is a 1
	std::vector<bool> nf(M()+1,true);
	
	auto i=this->begin(), e=this->end();
	while (i!=e) {
		nf.back()=true;
		for (size_t m=0; m<M(); ++i,++m) {
			if (nf[m] && lm__::ops::eq_s(*i,FT(1.0))) {
				nf[m]=nf.back()=false;
				continue;
			}
			if (lm__::ops::nz_s(*i))
				return false;
		}
		if (nf.back()) return false;
	}
	return std::none_of(nf.begin(),nf.end(),[](const bool i){return i;});
}


// printing
template<class TT, class FT, class CT>
void lm_tMat<TT,FT,CT>::writeToFile(const std::string& fileName, const bool noheader) const {
	// open file
	std::ofstream file;
	file.open(fileName,std::ios::binary);
	if (!file.good())
		throw(std::invalid_argument("open file \'"+fileName+"\' failed"));

	// write header
	if (!noheader) {
		const std::string hdr = this->cpx() ? "cMat": "fMat";
		file.write(hdr.c_str(),5);
	}

	// write dimensions
	uint64_t dat = 2;
	file.write((char*) &dat, sizeof(uint64_t));
	dat = this->M();
	file.write((char*) &dat, sizeof(uint64_t));
	dat = this->N();
	file.write((char*) &dat, sizeof(uint64_t));

	// write matrix data
	file.write((char*) data(), this->L()*sizeof(TT));
	
	file.flush();
	file.close();
}


// member functions
template<class TT, class FT, class CT>
void lm_tMat<TT,FT,CT>::readText_(const std::string& fileName) {
	// open file
	std::ifstream file;
	file.open(fileName);
	if (!file.good()) throw(std::invalid_argument("file \'"+fileName+"\' not found"));

	// regex
	const std::regex r("[\\s,;]+");

	// line and word counters
	size_t Lc=0, Wc=0;

	// run through file to find line and word count
	std::string line;
	while (std::getline(file,line)) {
		// tokenize
		std::sregex_token_iterator i(line.begin(),line.end(),r,-1), e;
		if (!i->length()) ++i; // line starts with delimiter

		// check word count
		if (!Lc) Wc = std::distance(i,e);
		else {
			if (size_t(std::distance(i,e))!=Wc)
				throw(std::invalid_argument("file \'"+fileName+"\', line "+
					std::to_string(Lc+1)+": found "+std::to_string(std::distance(i,e))+
					" entries, need "+std::to_string(Wc)+" entries"));
		}
		
		++Lc;
	}
	file.close();
	if (!Lc) return;

	// all good, allocate memory and read data
	file.open(fileName);
	try {
		*this = tMat(Lc,Wc);
		parse_(file);
	} catch (const std::exception& err) {
		file.close();
		throw(std::invalid_argument("file \'"+fileName+"\', "+err.what()));
	}
	file.close();
}
template<class TT, class FT, class CT>
void lm_tMat<TT,FT,CT>::readData_(std::ifstream& file, const std::function<TT(std::string)> pnfw) {
	assert(file.good());
	assert(!this->empty());

	// regex
	const std::regex rl("[\\s,;]+");
	
	// read data line by line
	std::string line;
	for (size_t m=0; m!=M(); ++m) {
		// read line, check for eof
		std::getline(file,line);
		if (!file) throw(std::invalid_argument("unexpected eof"));
		
		// tokenize line
		std::sregex_token_iterator i(line.begin(),line.end(),rl,-1), e;
		if (!i->length()) ++i; // line starts with delimiter
		if (size_t(std::distance(i,e))<N())
			throw(std::invalid_argument("number of words on line "+std::to_string(m+1)+" is "+
				std::to_string(std::distance(i,e))+", expected >= "+std::to_string(N())));
		
		// read numbers from line
		for (size_t n=0; n!=N(); ++n) {
			const std::string word = *i++;
			
			// check word for alien symbols
			size_t bpos = word.find_first_not_of("-+0123456789.ei ,;");
			if (bpos<std::string::npos)
				throw(std::invalid_argument("parsing error, found alien symbol \'"+
				word.substr(bpos,1)+"\'"));
			
			// parse number from word
			(*this)(m,n)=pnfw(word);
		}
	}
}
