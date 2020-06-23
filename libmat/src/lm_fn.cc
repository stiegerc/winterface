// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "lm_fn.h"
#include "lm_types.h"
#include "lm_cpxItr.h"
#include "lm_ops.h"
#include "blas.h"
#include <algorithm>
#include <random>
#include <cstring>
#include <cfloat>

using namespace lm__;

// tolerance default
#ifndef NTOLERANT__
#ifndef NVARTOL__
// variable tolerance implemented through stack(s)
#ifdef _OPENMP
std::vector<std::stack<RE__>> lm__::hide::mtol_st_(1,
	std::stack<RE__>(std::deque<RE__>{MTOL__}));
#else
std::stack<RE__> lm__::hide::mtol_st_(std::deque<RE__>{MTOL__});
#endif
#else
// no variable tolerance
constexpr RE__ lm__::hide::mtol_ = MTOL__;
#endif
#endif



// math functions
template<class TT, class FT, class CT>
fMat lm__::abs(const lm_tArray<TT,FT,CT>& inp) noexcept {
	fMat res(msize(inp));
	
	auto j=inp.cbegin();
	for (auto& i: res)
		i = std::abs(*j++);
	return res;
}
template fMat lm__::abs(const fArray& inp) noexcept;
template fMat lm__::abs(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
lm_tArray<TT,FT,CT>& lm__::ceilEq(lm_tArray<TT,FT,CT>& inp) noexcept {
	std::for_each(inp.begin(),inp.end(),[](TT& i){i=ops::ceil(i);});
	return inp;
}
template fArray& lm__::ceilEq(fArray& inp) noexcept;
template cArray& lm__::ceilEq(cArray& inp) noexcept;

template<class TT, class FT, class CT>
lm_tArray<TT,FT,CT>& lm__::floorEq(lm_tArray<TT,FT,CT>& inp) noexcept {
	std::for_each(inp.begin(),inp.end(),[](TT& i){i=ops::floor(i);});
	return inp;
}
template fArray& lm__::floorEq(fArray& inp) noexcept;
template cArray& lm__::floorEq(cArray& inp) noexcept;

template<class TT, class FT, class CT>
lm_tArray<TT,FT,CT>& lm__::roundEq(lm_tArray<TT,FT,CT>& inp) noexcept {
	std::for_each(inp.begin(),inp.end(),[](TT& i){i=ops::round(i);});
	return inp;
}
template fArray& lm__::roundEq(fArray& inp) noexcept;
template cArray& lm__::roundEq(cArray& inp) noexcept;

template<class TT, class FT, class CT>
lm_tArray<TT,FT,CT>& lm__::signEq(lm_tArray<TT,FT,CT>& inp) noexcept {
	std::for_each(inp.begin(),inp.end(),[](TT& i){i=ops::sign(i);});
	return inp;
}
template fArray& lm__::signEq(fArray& inp) noexcept;
template cArray& lm__::signEq(cArray& inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::ceil(const lm_tArray<TT,FT,CT>& inp) noexcept {
	lm_tMat<TT,FT,CT> res(inp);
	return ceilEq(res);
}
template fMat lm__::ceil(const fArray& inp) noexcept;
template cMat lm__::ceil(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::floor(const lm_tArray<TT,FT,CT>& inp) noexcept {
	lm_tMat<TT,FT,CT> res(inp);
	return floorEq(res);
}
template fMat lm__::floor(const fArray& inp) noexcept;
template cMat lm__::floor(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::round(const lm_tArray<TT,FT,CT>& inp) noexcept {
	lm_tMat<TT,FT,CT> res(inp);
	return roundEq(res);
}
template fMat lm__::round(const fArray& inp) noexcept;
template cMat lm__::round(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::sign(const lm_tArray<TT,FT,CT>& inp) noexcept {
	lm_tMat<TT,FT,CT> res(inp);
	return signEq(res);
}
template fMat lm__::sign(const fArray& inp) noexcept;
template cMat lm__::sign(const cArray& inp) noexcept;


// eigenvalue computation
cMat lm__::eig(fMat inp) noexcept {
	assert(inp.square());
	assert(inp.M()>=2);

	int info;
	fMat work(4*inp.N(),1);

	fMat wr(inp.M(),1), wi(inp.M(),1);
	c_xgeev('N','N',inp.M(),inp.data(),inp.M(),wr.data(),wi.data(),
		nullptr,1,nullptr,1,work.data(),work.M(),&info);
	cMat w(wr,wi);

	assert(!info);
	return w;
}
cMat lm__::eig(cMat inp) noexcept {
	assert(inp.square());
	assert(inp.M()>=2);

	int info;
	cMat work(4*inp.N(),1);

	cMat w(inp.M(),1);
	fMat rwork(2*inp.M(),1);
	c_xgeev('N','N',inp.M(),inp.data(),inp.M(),w.data(),nullptr,1,nullptr,
		1,work.data(),work.M(),rwork.data(),&info);

	assert(!info);
	return w;
}
fevr lm__::eigr(fMat inp) noexcept {
	assert(inp.square());
	assert(inp.M()>=2);

	// lapack
	int info;
	fMat work(4*inp.N(),1);
	fMat wr(inp.M(),1), wi(inp.M(),1);
	fMat Vr(inp.M(),inp.N());
	c_xgeev('N','V',inp.M(),inp.data(),inp.M(),wr.data(),wi.data(),
		nullptr,1,Vr.data(),Vr.M(),work.data(),work.M(),&info);
	assert(!info);

	return {std::move(wr),std::move(wi),std::move(Vr)};
}
cevr lm__::eigr(cMat inp) noexcept {
	assert(inp.square());
	assert(inp.M()>=2);

	int info;
	cMat work(4*inp.N(),1);
	cMat Vr(inp.M(),inp.N());

	cMat w(inp.M(),1);
	fMat rwork(2*inp.M(),1);
	c_xgeev('N','V',inp.M(),inp.data(),inp.M(),w.data(),nullptr,1,
		Vr.data(),Vr.M(),work.data(),work.M(),rwork.data(),&info);

	assert(!info);
	return {std::move(w),std::move(Vr)};
}
fevl lm__::eigl(fMat inp) noexcept {
	assert(inp.square());
	assert(inp.M()>=2);

	// lapack
	int info;
	fMat work(4*inp.N(),1);
	fMat wr(inp.M(),1), wi(inp.M(),1);
	fMat Vl(inp.M(),inp.N());
	c_xgeev('V','N',inp.M(),inp.data(),inp.M(),wr.data(),wi.data(),
		Vl.data(),Vl.M(),nullptr,1,work.data(),work.M(),&info);
	assert(!info);

	return {std::move(wr),std::move(wi),std::move(Vl)};
}
cevl lm__::eigl(cMat inp) noexcept {
	assert(inp.square());
	assert(inp.M()>=2);

	int info;
	cMat work(4*inp.N(),1);
	cMat Vl(inp.M(),inp.N());

	cMat w(inp.M(),1);
	fMat rwork(2*inp.M(),1);
	c_xgeev('V','N',inp.M(),inp.data(),inp.M(),w.data(),Vl.data(),Vl.M(),
		nullptr,1,work.data(),work.M(),rwork.data(),&info);
	assert(!info);

	return {std::move(w),std::move(Vl)};
}
fevrl lm__::eigrl(fMat inp) noexcept {
	assert(inp.square());
	assert(inp.M()>=2);

	int info;
	fMat work(4*inp.N(),1);
	fMat wr(inp.M(),1), wi(inp.M(),1);
	fMat Vr(inp.M(),inp.N());
	fMat Vl(inp.M(),inp.N());
	c_xgeev('V','V',inp.M(),inp.data(),inp.M(),wr.data(),wi.data(),
		Vl.data(),Vl.M(),Vr.data(),Vr.M(),work.data(),work.M(),&info);
	assert(!info);

	return {std::move(wr),std::move(wi),std::move(Vr),std::move(Vl)};
}
cevrl lm__::eigrl(cMat inp) noexcept {
	assert(inp.square());
	assert(inp.M()>=2);

	int info;
	cMat work(4*inp.N(),1);
	cMat Vl(inp.M(),inp.N());
	cMat Vr(inp.M(),inp.N());

	cMat w(inp.M(),1);
	fMat rwork(2*inp.M(),1);
	c_xgeev('V','V',inp.M(),inp.data(),inp.M(),w.data(),Vl.data(),Vl.M(),
			Vr.data(),Vr.M(),work.data(),work.M(),rwork.data(),&info);
	assert(!info);
	return {std::move(w),std::move(Vr),std::move(Vl)};
}
env lm__::eGet(const fMat& Er, const fMat& Ei, const fMat& V, const size_t i) noexcept {
	assert(i<V.N());

	if (i+1<V.N() && Er[i]==Er[i+1] && Ei[i]==-Ei[i+1])
		return {{Er[i],Ei[i]},cMat(V.cAt(i),V.cAt(i+1))};
	if (i>0 && Er[i]==Er[i-1] && Ei[i]==-Ei[i-1])
		return {{Er[i],Ei[i]},cMat(V.cAt(i-1),-V.cAt(i))};
	return {{Er[i],Ei[i]},cMat(V.cAt(i))};
}
fMat lm__::eigh(fMat inp) noexcept {
	assert(!inp.empty());
	assert(inp.M()>=2);

	fMat w(inp.M(),1);
	fMat work(4*inp.M(),1);
	int info;

	c_xsyev('N','U',inp.M(),inp.data(),inp.M(),w.data(),work.data(),work.M(),&info);

	assert(!info);
	return w;
}
fMat lm__::eigh(cMat inp) noexcept {
	assert(!inp.empty());
	assert(inp.M()>=2);

	fMat w(inp.M(),1);
	cMat work(4*inp.M(),1);
	int info;

	fMat rwork(3*inp.M()-2,1);
	c_xheev('N','U',inp.M(),inp.data(),inp.M(),w.data(),work.data(),work.M(),rwork.data(),&info);

	assert(!info);
	return w;
}
fehv lm__::eighv(fMat inp) noexcept {
	assert(!inp.empty());
	assert(inp.M()>=2);

	fMat w(inp.M(),1);
	fMat work(4*inp.M(),1);
	int info;

	c_xsyev('V','U',inp.M(),inp.data(),inp.M(),w.data(),work.data(),work.M(),&info);

	assert(!info);
	return {std::move(w),std::move(inp)};
}
cehv lm__::eighv(cMat inp) noexcept {
	assert(!inp.empty());
	assert(inp.M()>=2);

	fMat w(inp.M(),1);
	cMat work(4*inp.M(),1);
	int info;

	fMat rwork(3*inp.M()-2,1);
	c_xheev('V','U',inp.M(),inp.data(),inp.M(),w.data(),work.data(),work.M(),rwork.data(),&info);

	assert(!info);
	return {std::move(w),std::move(inp)};
}


// matrix generation functions
template<class MT>
MT lm__::ones(const size_t M) noexcept { return ones<MT>(M,M); }
template<>
fMat lm__::ones(const size_t M, const size_t N) noexcept {
	fMat res(M,N);
	std::fill(res.begin(),res.end(),RE__(1.0));
	return res;
}
template<>
cMat lm__::ones(const size_t M, const size_t N) noexcept {
	cMat res(M,N);
	std::fill(res.begin(),res.end(),CPX__(1.0,0.0));
	return res;
}
template fMat lm__::ones(const size_t M) noexcept;
template cMat lm__::ones(const size_t M) noexcept;
template<>
fMat lm__::ones(const lm_size& S) noexcept { return ones<fMat>(S.M,S.N); }
template<>
cMat lm__::ones(const lm_size& S) noexcept { return ones<cMat>(S.M,S.N); }
template<class MT>
MT lm__::zeros(const size_t M) noexcept { return zeros<MT>(M,M); }
template<>
fMat lm__::zeros(const size_t M, const size_t N) noexcept {
	fMat res(M,N);
	std::fill(res.begin(),res.end(),RE__(0.0));
	return res;
}
template<>
cMat lm__::zeros(const size_t M, const size_t N) noexcept {
	cMat res(M,N);
	std::fill(res.begin(),res.end(),CPX__(0.0,0.0));
	return res;
}
template fMat lm__::zeros(const size_t M) noexcept;
template cMat lm__::zeros(const size_t M) noexcept;
template<>
fMat lm__::zeros(const lm_size& S) noexcept { return zeros<fMat>(S.M,S.N); }
template<>
cMat lm__::zeros(const lm_size& S) noexcept { return zeros<cMat>(S.M,S.N); }
template<class MT>
MT lm__::eye(const size_t M) noexcept { return eye<MT>(M,M); }
template<>
fMat lm__::eye(const size_t M, const size_t N) noexcept {
	auto res = zeros<fMat>(M,N);
	for (auto i=res.dbegin(),e=res.dend(); i!=e; ++i)
		*i = RE__(1.0);
	return res;
}
template<>
cMat lm__::eye(const size_t M, const size_t N) noexcept {
	auto res = zeros<cMat>(M,N);
	for (auto i=res.dbegin(),e=res.dend(); i!=e; ++i)
		*i = CPX__(1.0,0.0);
	return res;
}
template fMat lm__::eye(const size_t M) noexcept;
template cMat lm__::eye(const size_t M) noexcept;
template<>
fMat lm__::eye(const lm_size& S) noexcept { return eye<fMat>(S.M,S.N); }
template<>
cMat lm__::eye(const lm_size& S) noexcept { return eye<cMat>(S.M,S.N); }
template<> 
fMat lm__::rId(const size_t N, const size_t i) noexcept {
	assert(i<N);
	auto res = zeros<fMat>(1,N); res[i]=RE__(1.0);
	return res;
}
template<> 
cMat lm__::rId(const size_t N, const size_t i) noexcept {
	assert(i<N);
	auto res = zeros<cMat>(1,N); res[i]=CPX__(1.0);
	return res;
}
template<>
fMat lm__::cId(const size_t M, const size_t i) noexcept {
	assert(i<M);
	auto res = zeros<fMat>(M,1); res[i]=RE__(1.0);
	return res;
}
template<>
cMat lm__::cId(const size_t M, const size_t i) noexcept {
	assert(i<M);
	auto res = zeros<cMat>(M,1); res[i]=CPX__(1.0);
	return res;
}
template<class MT>
MT lm__::rand(const size_t M) noexcept {
	return rand<MT>(M,M,RE__(0.0),RE__(1.0));
}
template<>
fMat lm__::rand(const size_t M, const size_t N, const RE__ rand_min, const RE__ rand_max) noexcept {
	fMat res(M,N);
	std::random_device rd;
#ifdef DOUBLE__
	std::mt19937_64 gen(rd());
#endif
#ifdef SINGLE__
	std::mt19937 gen(rd());
#endif
	std::uniform_real_distribution<RE__> dis(rand_min,rand_max);
	for (auto& i: res)
		i=dis(gen);
	return res;
}
template<>
cMat lm__::rand(const size_t M, const size_t N, const RE__ rand_min, const RE__ rand_max) noexcept {
	cMat res(M,N);
	std::random_device rd;
#ifdef DOUBLE__
	std::mt19937_64 gen(rd());
#endif
#ifdef SINGLE__
	std::mt19937 gen(rd());
#endif
	std::uniform_real_distribution<RE__> dis(rand_min,rand_max);
	for (auto& i: res)
		i={dis(gen),dis(gen)};
	return res;
}
template fMat lm__::rand(const size_t M) noexcept;
template cMat lm__::rand(const size_t M) noexcept;
template<>
fMat lm__::rand(const lm_size& S, const RE__ rand_min, const RE__ rand_max) noexcept {
	return rand<fMat>(S.M,S.N,rand_min,rand_max);
}
template<>
cMat lm__::rand(const lm_size& S, const RE__ rand_min, const RE__ rand_max) noexcept {
	return rand<cMat>(S.M,S.N,rand_min,rand_max);
}
template<>
fMat lm__::rand(const fMat& l, const fMat& u) noexcept {
	assert(msize(l)==msize(u));
	assert(all(l.leq(u)));
	
	fMat res(l.M(),l.N());
	std::random_device rd;
#ifdef DOUBLE__
	std::mt19937_64 gen(rd());
#endif
#ifdef SINGLE__
	std::mt19937 gen(rd());
#endif
	auto il=l.begin(), iu=u.begin();
	for (auto& i: res)
		i = std::uniform_real_distribution<RE__>(*il++,*iu++)(gen);
	return res;
}
template<>
cMat lm__::rand(const fMat& l, const fMat& u) noexcept {
	assert(msize(l)==msize(u));
	assert(all(l.leq(u)));
	
	cMat res(l.M(),l.N());
	std::random_device rd;
#ifdef DOUBLE__
	std::mt19937_64 gen(rd());
#endif
#ifdef SINGLE__
	std::mt19937 gen(rd());
#endif
	auto il=l.begin(), iu=u.begin();
	for (auto& i: res) {
		std::uniform_real_distribution<RE__> dis(*il++,*iu++);
		i = {dis(gen),dis(gen)};
	}
	return res;
}

template<class MT>
MT lm__::randi(const size_t M) noexcept { return randi<MT>(M,M,long(0),long(100)); }
template<>
fMat lm__::randi(const size_t M, const size_t N, const long rand_min, const long rand_max) noexcept {
	fMat res(M,N);

	std::random_device rd;
#ifdef DOUBLE__
	std::mt19937_64 gen(rd());
#endif
#ifdef SINGLE__
	std::mt19937 gen(rd());
#endif
	std::uniform_int_distribution<long> dis(rand_min,rand_max);
	for (auto& i: res)
		i=RE__(dis(gen));
	return res;
}
template<>
cMat lm__::randi(const size_t M, const size_t N, const long rand_min, const long rand_max) noexcept {
	cMat res(M,N);

	std::random_device rd;
#ifdef DOUBLE__
	std::mt19937_64 gen(rd());
#endif
#ifdef SINGLE__
	std::mt19937 gen(rd());
#endif
	std::uniform_int_distribution<long> dis(rand_min,rand_max);
	for (auto& i: res)
		i={RE__(dis(gen)),RE__(dis(gen))};
	return res;
}
template fMat lm__::randi(const size_t M) noexcept;
template cMat lm__::randi(const size_t M) noexcept;
template<>
fMat lm__::randi(const lm_size& S, const long rand_min, const long rand_max) noexcept {
	return randi<fMat>(S.M,S.N,rand_min,rand_max);
}
template<>
cMat lm__::randi(const lm_size& S, const long rand_min, const long rand_max) noexcept {
	return randi<cMat>(S.M,S.N,rand_min,rand_max);
}
template<>
fMat lm__::randi(const fMat& l, const fMat& u) noexcept {
	set_mtol(0.0);
	assert(msize(l)==msize(u));
	assert(l==(round(l)));
	assert(u==(round(u)));
	assert(all(l.leq(u)));
	reset_mtol();
	
	fMat res(l.M(),l.N());

	std::random_device rd;
#ifdef DOUBLE__
	std::mt19937_64 gen(rd());
#endif
#ifdef SINGLE__
	std::mt19937 gen(rd());
#endif
	auto il=l.begin(), iu=u.begin();
	for (auto& i: res)
		i = std::uniform_int_distribution<long>(long(*il++),long(*iu++))(gen);
	return res;
}
template<>
cMat lm__::randi(const fMat& l, const fMat& u) noexcept {
	set_mtol(0.0);
	assert(msize(l)==msize(u));
	assert(l==(round(l)));
	assert(u==(round(u)));
	assert(all(l.leq(u)));
	reset_mtol();
	
	cMat res(l.M(),l.N());

	std::random_device rd;
#ifdef DOUBLE__
	std::mt19937_64 gen(rd());
#endif
#ifdef SINGLE__
	std::mt19937 gen(rd());
#endif
	auto il=l.begin(), iu=u.begin();
	for (auto& i: res) {
		std::uniform_int_distribution<long> dis(long(*il++),long(*iu++));
		i = {RE__(dis(gen)),RE__(dis(gen))};
	}
	return res;
}
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::mcat(const lm_tMat<TT,FT,CT>& lhs, const lm_tMat<TT,FT,CT>& rhs) noexcept {
	assert(lhs.N()==rhs.N());
	if (lhs.empty()) return rhs;
	if (rhs.empty()) return lhs;

	lm_tMat<TT,FT,CT> res(lhs.M()+rhs.M(),rhs.N());
	for (size_t n=0; n!=lhs.N(); ++n) {
		const auto dst = res.data()+n*res.M();
		memcpy(dst,lhs.data()+n*lhs.M(),lhs.M()*sizeof(TT));
		memcpy(dst+lhs.M(),rhs.data()+n*rhs.M(),rhs.M()*sizeof(TT));
	}
	return res;
}
template fMat lm__::mcat(const fMat& lhs, const fMat& rhs) noexcept;	
template cMat lm__::mcat(const cMat& lhs, const cMat& rhs) noexcept;	

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::mcat(const lm_tMat<TT,FT,CT>& lhs, const lm_tRow<TT,FT,CT>& rhs) noexcept {
	assert(lhs.N()==rhs.N());
	if (lhs.empty()) return rhs;

	lm_tMat<TT,FT,CT> res(lhs.M()+1,lhs.N());
	auto ii = rhs.begin();
	for (size_t n=0; n!=lhs.N(); ++n) {
		const auto dst = res.data()+n*res.M();
		memcpy(dst,lhs.data()+n*lhs.M(),lhs.M()*sizeof(TT));
		*(dst+lhs.M()) = *ii;
		++ii;
	}
	return res;
}
template fMat lm__::mcat(const fMat& lhs, const fRow& rhs) noexcept;	
template cMat lm__::mcat(const cMat& lhs, const cRow& rhs) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::mcat(const lm_tRow<TT,FT,CT>& lhs, const lm_tMat<TT,FT,CT>& rhs) noexcept {
	assert(lhs.N()==rhs.N());
	if (rhs.empty()) return lhs;

	lm_tMat<TT,FT,CT> res(rhs.M()+1,lhs.N());
	auto ii = lhs.begin();
	for (size_t n=0; n!=lhs.N(); ++n) {
		const auto dst = res.data()+n*res.M();
		memcpy(dst+1,rhs.data()+n*rhs.M(),rhs.M()*sizeof(TT));
		*dst = *ii;
		++ii;
	}
	return res;	
}
template fMat lm__::mcat(const fRow& lhs, const fMat& rhs) noexcept;	
template cMat lm__::mcat(const cRow& lhs, const cMat& rhs) noexcept;	

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::mcat(const lm_tRow<TT,FT,CT>& lhs, const lm_tRow<TT,FT,CT>& rhs) noexcept {
	assert(lhs.N()==rhs.N());

	lm_tMat<TT,FT,CT> res(2,lhs.N());
	auto ii_lhs = lhs.begin(), ii_rhs = rhs.begin();
	for (auto i=res.begin(), e=res.end(); i!=e; ++ii_lhs, ++ii_rhs) {
		*i++ = *ii_lhs;
		*i++ = *ii_rhs;
	}
	return res;
}
template fMat lm__::mcat(const fRow& lhs, const fRow& rhs) noexcept;	
template cMat lm__::mcat(const cRow& lhs, const cRow& rhs) noexcept;	

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::ncat(const lm_tArray<TT,FT,CT>& lhs, const lm_tArray<TT,FT,CT>& rhs) noexcept {
	assert(lhs.M()==rhs.M());
	if (lhs.empty()) return rhs;
	if (rhs.empty()) return lhs;

	lm_tMat<TT,FT,CT> res(lhs.M(),0); res.reserve(rhs.N()+lhs.N());
	return res << lhs << rhs;
}
template fMat lm__::ncat(const fArray& lhs, const fArray& rhs) noexcept;
template cMat lm__::ncat(const cArray& lhs, const cArray& rhs) noexcept;


// orthogonalisation
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT>& lm__::gsorth(lm_tMat<TT,FT,CT>& inp) noexcept {
	assert(inp.square());
	assert(rank(inp)==inp.M());

	for (auto i=inp.cBegin(),e=inp.cEnd(); i!=e; ++i) {
		for (auto j=inp.ccBegin(); j!=i; ++j)
			sum(*i,*j,-dot(*j,*i));
		*i /= norm(*i);
	}
	return inp;
}
template fMat& lm__::gsorth(fMat& inp) noexcept;
template cMat& lm__::gsorth(cMat& inp) noexcept;
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::gsorth(lm_tMat<TT,FT,CT>&& inp) noexcept {
	return gsorth(inp);
}
template fMat lm__::gsorth(fMat&& inp) noexcept;
template cMat lm__::gsorth(cMat&& inp) noexcept;
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::complement(lm_tMat<TT,FT,CT> inp) noexcept {
	assert(inp.M()>=inp.N());
	assert(rank(inp)==inp.N());
	if (inp.empty() || inp.square()) return lm_tMat<TT,FT,CT>(inp.M(),0);
	
	// buffers
	lm_tMat<TT,FT,CT> tau(std::min(inp.M(),inp.N()),1);
	lm_tMat<TT,FT,CT> work(1,1);
	int info;

	// query lwork
	c_xgeqrf(inp.M(),inp.N(),inp.data(),inp.M(),tau.data(),work.data(),-1,&info);
	assert(!info);

	// resize work and compute qr decomposition
	work.resize(size_t(std::real(work[0])),1);
	c_xgeqrf(inp.M(),inp.N(),inp.data(),inp.M(),tau.data(),work.data(),work.size(),&info);
	assert(!info);

	// compute Q
	c_xungqr(inp.M(),inp.N(),tau.L(),inp.data(),inp.M(),tau.data(),work.data(),work.size(),&info);
	assert(!info);
	
	// find smallest entries to get ideal complement extension
	const size_t N_ = inp.N();
	inp.reserve(inp.M());
	auto S = nsum(abs(inp));
	while (inp.N()<inp.ccap()) {
		const auto i = mmin(S).pos[0]; S[i]=DBL_MAX;
		inp.push_back(cId<lm_tMat<TT,FT,CT>>(inp.M(),i));
	}

	// gsorth, return extended vectors
	return gsorth(inp).cGet(N_,inp.M()-N_);
}
template fMat lm__::complement(fMat inp) noexcept;
template cMat lm__::complement(cMat inp) noexcept;


// functionals
template<class TT, class FT, class CT>
TT lm__::sum(const lm_tArray<TT,FT,CT>& inp) noexcept {
	return std::accumulate(inp.begin(),inp.end(),TT(0.0));
}
template RE__ lm__::sum(const fArray& inp) noexcept;
template CPX__ lm__::sum(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
TT lm__::prod(const lm_tArray<TT,FT,CT>& inp) noexcept {
	return std::accumulate(inp.begin(),inp.end(),TT(1.0),std::multiplies<TT>());
}
template RE__ lm__::prod(const fArray& inp) noexcept;
template CPX__ lm__::prod(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
TT lm__::min(const lm_tArray<TT,FT,CT>& inp) noexcept {
	assert(!inp.empty());
	return *std::min_element(inp.begin(),inp.end(),
			[](const TT i, const TT j){return std::real(i)<std::real(j);});
}
template RE__ lm__::min(const fArray& inp) noexcept;
template CPX__ lm__::min(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
TT lm__::max(const lm_tArray<TT,FT,CT>& inp) noexcept {
	assert(!inp.empty());
	return *std::max_element(inp.begin(),inp.end(),
			[](const TT i, const TT j){return std::real(i)<std::real(j);});
}
template RE__ lm__::max(const fArray& inp) noexcept;
template CPX__ lm__::max(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
TT lm__::mean(const lm_tArray<TT,FT,CT>& inp) noexcept {
	assert(!inp.empty());
	return sum(inp)/FT(inp.L());
}
template RE__ lm__::mean(const fArray& inp) noexcept;
template CPX__ lm__::mean(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
FT lm__::stdd(const lm_tArray<TT,FT,CT>& inp) noexcept {
	return lm__::stdd(inp,mean(inp));
}
template RE__ lm__::stdd(const fArray& inp) noexcept;
template RE__ lm__::stdd(const cArray& inp) noexcept;
template<class TT, class FT, class CT>
FT lm__::stdd(const lm_tArray<TT,FT,CT>& inp, const TT& m) noexcept {
	if (inp.empty()) return FT(0.0);
	const FT res = std::accumulate(inp.cbegin(),inp.cend(),FT(0.0),
			[&m](const FT s, const TT& i) -> FT
			{ const FT tmp = std::abs(i-m); return s + tmp*tmp;});
	return std::sqrt(res/FT(inp.L()-1));
}
template RE__ lm__::stdd(const fArray& inp, const RE__& m) noexcept;
template RE__ lm__::stdd(const cArray& inp, const CPX__& m) noexcept;

RE__ lm__::normsq(const fArray& inp) noexcept {
	return dot(inp,inp);
}
RE__ lm__::normsq(const cArray& inp) noexcept {
	RE__ res(0.0);
	c_imItr ii = inp.begin();
	for (c_reItr ir=inp.begin(),er=inp.end(); ir!=er; ++ir,++ii)
		res += (*ir) * (*ir), res += (*ii) * (*ii);
	return res;
}

template<class TT, class FT, class CT>
FT lm__::norm(const lm_tArray<TT,FT,CT>& inp) noexcept {
	return std::sqrt(normsq(inp));
}
template RE__ lm__::norm(const fArray& inp) noexcept;
template RE__ lm__::norm(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
TT lm__::trace(const lm_tMat<TT,FT,CT>& inp) noexcept {
	assert(!inp.empty());
	assert(inp.square());

	TT res(0.0);
	for (auto i=inp.begin(),e=inp.end(); i<e; i+=inp.M()+1)
		res+=*i;
	return res;
}
template RE__ lm__::trace(const fMat& inp) noexcept;
template CPX__ lm__::trace(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
TT lm__::det(lm_tMat<TT,FT,CT> inp) noexcept {
	assert(!inp.empty());
	assert(inp.square());

	int* ipiv = new int [inp.N()]; int info;
	c_xgetrf(inp.M(),inp.M(),inp.data(),inp.M(),ipiv,&info);

	TT res(1.0); auto j=inp.begin();
	for (size_t i=0; i!=inp.N(); ++i) {
		if (i+1==(size_t)ipiv[i]) res *= *j;
		else res *= -*j;
		j+=inp.M()+1;
	}
	delete[] ipiv;
	return res;
}
template RE__ lm__::det(const fMat inp) noexcept;
template CPX__ lm__::det(const cMat inp) noexcept;


// dot products
RE__ lm__::dot(const fArray& lhs, const fArray& rhs) noexcept {
	assert(lhs.L()==rhs.L());
	return c_xdot(lhs.L(),lhs.data(),lhs.incr(),rhs.data(),rhs.incr());
}
CPX__ lm__::dot(const fArray& lhs, const cArray& rhs) noexcept {
	assert(lhs.L()==rhs.L());
	
	const c_reItr ir = rhs.begin();
	const c_imItr ii = rhs.begin();
	return { c_xdot(lhs.L(),lhs.data(),lhs.incr(),ir.data(),ir.incr()),
		 c_xdot(lhs.L(),lhs.data(),lhs.incr(),ii.data(),ii.incr()) };
}
CPX__ lm__::dot(const cArray& lhs, const fArray& rhs) noexcept {
	assert(lhs.L()==rhs.L());
	
	const c_reItr ir = lhs.begin();
	const c_imItr ii = lhs.begin();
	return { c_xdot(lhs.L(),ir.data(),ir.incr(),rhs.data(),rhs.incr()),
		-c_xdot(lhs.L(),ii.data(),ii.incr(),rhs.data(),rhs.incr()) };
}
CPX__ lm__::dot(const cArray& lhs, const cArray& rhs) noexcept {
	assert(lhs.L()==rhs.L());
	return c_xdotc(lhs.L(),lhs.data(),lhs.incr(),rhs.data(),rhs.incr());
}
RE__ lm__::dotu(const fArray& lhs, const fArray& rhs) noexcept {
	return dot(lhs,rhs);
}
CPX__ lm__::dotu(const fArray& lhs, const cArray& rhs) noexcept {
	return dot(lhs,rhs);
}
CPX__ lm__::dotu(const cArray& lhs, const fArray& rhs) noexcept {
	assert(lhs.L()==rhs.L());
	
	const c_reItr ir = lhs.begin();
	const c_imItr ii = lhs.begin();
	return { c_xdot(lhs.L(),ir.data(),ir.incr(),rhs.data(),rhs.incr()),
		 c_xdot(lhs.L(),ii.data(),ii.incr(),rhs.data(),rhs.incr()) };
}
CPX__ lm__::dotu(const cArray& lhs, const cArray& rhs) noexcept {
	assert(lhs.L()==rhs.L());
	return c_xdotu(lhs.L(),lhs.data(),lhs.incr(),rhs.data(),rhs.incr());
}


// vector sums
fRow& lm__::sum(fRow& lhs, const fRow& rhs, const RE__ f) noexcept {
	assert(lhs.L()==rhs.L());
	
	c_xaxpy(lhs.L(),f,rhs.data(),rhs.incr(),lhs.data(),lhs.incr());
	return lhs;
}
fRow& lm__::sum(fRow& lhs, const cRow& rhs, const RE__ f) noexcept {
	assert(lhs.L()==rhs.L());

	const c_reItr ir = rhs.begin();
	c_xaxpy(lhs.L(),f,ir.data(),ir.incr(),lhs.data(),lhs.incr());
	return lhs;
}
fRow& lm__::sum(fRow& lhs, const cRow& rhs, const CPX__ f) noexcept {
	assert(lhs.L()==rhs.L());
	
	const c_reItr ir = rhs.begin();
	const c_imItr ii = rhs.begin();
	c_xaxpy(lhs.L(),std::real(f),ir.data(),ir.incr(),lhs.data(),lhs.incr());
	c_xaxpy(lhs.L(),-std::imag(f),ii.data(),ii.incr(),lhs.data(),lhs.incr());
	return lhs;
}
cRow& lm__::sum(cRow& lhs, const fRow& rhs, const RE__ f) noexcept {
	assert(lhs.L()==rhs.L());

	const c_reItr ir = lhs.begin();
	c_xaxpy(lhs.L(),f,rhs.data(),rhs.incr(),ir.data(),ir.incr());
	return lhs;
}
cRow& lm__::sum(cRow& lhs, const fRow& rhs, const CPX__ f) noexcept {
	assert(lhs.L()==rhs.L());

	const c_reItr ir = lhs.begin();
	const c_imItr ii = lhs.begin();
	c_xaxpy(lhs.L(),std::real(f),rhs.data(),rhs.incr(),ir.data(),ir.incr());
	c_xaxpy(lhs.L(),std::imag(f),rhs.data(),rhs.incr(),ii.data(),ii.incr());
	return lhs;
}
cRow& lm__::sum(cRow& lhs, const cRow& rhs, const CPX__ f) noexcept {
	assert(lhs.L()==rhs.L());
	c_xaxpy(lhs.L(),f,rhs.data(),rhs.incr(),lhs.data(),lhs.incr());
	return lhs;
}
fCol& lm__::sum(fCol& lhs, const fCol& rhs, const RE__ f) noexcept {
	assert(lhs.L()==rhs.L());
	c_xaxpy(lhs.L(),f,rhs.data(),rhs.incr(),lhs.data(),lhs.incr());
	return lhs;
}
fCol& lm__::sum(fCol& lhs, const cCol& rhs, const RE__ f) noexcept {
	assert(lhs.L()==rhs.L());

	const c_reItr ir = rhs.begin();
	c_xaxpy(lhs.L(),f,ir.data(),ir.incr(),lhs.data(),lhs.incr());
	return lhs;
}
fCol& lm__::sum(fCol& lhs, const cCol& rhs, const CPX__ f) noexcept {
	assert(lhs.L()==rhs.L());

	const c_reItr ir = rhs.begin();
	const c_imItr ii = rhs.begin();
	c_xaxpy(lhs.L(),std::real(f),ir.data(),ir.incr(),lhs.data(),lhs.incr());
	c_xaxpy(lhs.L(),-std::imag(f),ii.data(),ii.incr(),lhs.data(),lhs.incr());
	return lhs;
}
cCol& lm__::sum(cCol& lhs, const fCol& rhs, const RE__ f) noexcept {
	assert(lhs.L()==rhs.L());

	const c_reItr ir = lhs.begin();
	c_xaxpy(lhs.L(),f,rhs.data(),rhs.incr(),ir.data(),ir.incr());
	return lhs;
}
cCol& lm__::sum(cCol& lhs, const fCol& rhs, const CPX__ f) noexcept {
	assert(lhs.L()==rhs.L());

	const c_reItr ir = lhs.begin();
	const c_imItr ii = lhs.begin();
	c_xaxpy(lhs.L(),std::real(f),rhs.data(),rhs.incr(),ir.data(),ir.incr());
	c_xaxpy(lhs.L(),std::imag(f),rhs.data(),rhs.incr(),ii.data(),ii.incr());
	return lhs;
}
cCol& lm__::sum(cCol& lhs, const cCol& rhs, const CPX__ f) noexcept {
	assert(lhs.L()==rhs.L());
	c_xaxpy(lhs.L(),f,rhs.data(),rhs.incr(),lhs.data(),lhs.incr());
	return lhs;
}
fRow&& lm__::sum(fRow&& lhs, const fRow& rhs, const RE__ f) noexcept {
	return std::move(sum(lhs,rhs,f));
}
fRow&& lm__::sum(fRow&& lhs, const cRow& rhs, const RE__ f) noexcept {
	return std::move(sum(lhs,rhs,f));
}
fRow&& lm__::sum(fRow&& lhs, const cRow& rhs, const CPX__ f) noexcept {
	return std::move(sum(lhs,rhs,f));
}
cRow&& lm__::sum(cRow&& lhs, const fRow& rhs, const RE__ f) noexcept {
	return std::move(sum(lhs,rhs,f));
}
cRow&& lm__::sum(cRow&& lhs, const fRow& rhs, const CPX__ f) noexcept {
	return std::move(sum(lhs,rhs,f));
}
cRow&& lm__::sum(cRow&& lhs, const cRow& rhs, const CPX__ f) noexcept {
	return std::move(sum(lhs,rhs,f));
}
fCol&& lm__::sum(fCol&& lhs, const fCol& rhs, const RE__ f) noexcept {
	return std::move(sum(lhs,rhs,f));
}
fCol&& lm__::sum(fCol&& lhs, const cCol& rhs, const RE__ f) noexcept {
	return std::move(sum(lhs,rhs,f));
}
fCol&& lm__::sum(fCol&& lhs, const cCol& rhs, const CPX__ f) noexcept {
	return std::move(sum(lhs,rhs,f));
}
cCol&& lm__::sum(cCol&& lhs, const fCol& rhs, const RE__ f) noexcept {
	return std::move(sum(lhs,rhs,f));
}
cCol&& lm__::sum(cCol&& lhs, const fCol& rhs, const CPX__ f) noexcept {
	return std::move(sum(lhs,rhs,f));
}
cCol&& lm__::sum(cCol&& lhs, const cCol& rhs, const CPX__ f) noexcept {
	return std::move(sum(lhs,rhs,f));
}


// R^3 only
fMat lm__::cross(const fArray& lhs, const fArray& rhs) noexcept {
	assert(lhs.L()==3);
	assert(rhs.L()==3);

	const size_t M = (lhs.N()==3&&rhs.N()==3) ? 1: 3;
	const size_t N = (lhs.N()==3&&rhs.N()==3) ? 3: 1;

	return fMat({lhs[1]*rhs[2] - lhs[2]*rhs[1],
		     lhs[2]*rhs[0] - lhs[0]*rhs[2],
		     lhs[0]*rhs[1] - lhs[1]*rhs[0]},M,N);
}
fMat lm__::crossn(const fArray& lhs, const fArray& rhs) noexcept {
	auto res = cross(lhs,rhs);
	return res/norm(res);
}
fMat lm__::getR(const RE__ phi, const fArray& n) noexcept {
	assert(n.L()==3);
	assert(lm__::ops::eq(norm(n),1.0));

	const fMat W({0.0, n[2],-n[1],
		    -n[2],  0.0, n[0],
		     n[1], -n[0], 0.0},3,3);

	return eye<fMat>(3,3) + std::sin(phi)*W + (RE__(1.0)-std::cos(phi))*W.prod(W);
}
fMat lm__::getR(const fArray& v, const fArray& w) noexcept {
	assert(v.L()==3);
	assert(w.L()==3);
	return getR(std::acos(dot(v,w)/norm(v)/norm(w)),crossn(v,w));
}


// matlab like conversion functions
template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::msum(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return lm_tMat<TT,FT,CT>(0,inp.N());

	lm_tMat<TT,FT,CT> res(1,inp.N());
	auto ii = inp.cBegin();
	for (auto& i: res)
		i = sum(*ii++);
	return res;
}
template fMat lm__::msum(const fMat& inp) noexcept;
template cMat lm__::msum(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::nsum(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return lm_tMat<TT,FT,CT>(inp.M(),0);

	lm_tMat<TT,FT,CT> res(inp.M(),1);
	auto ii = inp.rBegin();
	for (auto& i: res)
		i = sum(*ii++);
	return res;
}
template fMat lm__::nsum(const fMat& inp) noexcept;
template cMat lm__::nsum(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::mprod(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return lm_tMat<TT,FT,CT>(0,inp.N());

	lm_tMat<TT,FT,CT> res(1,inp.N());
	auto ii = inp.ccBegin();
	for (auto& i: res)
		i = prod(*ii++);
	return res;
}
template fMat lm__::mprod(const fMat& inp) noexcept;
template cMat lm__::mprod(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::nprod(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return lm_tMat<TT,FT,CT>(inp.M(),0);

	lm_tMat<TT,FT,CT> res(inp.M(),1);
	auto ii = inp.crBegin();
	for (auto& i: res)
		i = prod(*ii++);
	return res;
}
template fMat lm__::nprod(const fMat& inp) noexcept;
template cMat lm__::nprod(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
mnp<TT,FT,CT> lm__::mmin(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return {lm_tMat<TT,FT,CT>(0,inp.N()),{}};

	mnp<TT,FT,CT> res = {lm_tMat<TT,FT,CT>(1,inp.N()),{}};
	res.pos.reserve(inp.N());

	auto ri=res.mat.begin();
	for (auto ii=inp.cBegin(),ie=inp.cEnd(); ii!=ie; ++ii,++ri) {
		const auto p = std::min_element(ii->begin(),ii->end(),
			[](const TT i, const TT j){return ops::lt_s(i,j);});
		*ri = *p; res.pos.push_back(p-ii->begin());
	}

	return res;
}
template mnp<RE__,RE__,CPX__> lm__::mmin(const fMat& inp) noexcept;
template mnp<CPX__,RE__,CPX__> lm__::mmin(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
mnp<TT,FT,CT> lm__::nmin(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return {lm_tMat<TT,FT,CT>(inp.M(),0),{}};

	mnp<TT,FT,CT> res = {lm_tMat<TT,FT,CT>(inp.M(),1),{}};
	res.pos.reserve(inp.M());

	auto ri=res.mat.begin();
	for (auto ii=inp.rBegin(),ie=inp.rEnd(); ii!=ie; ++ii,++ri) {
		const auto p = std::min_element(ii->begin(),ii->end(),
			[](const TT i, const TT j){return ops::lt_s(i,j);});
		*ri = *p; res.pos.push_back(p-ii->begin());
	}

	return res;
}
template mnp<RE__,RE__,CPX__> lm__::nmin(const fMat& inp) noexcept;
template mnp<CPX__,RE__,CPX__> lm__::nmin(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
mnp<TT,FT,CT> lm__::mmax(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return {lm_tMat<TT,FT,CT>(0,inp.N()),{}};

	mnp<TT,FT,CT> res = {lm_tMat<TT,FT,CT>(1,inp.N()),{}};
	res.pos.reserve(inp.N());

	auto ri=res.mat.begin();
	for (auto ii=inp.cBegin(),ie=inp.cEnd(); ii!=ie; ++ii,++ri) {
		const auto p = std::max_element(ii->begin(),ii->end(),
			[](const TT i, const TT j){return ops::lt_s(i,j);});
		*ri = *p; res.pos.push_back(p-ii->begin());
	}

	return res;
}
template mnp<RE__,RE__,CPX__> lm__::mmax(const fMat& inp) noexcept;
template mnp<CPX__,RE__,CPX__> lm__::mmax(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
mnp<TT,FT,CT> lm__::nmax(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return {lm_tMat<TT,FT,CT>(inp.M(),0),{}};

	mnp<TT,FT,CT> res = {lm_tMat<TT,FT,CT>(inp.M(),1),{}};
	res.pos.reserve(inp.M());

	auto ri=res.mat.begin();
	for (auto ii=inp.rBegin(),ie=inp.rEnd(); ii!=ie; ++ii,++ri) {
		const auto p = std::max_element(ii->begin(),ii->end(),
			[](const TT i, const TT j){return ops::lt_s(i,j);});
		*ri = *p; res.pos.push_back(p-ii->begin());
	}

	return res;
}
template mnp<RE__,RE__,CPX__> lm__::nmax(const fMat& inp) noexcept;
template mnp<CPX__,RE__,CPX__> lm__::nmax(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::mmean(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return lm_tMat<TT,FT,CT>(0,inp.N());
	return msum(inp)/FT(inp.M());
}
template fMat lm__::mmean(const fMat& inp) noexcept;
template cMat lm__::mmean(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::nmean(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return lm_tMat<TT,FT,CT>(inp.M(),0);
	return nsum(inp)/FT(inp.N());
}
template fMat lm__::nmean(const fMat& inp) noexcept;
template cMat lm__::nmean(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
fMat lm__::mnormsq(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return fMat(0,inp.N());
	fMat res(1,inp.N());
	std::transform(inp.ccBegin(),inp.ccEnd(),res.begin(),
		[](const lm_tCol<TT,FT,CT>& i){return normsq(i);});
	return res;
}
template fMat lm__::mnormsq(const fMat& inp) noexcept;
template fMat lm__::mnormsq(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
fMat lm__::nnormsq(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return fMat(inp.M(),0);
	fMat res(inp.M(),1);
	std::transform(inp.crBegin(),inp.crEnd(),res.begin(),
		[](const lm_tRow<TT,FT,CT>& i){return normsq(i);});
	return res;
}
template fMat lm__::nnormsq(const fMat& inp) noexcept;
template fMat lm__::nnormsq(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
fMat lm__::mnorm(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return fMat(0,inp.N());
	fMat res(1,inp.N());
	std::transform(inp.ccBegin(),inp.ccEnd(),res.begin(),
		[](const lm_tCol<TT,FT,CT>& i){return norm(i);});
	return res;
}
template fMat lm__::mnorm(const fMat& inp) noexcept;
template fMat lm__::mnorm(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
fMat lm__::nnorm(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return fMat(inp.M(),0);
	fMat res(inp.M(),1);
	std::transform(inp.crBegin(),inp.crEnd(),res.begin(),
		[](const lm_tRow<TT,FT,CT>& i){return norm(i);});
	return res;
}
template fMat lm__::nnorm(const fMat& inp) noexcept;
template fMat lm__::nnorm(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::mdot(const lm_tMat<TT,FT,CT>& lhs, const fMat& rhs) noexcept {
	assert(msize(lhs)==msize(rhs));
	if (lhs.empty()) return lm_tMat<TT,FT,CT>(0,lhs.N());
	
	lm_tMat<TT,FT,CT> res(1,lhs.N());
	std::transform(lhs.ccBegin(),lhs.ccEnd(),rhs.ccBegin(),res.begin(),
		[](const lm_tCol<TT,FT,CT>& i, const fCol& j){return dot(i,j);});
	return res;
}
template fMat lm__::mdot(const fMat& lhs, const fMat& rhs) noexcept;
template cMat lm__::mdot(const cMat& lhs, const fMat& rhs) noexcept;

template<class TT, class FT, class CT>
cMat lm__::mdot(const lm_tMat<TT,FT,CT>& lhs, const cMat& rhs) noexcept {
	assert(msize(lhs)==msize(rhs));
	if (lhs.empty()) return cMat(0,lhs.N());

	cMat res(1,lhs.N());
	std::transform(lhs.ccBegin(),lhs.ccEnd(),rhs.ccBegin(),res.begin(),
		[](const lm_tCol<TT,FT,CT>& i, const cCol& j){return dot(i,j);});
	return res;
}
template cMat lm__::mdot(const fMat& lhs, const cMat& rhs) noexcept;
template cMat lm__::mdot(const cMat& lhs, const cMat& rhs) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::ndot(const lm_tMat<TT,FT,CT>& lhs, const fMat& rhs) noexcept {
	assert(msize(lhs)==msize(rhs));
	if (lhs.empty()) return lm_tMat<TT,FT,CT>(lhs.M(),0);

	lm_tMat<TT,FT,CT> res(lhs.M(),1);
	std::transform(lhs.crBegin(),lhs.crEnd(),rhs.crBegin(),res.begin(),
		[](const lm_tRow<TT,FT,CT>& i, const fRow& j){return dot(i,j);});
	return res;
}
template fMat lm__::ndot(const fMat& lhs, const fMat& rhs) noexcept;
template cMat lm__::ndot(const cMat& lhs, const fMat& rhs) noexcept;

template<class TT, class FT, class CT>
cMat lm__::ndot(const lm_tMat<TT,FT,CT>& lhs, const cMat& rhs) noexcept {
	assert(msize(lhs)==msize(rhs));
	if (lhs.empty()) return cMat(lhs.M(),0);

	cMat res(lhs.M(),1);
	std::transform(lhs.crBegin(),lhs.crEnd(),rhs.crBegin(),res.begin(),
		[](const lm_tRow<TT,FT,CT>& i, const cRow& j){return dot(i,j);});
	return res;
}
template cMat lm__::ndot(const fMat& lhs, const cMat& rhs) noexcept;
template cMat lm__::ndot(const cMat& lhs, const cMat& rhs) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::mdotu(const lm_tMat<TT,FT,CT>& lhs, const fMat& rhs) noexcept {
	assert(msize(lhs)==msize(rhs));
	if (lhs.empty()) return lm_tMat<TT,FT,CT>(0,lhs.N());

	lm_tMat<TT,FT,CT> res(1,lhs.N());
	std::transform(lhs.ccBegin(),lhs.ccEnd(),rhs.ccBegin(),res.begin(),
		[](const lm_tCol<TT,FT,CT>& i, const fCol& j){return dotu(i,j);});
	return res;
}
template fMat lm__::mdotu(const fMat& lhs, const fMat& rhs) noexcept;
template cMat lm__::mdotu(const cMat& lhs, const fMat& rhs) noexcept;

template<class TT, class FT, class CT>
cMat lm__::mdotu(const lm_tMat<TT,FT,CT>& lhs, const cMat& rhs) noexcept {
	assert(msize(lhs)==msize(rhs));
	if (lhs.empty()) return cMat(0,lhs.N());
	
	cMat res(1,lhs.N());
	std::transform(lhs.ccBegin(),lhs.ccEnd(),rhs.ccBegin(),res.begin(),
		[](const lm_tCol<TT,FT,CT>& i, const cCol& j){return dotu(i,j);});
	return res;
}
template cMat lm__::mdotu(const fMat& lhs, const cMat& rhs) noexcept;
template cMat lm__::mdotu(const cMat& lhs, const cMat& rhs) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::ndotu(const lm_tMat<TT,FT,CT>& lhs, const fMat& rhs) noexcept {
	assert(msize(lhs)==msize(rhs));
	if (lhs.empty()) return lm_tMat<TT,FT,CT>(lhs.M(),0);

	lm_tMat<TT,FT,CT> res(lhs.M(),1);
	std::transform(lhs.crBegin(),lhs.crEnd(),rhs.crBegin(),res.begin(),
		[](const lm_tRow<TT,FT,CT>& i, const fRow& j){return dotu(i,j);});
	return res;
}
template fMat lm__::ndotu(const fMat& lhs, const fMat& rhs) noexcept;
template cMat lm__::ndotu(const cMat& lhs, const fMat& rhs) noexcept;

template<class TT, class FT, class CT>
cMat lm__::ndotu(const lm_tMat<TT,FT,CT>& lhs, const cMat& rhs) noexcept {
	assert(msize(lhs)==msize(rhs));
	if (lhs.empty()) return cMat(lhs.M(),0);

	cMat res(lhs.M(),1);
	std::transform(lhs.crBegin(),lhs.crEnd(),rhs.crBegin(),res.begin(),
		[](const lm_tRow<TT,FT,CT>& i, const cRow& j){return dotu(i,j);});
	return res;
}
template cMat lm__::ndotu(const fMat& lhs, const cMat& rhs) noexcept;
template cMat lm__::ndotu(const cMat& lhs, const cMat& rhs) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::inv(lm_tMat<TT,FT,CT> inp) noexcept {
	return inp.inv();
}
template fMat lm__::inv(fMat inp) noexcept;
template cMat lm__::inv(cMat inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::diag(const lm_tArray<TT,FT,CT>& inp, const size_t os, const bool m) noexcept {
	if (inp.empty()) return lm_tMat<TT,FT,CT>(0,0);
	if (inp.M()==1 || inp.N()==1) {
		lm_tMat<TT,FT,CT> res(inp.L()+os);
		
		std::fill(res.begin(),res.end(),TT(0.0));
		std::copy(inp.cbegin(),inp.cend(),res.dbegin(os,m));
		return res;
	}

	const auto i=inp.cdbegin(os,m), e=inp.cdend(os,m);
	lm_tMat<TT,FT,CT> res(e-i,1);
	std::copy(i,e,res.begin());
	
	return res;
}
template fMat lm__::diag(const fArray& inp, const size_t os, const bool m) noexcept;
template cMat lm__::diag(const cArray& inp, const size_t os, const bool m) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::R(const lm_tArray<TT,FT,CT>& inp) noexcept {
	return lm_tMat<TT,FT,CT>(inp).R();
}
template fMat lm__::R(const fArray& inp) noexcept;
template cMat lm__::R(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::C(const lm_tArray<TT,FT,CT>& inp) noexcept {
	return lm_tMat<TT,FT,CT>(inp).C();
}
template fMat lm__::C(const fArray& inp) noexcept;
template cMat lm__::C(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
lm_tMat<TT,FT,CT> lm__::T(const lm_tArray<TT,FT,CT>& inp) noexcept {
	return lm_tMat<TT,FT,CT>(inp).T();
}
template fMat lm__::T(const fArray& inp) noexcept;
template cMat lm__::T(const cArray& inp) noexcept;	

template<class TT, class FT, class CT>
std::vector<Ipair> lm__::Iz(const lm_tArray<TT,FT,CT>& inp) noexcept {
	std::vector<Ipair> res; res.reserve(inp.L());
	for (size_t n=0; n!=inp.N(); ++n)
		for (size_t m=0; m!=inp.M(); ++m)
			if (lm__::ops::z(inp(m,n)))
				res.push_back({m,n});
	res.shrink_to_fit();
	return res;
}
template std::vector<Ipair> lm__::Iz(const fArray& inp) noexcept;
template std::vector<Ipair> lm__::Iz(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
std::vector<Ipair> lm__::Inz(const lm_tArray<TT,FT,CT>& inp) noexcept {
	std::vector<Ipair> res; res.reserve(inp.L());
	for (size_t n=0; n!=inp.N(); ++n)
		for (size_t m=0; m!=inp.M(); ++m)
			if (lm__::ops::nz(inp(m,n)))
				res.push_back({m,n});
	res.shrink_to_fit();
	return res;
}
template std::vector<Ipair> lm__::Inz(const fArray& inp) noexcept;
template std::vector<Ipair> lm__::Inz(const cArray& inp) noexcept;

fMat lm__::conj(const fArray& inp) noexcept {
	return inp;
}
cMat lm__::conj(const cArray& inp) noexcept {
	cMat res(msize(inp));
	std::transform(inp.cbegin(),inp.cend(),res.begin(),
		[](const CPX__& i){return std::conj(i);});
	return res;
}
fMat lm__::real(const fArray& inp) noexcept {
	return inp;
}
fMat lm__::real(const cArray& inp) noexcept {
	return fMat(inp);
}
fMat lm__::imag(const fArray& inp) noexcept {
	return zeros<fMat>(inp.M(),inp.N());
}
fMat lm__::imag(const cArray& inp) noexcept {
	fMat res(msize(inp));

	c_imItr ir = inp.begin();
	for (auto& i: res)
		i = *ir, ir+=inp.incr();
	return res;
}


// information functions
template<class TT, class FT, class CT>
fMat lm__::isnan(const lm_tArray<TT,FT,CT>& inp) noexcept {
	fMat res(msize(inp));	
	std::transform(inp.cbegin(),inp.cend(),res.begin(),
		[](const TT& i){return std::isnan(i) ? FT(1.0): FT(0.0);});
	return res;
}
template fMat lm__::isnan(const fArray& inp) noexcept;
template fMat lm__::isnan(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
size_t lm__::length(const lm_tArray<TT,FT,CT>& inp) noexcept {
	return inp.M()>inp.N() ? inp.M(): inp.N();
}
template size_t lm__::length(const fArray& inp) noexcept;
template size_t lm__::length(const cArray& inp) noexcept;

size_t lm__::rank(fMat inp, const RE__ min_sv) noexcept {
	if (inp.empty()) return 0;

	const size_t min_rc = inp.M()<inp.N() ? inp.M(): inp.N();
	RE__* s = new RE__ [min_rc];

	// compute ideal lwork
	RE__ tmp; int info;
	c_xgesvd('N','N',inp.M(),inp.N(),nullptr,inp.M(),nullptr,nullptr,
			inp.M(),nullptr,min_rc,&tmp,-1,&info);
	size_t lwork=tmp;

	// compute singular values
	RE__* work = new RE__ [lwork];
	c_xgesvd('N','N',inp.M(),inp.N(),inp.data(),inp.M(),s,nullptr,
			inp.M(),nullptr,min_rc,work,lwork,&info);

	size_t res=0;
	for (size_t i=0; i!=min_rc; ++i)
		res+=(std::abs(s[i])>min_sv);
	delete[] work;
	delete[] s;

	return res;
}
size_t lm__::rank(cMat inp, const RE__ min_sv) noexcept {
	if (inp.empty()) return 0;

	const size_t min_rc = inp.M()<inp.N() ? inp.M(): inp.N();
	RE__* s = new RE__ [min_rc];

	// compute ideal lwork
	CPX__ tmp; int info;
	c_xgesvd('N','N',inp.M(),inp.N(),nullptr,inp.M(),nullptr,nullptr,
			inp.M(),nullptr,min_rc,&tmp,-1,nullptr,&info);
	size_t lwork=std::real(tmp);

	// compute singular values
	CPX__* work = new CPX__ [lwork];
	RE__* rwork = new RE__ [5*min_rc];
	c_xgesvd('N','N',inp.M(),inp.N(),inp.data(),inp.M(),s,nullptr,
			inp.M(),nullptr,min_rc,work,lwork,rwork,&info);
	delete[] rwork;

	size_t res=0;
	for (size_t i=0; i<min_rc; ++i)
		res+=(std::abs(s[i])>min_sv);
	delete[] work;
	delete[] s;

	return res;
}


// matlab like comparison functions
template<class TT, class FT, class CT>
bool lm__::all(const lm_tArray<TT,FT,CT>& inp) noexcept {
	return std::all_of(inp.begin(), inp.end(), [](const TT& i){return ops::nz(i);});
}
template bool lm__::all(const fArray& inp) noexcept;
template bool lm__::all(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
fMat lm__::mall(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return fMat(0,inp.N());
	fMat res(1,inp.N());
	std::transform(inp.ccBegin(),inp.ccEnd(),res.begin(),
		[](const lm_tCol<TT,FT,CT>& i){return all(i);});
	return res;
}
template fMat lm__::mall(const fMat& inp) noexcept;
template fMat lm__::mall(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
fMat lm__::nall(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return fMat(inp.M(),0);
	fMat res(inp.M(),1);
	std::transform(inp.crBegin(),inp.crEnd(),res.begin(),
		[](const lm_tRow<TT,FT,CT>& i){return all(i);});
	return res;
}
template fMat lm__::nall(const fMat& inp) noexcept;
template fMat lm__::nall(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
bool lm__::any(const lm_tArray<TT,FT,CT>& inp) noexcept {
	return std::any_of(inp.begin(), inp.end(), [](const TT& i){return ops::nz(i);});
}
template bool lm__::any(const fArray& inp) noexcept;
template bool lm__::any(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
fMat lm__::many(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return fMat(0,inp.N());
	fMat res(1,inp.N());
	std::transform(inp.ccBegin(),inp.ccEnd(),res.begin(),
		[](const lm_tCol<TT,FT,CT>& i){return any(i);});
	return res;
}
template fMat lm__::many(const fMat& inp) noexcept;
template fMat lm__::many(const cMat& inp) noexcept;

template<class TT, class FT, class CT>
fMat lm__::nany(const lm_tMat<TT,FT,CT>& inp) noexcept {
	if (inp.empty()) return fMat(inp.M(),0);
	fMat res(inp.M(),1);
	std::transform(inp.crBegin(),inp.crEnd(),res.begin(),
		[](const lm_tRow<TT,FT,CT>& i){return any(i);});
	return res;
}
template fMat lm__::nany(const fMat& inp) noexcept;
template fMat lm__::nany(const cMat& inp) noexcept;	

template<class TT, class FT, class CT>
std::vector<size_t> lm__::find(const lm_tArray<TT,FT,CT>& inp) noexcept {
	std::vector<size_t> res; res.reserve(inp.L());

	const auto j=inp.begin();
	for (auto i=inp.begin(),e=inp.end(); i!=e; ++i)
		if (ops::nz(*i))
			res.push_back(i-j);
	res.shrink_to_fit();
	
	return res;
}
template std::vector<size_t> lm__::find(const fArray& inp) noexcept;
template std::vector<size_t> lm__::find(const cArray& inp) noexcept;

template<class TT, class FT, class CT>
size_t lm__::nnz(const lm_tArray<TT,FT,CT>& inp) noexcept {
	return std::accumulate(inp.cbegin(),inp.cend(),size_t(0),
		[](const size_t s, const TT& i)->size_t{return s+lm__::ops::nz(i);});
}
template size_t lm__::nnz(const fArray& inp) noexcept;
template size_t lm__::nnz(const cArray& inp) noexcept;
