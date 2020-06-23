// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "ll_fn.h"
#include "ll_defs.h"
#include "aux_sort.h"
#include "lm_fn.h"
#include "lm_lambda.h"
#include "ll_cell.h"
#include "ll_io.h"
#include "ll_hio.h"
#include "ll_mesh.h"
#include <list>
#include <cassert>
#include <string>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <cfloat>

using namespace ll__;


// next neighbour matrix
const fMat ll__::NN({
		-1.0,-1.0,-1.0,
		-1.0,-1.0, 0.0,
		-1.0,-1.0, 1.0,
		-1.0, 0.0,-1.0,
		-1.0, 0.0, 0.0,
		-1.0, 0.0, 1.0,
		-1.0, 1.0,-1.0,
		-1.0, 1.0, 0.0,
		-1.0, 1.0, 1.0,
		 0.0,-1.0,-1.0,
		 0.0,-1.0, 0.0,
		 0.0,-1.0, 1.0,
		 0.0, 0.0,-1.0,
		 0.0, 0.0, 0.0,
		 0.0, 0.0, 1.0,
		 0.0, 1.0,-1.0,
		 0.0, 1.0, 0.0,
		 0.0, 1.0, 1.0,
		 1.0,-1.0,-1.0,
		 1.0,-1.0, 0.0,
		 1.0,-1.0, 1.0,
		 1.0, 0.0,-1.0,
		 1.0, 0.0, 0.0,
		 1.0, 0.0, 1.0,
		 1.0, 1.0,-1.0,
		 1.0, 1.0, 0.0,
		 1.0, 1.0, 1.0},DIM__);


// basis finding
fMat ll__::findBasisExpansion(const fMat& B, const fMat& V, const double lim, const double tol) {
	assert(ll__::volnorm(B)>mtol());
	assert(ll__::volnorm(V)>mtol());
	assert(lim>=1.0);
	assert(tol>=0.0);

	auto findDenominator = [lim,tol](double f) -> double {
		f=std::fmod(std::abs(f),1.0);
		double d=1.0;
		while (f>tol) {
			f=1.0/f; d*=f; f=std::fmod(f,1.0);
			if (d>lim) throw(0);
		}
		return std::round(d);
	};
	auto makeIntegers = [&findDenominator](fCol& inp) -> fCol& {
		inp/=max(abs(inp));
		for (auto i: inp)
			inp*=findDenominator(i);
		return inp;			
	};

	auto Vb = B.leftDivide(V);
	try {
		for (auto i=Vb.cBegin(),e=Vb.cEnd(); i!=e; ++i)
			makeIntegers(*i);
	} catch(...) {
		throw std::runtime_error("no valid basis found, bad input");
	}
	return round(Vb);
}
fMat ll__::orthogonalize(fMat B, const rv& r) noexcept {
	const auto i = inds(r);
	if (i.empty()) return B;

	const auto m = mnorm(B);
	const auto V = i.size()<B.M() ? complement(B.get({},ninds(r))): gsorth(B);
	for (size_t j=0; j!=i.size(); ++j)
		B.cAt(i[j]) = V.cAt(j)*m[i[j]];

	return B;
}


// linear path and mesh generation
p_p ll__::genPath(const fMat& P, std::vector<size_t> Nps) {
	if (P.empty()) return {fMat(P.M(),P.N()),{}};
	if (P.N()==1) return {P,{0.0}};

	if (Nps.size()!=P.N()-1)
		throw(std::invalid_argument(std::to_string(P.N())+" points provided but only "+
			std::to_string(Nps.size())+" number of steps, need "+std::to_string(P.N()-1)));
	if (std::any_of(Nps.begin(), Nps.end(), [](size_t i){return i<2;}))
		throw(std::invalid_argument("found number of steps < 2, need all steps >= 2"));
	
	fMat path(P.M(),std::accumulate(Nps.begin(),Nps.end(),0));
	fMat pos(1,path.N());
	double l=0.0;

	auto i=pos.begin(); auto j=path.cBegin(); auto p=P.ccBegin();
	for (const auto Np: Nps) {
		const auto v = *(p+1)-*p;
		const auto nv = lm__::norm(v);
		for (size_t n=0; n!=Np; ++n) {
			const double f = double(n)/double(Np-1);
			*j = *p+f*v;
			*i = l+f*nv;
			++i,++j;
		}
		++p; l+=nv;
	}

	return {std::move(path),std::move(pos),std::move(Nps)};
}
p_p ll__::genPath(fMat P, const size_t Np, const fMat& B) {
	if (P.empty() || !Np) return {fMat(P.M(),0),{}};
	if (P.N()==1) return {P,{0.0}};
	if (!B.empty()) P = B.prod(P);
	
	// remove consecutive points
	for (auto i=P.ccBegin(); i!=P.ccEnd()-1; ++i)
		if (*i==*(i+1)) P.cRm(i);

	double l=0.0;
	for (auto i=P.ccBegin(),ie=P.ccEnd()-1; i!=ie; ++i)
		l+=norm(*(i+1) - *i);
	
	std::vector<size_t> Nps; Nps.reserve(P.N()-1);
	double pl=0.0; size_t st=0;
	for (auto i=P.ccBegin(),ie=P.ccEnd()-1; i!=ie; ++i) {
		pl+=norm(*(i+1) - *i);
		
		size_t pos = size_t(pl/l*double(Np));
		Nps.push_back(pos-st);
		
		if (Nps.back()<2)
			throw(std::invalid_argument(
			"found number of steps < 2, need all >= 2, set more steps or less dense points"));
		st = pos;
	}
	
	return genPath(B.empty() ? P: B.leftDivide(P),std::move(Nps));
}


// metrics and clustering
fMat ll__::genGrid(const fMat& bounds, std::vector<size_t> maj) noexcept {
	assert(bounds==lm__::round(bounds));
	assert(bounds.M()==maj.size());
	assert(bounds.N()==2);
	assert(all(bounds.cFront().leq(bounds.cBack())));

	std::vector<size_t> D(bounds.M());
	for (size_t i=0; i!=D.size(); ++i)
		D[i] = (size_t)std::round(bounds(i,1)-bounds(i,0))+1;

	return std::move(genMesh(bounds,std::move(maj),std::move(D)).base_);
}
fMat ll__::genGrid(const fMat& bounds) noexcept {
	std::vector<size_t> maj; maj.reserve(bounds.M());
	for (size_t d=0; d!=bounds.M(); ++d)
		maj.push_back(bounds.M()-1-d);
	return genGrid(bounds,std::move(maj));
}
fMat ll__::genNNmat(const rv& r, std::vector<size_t> maj) noexcept {
	if (r.empty()) return {};

	fMat bounds(r.size(),0); bounds.reserve(2);
	bounds.push_back(-fMat(!r)); bounds.push_back(-bounds.cBack());

	return genGrid(bounds,std::move(maj));
}
fMat ll__::genNNmat(const rv& r) noexcept {
	std::vector<size_t> maj; maj.reserve(r.size());
	for (size_t d=0; d!=r.size(); ++d)
		maj.push_back(r.size()-1-d);
	return genNNmat(r,std::move(maj));
}
double ll__::dist(const fMat& B, const fArray& p1, const fArray& p2, fMat NN) noexcept {
	assert(all(p1.geq(0.0)) && all(p1.lt(1.0)));
	assert(all(p2.geq(0.0)) && all(p2.lt(1.0)));
	assert(p1.M()==DIM__ && p1.N()==1);
	assert(p2.M()==DIM__ && p2.N()==1);
	assert(B.square() && B.M()==DIM__);

	return min(mnorm(B.prod(cadd(NN,p2-p1))));
}
double ll__::dist(const fMat& B, const fArray& p1, const fArray& p2, const double f, fMat NN) noexcept {
	assert(all(p1.geq(0.0)) && all(p1.lt(1.0)));
	assert(all(p2.geq(0.0)) && all(p2.lt(1.0)));
	assert(p1.M()==DIM__ && p1.N()==1);
	assert(p2.M()==DIM__ && p2.N()==1);
	assert(B.square() && B.M()==DIM__);
	
	auto n = mnorm(B.prod(cadd(NN,p2-p1)));
	if (n.L()==1) return n[0];
	
	// sort and find largest acceptable entry
	std::sort(n.begin(),n.end());
	auto s=n.begin(), m=n.begin()+1, e=n.end();
	if (f<0.0) { // factor mode
		for (; m<e; ++s, ++m)
			if (-(*s)*f<(*m)) return *s;
	} else { // direct cutoff mode
		for (; m<e; ++s, ++m)
			if (f<(*m)) return *s;
	}
	return *s;
}
d_b ll__::distb(const fMat& B, const fArray& p1, const fArray& p2, fMat NN) noexcept {
	assert(all(p1.geq(0.0)) && all(p1.lt(1.0)));
	assert(all(p2.geq(0.0)) && all(p2.lt(1.0)));
	assert(p1.M()==DIM__ && p1.N()==1);
	assert(p2.M()==DIM__ && p2.N()==1);
	assert(B.square() && B.M()==DIM__);

	const auto n = mnorm(B.prod(cadd(NN,p2-p1)));
	const auto val = min(n);
	return {val,NN.getl(ones<fMat>(NN.M(),1),n.eq(val))};
}
d_b ll__::distb(const fMat& B, const fArray& p1, const fArray& p2, const double f, fMat NN) noexcept {
	assert(all(p1.geq(0.0)) && all(p1.lt(1.0)));
	assert(all(p2.geq(0.0)) && all(p2.lt(1.0)));
	assert(p1.M()==DIM__ && p1.N()==1);
	assert(p2.M()==DIM__ && p2.N()==1);
	assert(B.square() && B.M()==DIM__);

	auto n = mnorm(B.prod(cadd(NN,p2-p1)));
	if (n.L()==1) return {n[0],NN};
	
	// find acceptable permissive bond length lambda	
	auto sn = n;
	std::sort(sn.begin(),sn.end());
	auto fapbl = [&sn,f]() -> double {
		auto s=sn.begin(), m=sn.begin()+1, e=sn.end();
		if (f<0.0) { // factor mode
			for (; m<e; ++s, ++m)
				if (-(*s)*f<(*m)) return *s;
		} else { // direct cutoff mode
			for (; m<e; ++s, ++m)
				if (f<(*m)) return *s;
		}
		return *s;
	};

	const double val = fapbl();
	return {val,NN.getl(ones<fMat>(NN.M(),1),n.leq(val))};
}
fMat ll__::genDmat(const fMat& B, const fMat& P, const fMat& NN) noexcept {
	fMat res(P.N());
	std::for_each(res.dbegin(),res.dend(),[](auto& i){i=0.0;});
	for (auto i=P.ccBegin(),e=P.ccEnd(); i!=e; ++i)
		for (auto j=i+1; j!=e; ++j)
			res(i.i(),j.i()) = res(j.i(),i.i()) = dist(B,*j,*i,NN);
	return res;
}
fMat ll__::genDmat(const fMat& B, const fMat& P, const double f, const fMat& NN) noexcept {
	assert(f<-1.0 || f>=0.0);
	fMat res(P.N());
	std::for_each(res.dbegin(),res.dend(),[](auto& i){i=0.0;});
	for (auto i=P.ccBegin(),e=P.ccEnd(); i!=e; ++i)
		for (auto j=i+1; j!=e; ++j)
			res(i.i(),j.i()) = res(j.i(),i.i()) = dist(B,*j,*i,f,NN);
	return res;
}
fMat ll__::genDmat(const fMat& B, const fMat& P1, const fMat& P2, const fMat& NN) noexcept {
	fMat res(P1.N(),P2.N());
	auto ir = res.begin();
	for (auto j=P2.ccBegin(),je=P2.ccEnd(); j!=je; ++j)
		for (auto i=P1.ccBegin(),ie=P1.ccEnd(); i!=ie; ++i,++ir)
			*ir = dist(B,*i,*j,NN);
	return res;
}
fMat ll__::genDmat(const fMat& B, const fMat& P1, const fMat& P2, const double f, const fMat& NN) noexcept {
	assert(f<-1.0 || f>=0.0);
	fMat res(P1.N(),P2.N());
	auto ir = res.begin();
	for (auto j=P2.ccBegin(),je=P2.ccEnd(); j!=je; ++j)
		for (auto i=P1.ccBegin(),ie=P1.ccEnd(); i!=ie; ++i,++ir)
			*ir = dist(B,*i,*j,f,NN);
	return res;
}
fMat ll__::com(const fArray& inp, const fMat& w) noexcept {
	assert(all(inp.geq(0.0)) && all(inp.lt(1.0)));
	assert(w.empty() || w.L()==inp.N());
	if (inp.empty()) return inp;

	auto cos_ = fMat(inp.M(),inp.N());
	std::transform(inp.begin(),inp.end(),cos_.begin(),[](auto i){return std::cos(2.0*M_PI*i);});
	auto sin_ = fMat(inp.M(),inp.N());
	std::transform(inp.begin(),inp.end(),sin_.begin(),[](auto i){return std::sin(2.0*M_PI*i);});

	auto c_ = w.empty() ? nmean(cos_): nsum(rprod(cos_,w))/sum(w);
	auto s_ = w.empty() ? nmean(sin_): nsum(rprod(sin_,w))/sum(w);

	fMat res(inp.M(),1);
	for (auto ir=res.begin(),ci=c_.begin(),si=s_.begin(),er=res.end(); ir!=er; ++ir,++ci,++si)
		*ir = std::atan2(*si,*ci)/(2.0*M_PI);
	return res%=1.0;
}
wi ll__::dbscan(const fMat& D, const size_t minpts, const double eps) {

	// bitching
	if (!D.square())
		throw(std::invalid_argument("size(D) is ("+std::to_string(D.M())+","+
			std::to_string(D.N())+"), need D square"));
	if (any(D.lt(0.0))) throw(std::invalid_argument("need D>=0.0"));
	if (eps<0.0)
		throw(std::invalid_argument("eps=="+std::to_string(eps)+", need eps>=0.0"));
			

	// assigned and visited vectors, maximum cluster size and list of clusters
	std::vector<bool> a(D.M(),false);
	std::vector<bool> v(D.M(),false);
	size_t Cm = a.size();
	wi res; res.reserve(a.size());
	
	// expand cluster lambda
	std::function<void(wTv&)> ec;
	ec = [&](wTv& C) -> void {
		// update visited	
		v[C.back()] = true;

		// region query
		std::vector<size_t> rg; rg.reserve(D.M());
		for (size_t i=0; i<D.M(); ++i)
			if (D(i,C.back())<=eps) rg.push_back(i);
	
		// recursive call
		if (rg.size()<minpts) return;			// edge point
		for (auto i: rg) {				// cluster point
			if (v[i]) continue;
			C.push_back(i);
			ec(C);
		}
	};


	// find clusters
	for (size_t i=0; i!=a.size(); ++i) {
		if (a[i]) continue;

		// prepare cluster input and try to expand
		wTv C; C.reserve(Cm); C.push_back(i);
		std::fill(v.begin(),v.end(),false);
		ec(C);

		// extend result if expand successful
		if (C.size()>1) {
			C.shrink_to_fit();		// trim cluster C
			for (auto i: C) a[i] = true;	// update assigned vector
			Cm-=C.size();			// adapt maximum cluster size

			res.push_back(std::move(C));	// extend result
		}
	}

	// sort clusters in result
	for (auto& i: res)
		std::sort(i.begin(),i.end());

	// add noise points as clusters of size 1
	for (size_t i=0; i!=a.size(); ++i)
		if (!a[i]) res.push_back({i});

	// trim result and return
	res.shrink_to_fit();
	return res;	
}


// generate wannier hamiltonian
R_H<> ll__::genHam(const k_U& inp, const fMat& E, const fMat& R) noexcept {
	assert(!R.empty() && !E.empty() && !inp.empty());
	assert(R.M()==DIM__ && inp.dim()==DIM__);
	assert(E.N()==inp.size());

	// check E vs U
	assert(E.M()==inp.Nb());

	// generate wannier hamiltonian
	constexpr std::complex<double> mitpi(0.0,-2.0*M_PI);
	
	cMatv H; H.reserve(R.N());
	for (auto r=R.ccBegin(),re=R.ccEnd(); r!=re; ++r) {
		H.push_back(zeros<cMat>(inp.Nw()));

		for (auto k=inp.ccBegin(),ke=inp.ccEnd(); k!=ke; ++k)
			H.back() += lm__::T(inp[size_t(k)]).prod(diag(E.cAt(k.i()))).prod(inp[size_t(k)])
					* std::exp(mitpi*dot(*r,*k));
		H.back()/=inp.N();
	}

	return {R,std::move(H)};
}
template<class MT, class WT>
R_H<MT> ll__::genHam(const fMat& B, const fMat& Ap, const idv& id,
		const WT& W, fMat R, const bool strict) noexcept {
	
	assert(B.square());
	assert(B.M()==W.dim());
	assert(Ap.M()==B.M());
	assert(R.M()==W.dim());
	assert(id.size()==Ap.N());
	if (Ap.empty() || R.empty()) return {{},{}};

	// get type and block size
	const auto T = W.ind(id);
	const size_t Nw = std::accumulate(T.cbegin(),T.cend(),size_t(0),
			[&W](const size_t s, const size_t t)->size_t
			{ return s+W.Norb(t); });
	const double IR = W.radius();
	if (W.empty()) return {zeros<fMat>(W.dim(),1),{zeros<MT>(Nw,Nw)}};

	// generate H
	R_H<MT> res(std::move(R),std::vector<MT>(R.N(),zeros<MT>(Nw,Nw)));
	hctor(W,B,Ap,T,ll_writerMEM<MT>(res),IR,!strict);

	return res;
}
template R_H<fMat> ll__::genHam(const fMat& B, const fMat& Ap, const idv& id,
		const ll_hbonds& W, fMat R, const bool strict) noexcept;
template R_H<cMat> ll__::genHam(const fMat& B, const fMat& Ap, const idv& id,
		const ll_hbonds& W, fMat R, const bool strict) noexcept;
template R_H<fMat> ll__::genHam(const fMat& B, const fMat& Ap, const idv& id,
		const ll_hbondss& W, fMat R, const bool strict) noexcept;
template R_H<cMat> ll__::genHam(const fMat& B, const fMat& Ap, const idv& id,
		const ll_hbondss& W, fMat R, const bool strict) noexcept;
template<class MT, class WT>
R_H<MT> ll__::genHam(const fMat& B, const fMat& Ap, const idv& id,
		const rv& r, const WT& W, const bool strict) noexcept {

	assert(B.square());
	assert(B.M()==W.dim());
	assert(Ap.M()==B.M());
	assert(id.size()==Ap.N());
	if (Ap.empty()) return {{},{}};

	// get type and block size
	const auto T = W.ind(id);
	const size_t Nw = std::accumulate(T.cbegin(),T.cend(),size_t(0),
			[&W](const size_t s, const size_t t)->size_t
			{ return s+W.Norb(t); });
	const double IR = W.radius();
	if (W.empty()) return {zeros<fMat>(W.dim(),1),{zeros<MT>(Nw,Nw)}};

	// find appropriate R grid
	const auto R = getConnectedGrid(W,B,Ap,T,r,IR,!strict);
	
	// generate H
	R_H<MT> res(std::move(R),std::vector<MT>(R.N(),zeros<MT>(Nw,Nw)));
	hctor(W,B,Ap,T,ll_writerMEM<MT>(res),IR,!strict);

	return res;
}
template R_H<fMat> ll__::genHam(const fMat& B, const fMat& Ap, const idv& id,
		const rv& r, const ll_hbonds& W, const bool strict) noexcept;
template R_H<cMat> ll__::genHam(const fMat& B, const fMat& Ap, const idv& id,
		const rv& r, const ll_hbonds& W, const bool strict) noexcept;
template R_H<fMat> ll__::genHam(const fMat& B, const fMat& Ap, const idv& id,
		const rv& r, const ll_hbondss& W, const bool strict) noexcept;
template R_H<cMat> ll__::genHam(const fMat& B, const fMat& Ap, const idv& id,
		const rv& r, const ll_hbondss& W, const bool strict) noexcept;


// wannier matching
wi ll__::matchCenters(const fMat& B, const fMat& Ap, const fMat& Wp, const fMat& NN) noexcept {
	if (Ap.empty() || Wp.empty())
		return {{}};
	
	assert(B.square());
	assert(B.M()==Ap.M());
	assert(Ap.M()==Wp.M());
	assert(NN.M()==B.M());
	
	// lambda to find position in Ap with least distance from given p
	auto gmdi = [&B,&Ap,&NN](const fMat& p) -> size_t {
		double l=DBL_MAX; size_t mi=NPOS__;
		for (auto i=Ap.ccBegin(),e=Ap.ccEnd(); i!=e; ++i) {
			const double nl = dist(B,*i,p,NN);
			if (l>nl) { l=nl; mi=(i-Ap.ccBegin()); }
		}
		return mi;
	};

	// generate result
	wi res;
	res.resize(Ap.N());
	for (auto i=Wp.ccBegin(),e=Wp.ccEnd(); i!=e; ++i)
		res[gmdi(*i%1.0)].push_back(i.i());
	
	return res;
}
void ll__::balance(wi& I, const ll_cell& cell, const fMat& Wp, const fMat& NN) noexcept {
	assert(cell.dim()==Wp.M());
	assert(Wp.M()==NN.M());
	assert(all(Wp.gt(0.0)) && all(Wp.lt(1.0)));

	// stripped ids
	idv id; id.reserve(cell.N());
	for (const auto& s: cell.id())
		id.push_back(cell.stripId(s));

	// unique stripped ids
	idv uid = id;
	std::sort(uid.begin(),uid.end());
	uid.resize(std::distance(uid.begin(),std::unique(uid.begin(),uid.end())));

	// get indices for each id
	std::vector<std::vector<size_t>> Iid; Iid.reserve(uid.size());
	for (const auto& s: uid) {
		Iid.push_back({});
		for (size_t i=0; i!=id.size(); ++i)
			if (s == id[i]) Iid.back().push_back(i);
	}

	// balance matchings in I if possible
	for (const auto& i: Iid) {
		const size_t N = std::accumulate(i.cbegin(),i.cend(),size_t(0),
				[&I](const size_t s, const size_t j){return s+I[j].size();});
		if (N%i.size()) continue;

		std::vector<size_t> J; J.reserve(N);
		for (const auto t: i)
			for (const auto j: I[t]) J.push_back(j);
		auto D = genDmat(cell.B(),cell.Ap().get({},i),Wp.get({},J),NN);

		// find smallest element successively until all centers are equally distributed
		wi NI(i.size()); for (auto& j: NI) j.reserve(N/i.size());
		size_t cnt=0;
		while (cnt!=N) {
			const auto nm = nmin(D);
			const auto mm = mmin(nm.mat);

			const size_t m = mm.pos[0], n = nm.pos[mm.pos[0]];
			if (NI[m].size()<NI[m].capacity())
				NI[m].push_back(n), ++cnt;
			D(m,n)=DBL_MAX;
		}
		std::for_each(NI.begin(),NI.end(),[](auto& i)->void
			{std::reverse(i.begin(),i.end());});

		// reindex NI according to J
		for (auto& t: NI)
		for (auto& j: t) j = J[j];

		// replace entries in I with NI
		for (size_t j=0; j!=i.size(); ++j)
			I[i[j]] = NI[j];
	}
}
wi ll__::clusterize(const fMat& B, const fMat& Wp,
		const double eps, const size_t minpts, const fMat& NN) noexcept {
	if (Wp.empty()) return {{}};
	assert(B.square());
	assert(B.M()==Wp.M());
	assert(NN.M()==B.M());
	
	return dbscan(genDmat(B,Wp%1.0,NN),minpts,eps);
}
fMat ll__::genCenters(const fMat& B, const fMat& Wp, const fMat& s, const wi& I) noexcept {
	assert(B.square());
	assert(B.M()==Wp.M());
	assert(std::all_of(I.begin(),I.end(),[](const wTv& i){return !i.empty();}));
	assert(std::all_of(I.begin(),I.end(),
		[&Wp](const wTv& i){
		return std::all_of(i.begin(),i.end(),[&Wp](const size_t j){return j<Wp.N();});
		}));

	fMat res(B.M(),0); res.reserve(I.size());
	for (const auto& i: I)
		res.push_back(com(Wp.get({},i)%=1.0,s.get({},i)));	
	return res;
}


// wannier tools
aTv ll__::wiToT(const wi& I) noexcept {
	aTv res(std::accumulate(I.begin(),I.end(),size_t(0),
			[](const size_t s, const wTv& i){return s+i.size();}));

	for (size_t i=0; i!=res.size(); ++i)
		for (size_t j=0; j!=I.size(); ++j)
			if (std::find(I[j].begin(),I[j].end(),i)!=I[j].end()) {
				res[i]=j; break;
			}
	return res;
}
wi ll__::TToWi(const aTv& T) noexcept {
	wi res(*std::max_element(T.begin(),T.end())+1);
	for (size_t i=0; i!=res.size(); ++i)
		res[i].reserve(std::count(T.begin(),T.end(),i));

	for (size_t i=0; i!=T.size(); ++i)
		res[T[i]].push_back(i);

	return res;
}
bool ll__::checkWi(const wi& I) noexcept {
	return checkWi(I,std::accumulate(I.begin(),I.end(),size_t(0),
		[](const size_t a, const aTv& b){return a+b.size();}));
}
bool ll__::checkWi(const wi& I, const size_t N) noexcept {
	std::vector<bool> ck(N,false);
	for (const auto& i: I)
		for (const auto j: i)
			if (j>=N) return false;
			else {
				if (ck[j]) return false;
				ck[j] = true;
			}
	return std::all_of(ck.begin(),ck.end(),[](const bool i){return i;});
}
fMat ll__::checkHr(const R_H<>& hr) noexcept {
	fMat res(2,hr.size());
	auto i = res.begin();
	for (const auto& H: hr) {
		*i++ = std::abs(*std::max_element(c_reItr(H.data()),c_reItr(H.data()+H.L()),
				[](const double i, const double j){return std::abs(i)<std::abs(j);}));
		*i++ = std::abs(*std::max_element(c_imItr(H.data()),c_imItr(H.data()+H.L()),
				[](const double i, const double j){return std::abs(i)<std::abs(j);}));
	}
	return res;
}


vb_cb ll__::findBandEdges(const double Ef, const fMat& E) {
	vb_cb res = {-DBL_MAX,DBL_MAX};
	for (const auto i: E) {
		if (i<Ef) { if(i>res.vb) res.vb=i; }
		else { if(i<res.cb) res.cb=i; }
	}
	return res;
}
	

// simple bandstructure calculation
template<class MT=cMat>
fMat ll__::calcBS(const R_H<MT>& hr, const fMat& k, const size_t Nthreads) noexcept {
	assert(k.M()==hr.dim());
	assert(Nthreads||true);
	
	fMat res(hr.Nw(),k.N());
	
	constexpr std::complex<double> tpi(0.0,2.0*M_PI);
	
	const auto ie=k.ccEnd();
	#pragma omp parallel for num_threads(Nthreads)
	for (auto i=k.ccBegin(); i<ie; ++i) {
		auto h = zeros<cMat>(hr.Nw());
		for (auto j=hr.ccBegin(),je=hr.ccEnd(); j!=je; ++j)	
			h += hr[size_t(j)]*std::exp(tpi*dot(*i,*j));
		res.cAt((size_t)i) = eigh(h);
	}
	return res;
}
template fMat ll__::calcBS(const R_H<fMat>& hr, const fMat& k, const size_t Nthreads) noexcept;
template fMat ll__::calcBS(const R_H<cMat>& hr, const fMat& k, const size_t Nthreads) noexcept;

template<class MT=cMat>
fMat ll__::calcFoldedBS(const R_H<MT>& hr, const fMat& k, const fMat& B,
		const fMat& Bp, const size_t Nthreads) noexcept {
	assert(!k.empty());
	assert(B.square() && B.M()==k.M());
	assert(size(B)==size(Bp));
	assert(Bp.leftDivide(B)==round(Bp.leftDivide(B)));
	assert(std::abs(det(B))>mtol());
	assert(std::abs(det(Bp))>mtol());

	// find k points and 'type' in primitive brillouin zone
	auto bzp = ll_cell(inv(B).T(),k%1.0,aCv(k.N(),1)).changeBasis(inv(Bp).T());

	// calculate energies for all k points
	auto Ep = calcBS<MT>(hr,bzp.Ap(),Nthreads);

	// combine energies into new bands
	size_t Nwf = std::round(std::abs(det(B)/det(Bp))) * hr.Nw();
	fMat res(Nwf,0); res.reserve(k.N());
	
	for (const auto t: bzp.types()) {
		auto e = Ep.get({},bzp.ind(t)).C();
		std::sort(e.begin(),e.end());
		res.push_back(e);
	}

	// return
	return res;
}
template fMat ll__::calcFoldedBS(const R_H<fMat>& hr, const fMat& k, const fMat& B,
		const fMat& Bp, const size_t Nthreads) noexcept;
template fMat ll__::calcFoldedBS(const R_H<cMat>& hr, const fMat& k, const fMat& B,
		const fMat& Bp, const size_t Nthreads) noexcept;


// bandstructure calculation including gradients and curvature tensors
template<class MT=cMat>
ll__::egc<fMat> ll__::calcBS_gc(const R_H<MT>& hr, const fMat& k,
		const fMat& B, const size_t Nthreads) noexcept {
	assert(k.M()==hr.dim());
	assert(Nthreads||true);

	egc<fMat> res = {
		fMat(hr.Nw(),k.N()),
		ll_mesh<>(hr.Nw(),{0,1},{k.M(),k.N()}),
		ll_mesh<>(hr.Nw(),{0,1,2},{k.M(),k.M(),k.N()})
	};

	constexpr std::complex<double> tpi(0.0,2.0*M_PI);
	constexpr std::complex<double> i_(0.0,1.0);
	const auto BR = B.prod(hr.R());
	
	constexpr double thr = 1e-6;

	const auto ie=k.ccEnd();
	#pragma omp parallel for num_threads(Nthreads)
	for (auto i=k.ccBegin(); i<ie; ++i) {

		// energies, degenerate subspaces, states U(k)
		std::vector<std::vector<size_t>> DSI;
		cMat U, UT;
		{
			auto h = zeros<cMat>(hr.front().M());
			auto jh = hr.cbegin();
			for (auto jr=hr.ccBegin(),je=hr.ccEnd(); jr!=je; ++jr,++jh)	
				h += (*jh)*std::exp(tpi*dot(*i,*jr));
			
			// diagonalize
			auto ev = eighv(h);
			const auto J = aux::sorted_order(ev.E.cbegin(),ev.E.cend());
			aux::reorder(ev.E.begin(),J); aux::reorder(ev.V.cBegin(),J);
			
			// find degenerate subspaces
			size_t cnt=0;
			for (auto j=ev.E.cbegin()+1, je=ev.E.cend(); j<je; ++j) {
				DSI.push_back({cnt++});
				for (;j!=je && std::abs(*j-*(j-1))<thr; ++j)
					DSI.back().push_back(cnt++);
				if (DSI.back().size()<2)
					DSI.pop_back();
			}
			
			res.E.cAt((size_t)i) = ev.E; U = std::move(ev.V); UT = T(U);
		}
		
		// dH(k), gradient
		auto jG = res.G.begin(0,{0,(size_t)i});
		std::vector<cMat> dh(k.M(),zeros<cMat>(hr.Nw()));
		for (size_t m=0; m!=k.M(); ++m,++jG) {
			auto jh = hr.cbegin();
			auto brm = BR.rAt(m).cbegin();
			for (auto jr=hr.ccBegin(),je=hr.ccEnd(); jr!=je; ++jr,++jh,++brm)
				dh[m] += i_*(*jh)*(*brm)*std::exp(tpi*dot(*i,*jr));
			
			dh[m] = UT.prod(dh[m]).prod(U);
			*jG = diag(dh[m]);
		}

		// curvatures
		for (size_t n=0; n!=k.M(); ++n) {

			// D matrix
			auto D = zeros<cMat>(hr.Nw());
			for (size_t n_=0; n_!=hr.Nw(); ++n_)
			for (size_t m_=0; m_!=hr.Nw(); ++m_)
				D(m_,n_) = (m_==n_) ? 0.0: dh[n](m_,n_)/
					(res.E(n_,(size_t)i)-res.E(m_,(size_t)i));
			for (const auto& I: DSI)
				for (const auto i: I)
				for (const auto j: I)
					D(i,j) = 0.0;
			
			auto jC = res.C.begin(0,{0,n,(size_t)i});
			for (size_t m=0; m!=k.M(); ++m,++jC) {
				auto ddh = zeros<cMat>(hr.Nw());
				
				auto jh = hr.cbegin();
				auto brm = BR.rAt(m).cbegin(), brn = BR.rAt(n).cbegin();
				for (auto jr=hr.ccBegin(),je=hr.ccEnd(); jr!=je; ++jr,++jh,++brn,++brm)
					ddh -= (*jh)*(*brn)*(*brm)*std::exp(tpi*dot(*i,*jr));
				*jC = diag(UT.prod(ddh).prod(U)) + 2.0*diag(dh[m].prod(D));
			}
		}
	}

	return res;
}
template ll__::egc<fMat> ll__::calcBS_gc(const R_H<fMat>& hr, const fMat& k,
		const fMat& B, const size_t Nthreads) noexcept;
template ll__::egc<fMat> ll__::calcBS_gc(const R_H<cMat>& hr, const fMat& k,
		const fMat& B, const size_t Nthreads) noexcept;

template<class MT=cMat>
egc<fMat> ll__::calcFoldedBS_gc(const R_H<MT>& hr, const fMat& k,
		const fMat& B, const fMat& Bp, const size_t Nthreads) noexcept {
	assert(!k.empty());
	assert(B.square() && B.M()==k.M());
	assert(size(B)==size(Bp));
	assert(Bp.leftDivide(B)==round(Bp.leftDivide(B)));
	assert(std::abs(det(B))>mtol());
	assert(std::abs(det(Bp))>mtol());


	// find k points and 'type' in primitive brillouin zone
	auto bzp = ll_cell(inv(B).T(),k%1.0,aCv(k.N(),1)).changeBasis(inv(Bp).T());

	// calculate energies for all k points
	auto EGCp = calcBS_gc<MT>(hr,bzp.Ap(),Bp,Nthreads);

	// create result container
	const size_t Nwf = std::round(std::abs(det(B)/det(Bp))) * hr.Nw();
	egc<fMat> res = {
		fMat(Nwf,k.N()),
		ll_mesh<>(Nwf,{0,1},{k.M(),k.N()}),
		ll_mesh<>(Nwf,{0,1,2},{k.M(),k.M(),k.N()})
	};
	
	// sort and combine energies into new bands
	for (size_t ki=0; ki!=k.N(); ++ki) {
		const auto KJ = bzp.ind(ki);
		
		// get sorted order for energies
		auto e = EGCp.E.get({},KJ).C();
		const auto J = aux::sorted_order(e.cbegin(),e.cend());
		
		// get energy
		aux::reorder(e.begin(),J); res.E.cAt(ki) = e;

		// get gradient
		for (size_t m=0; m!=k.M(); ++m) {
			fMat buff(EGCp.Nb(),0); buff.reserve(KJ.size());
			for (const auto kj: KJ)
				buff.push_back(EGCp.G.cAt({m,kj}));
			aux::reorder(buff.begin(),J);
			res.G.cAt({m,ki}) = buff.C();
		}

		// get curvature tensor
		for (size_t m=0; m!=k.M(); ++m)
		for (size_t n=0; n!=k.M(); ++n) {
			fMat buff(EGCp.Nb(),0); buff.reserve(KJ.size());
			for (const auto kj: KJ)
				buff.push_back(EGCp.C.cAt({m,n,kj}));
			aux::reorder(buff.begin(),J);
			res.C.cAt({m,n,ki}) = buff.C();
		}
	}

	return res;
}
template egc<fMat> ll__::calcFoldedBS_gc(const R_H<fMat>& hr, const fMat& k,
		const fMat& B, const fMat& Bp, const size_t Nthreads) noexcept;
template egc<fMat> ll__::calcFoldedBS_gc(const R_H<cMat>& hr, const fMat& k,
		const fMat& B, const fMat& Bp, const size_t Nthreads) noexcept;
	

// wbh related
template<class WT>
ll__::struct_report ll__::analyzeStructure(const fMat& Ap, const idv& id,
					   const WT& W, const double tol, const double f) noexcept {
	assert(Ap.N()==id.size());
	assert(Ap.M()==W.dim());

	// types
	typedef struct_report::iiN  iiN;
	typedef struct_report::iijN iijN;

	// type indices for id
	const auto I = W.ind(id);
	
	// prepare result
	struct_report res;
	res.exact.reserve(id.size()*id.size());
	res.approx.reserve(id.size()*id.size());
	res.unmatched = I; std::sort(res.unmatched.begin(),res.unmatched.end());

	// extend result lambda
	const auto extVec = [](auto& vec, const auto& inp) -> void {
		const auto itr = std::lower_bound(
			vec.begin(),vec.end(),inp,
			[](const auto& i, const auto& j) -> bool {
				if (i.i1<j.i1) return true;
				if (i.i1>j.i1) return false;
				return i.i2<j.i2;
			});
		if (itr == vec.end() || itr->i1!=inp.i1 || itr->i2!=inp.i2)
			vec.insert(itr,inp);
		else ++(itr->N);
	};
	const auto uptUnm = [&res](const auto& inp) -> void {
		for (auto i=inp.cbegin(),e=inp.cend(); i!=e; ++i) {
			const auto itr = std::lower_bound(
				res.unmatched.cbegin(),res.unmatched.cend(),*i);
			if (itr!=res.unmatched.cend() && *itr==*i)
				res.unmatched.erase(itr);
		}
	};

	// cutoff radius
	const double Rcut = W.radius()*f;

	// find matches
	for (auto p1=Ap.ccBegin(),e=Ap.ccEnd(); p1!=e; ++p1)
	for (auto p2=p1; p2!=e; ++p2) {
		
		const fMat b = *p2-*p1;
		if (norm(*p2-*p1)>Rcut) continue;
		
		const size_t i1 = I[(size_t)p1];
		const size_t i2 = I[(size_t)p2];
		
		// check direct match
		const cMat& Hdirect = W.getInteraction({i1,i2},b);
		if (Hdirect!=W.eH() && any(abs(Hdirect).gt(tol))) {
			extVec(res.exact,iiN{{i1,i2,0}});
			extVec(res.exact,iiN{{i2,i1,0}});
			uptUnm(std::array<size_t,2>{{i1,i2}});
			continue;
		}

		// check approximate match
		const auto Happrox = W.getApproximateInteraction({i1,i2},b);
		if (Happrox.H!=W.eH() && any(abs(Happrox.H).gt(tol))) {
			if (Happrox.pi!=NPOS__)
				extVec(res.approx,iijN{{i1,i2,Happrox.pi}});
			if (Happrox.mi!=NPOS__)
				extVec(res.approx,iijN{{i2,i1,Happrox.mi}});
			uptUnm(std::array<size_t,3>{{i1,Happrox.pi,Happrox.mi}});
			continue;
		}
	}

	return res;
}
template ll__::struct_report ll__::analyzeStructure(const fMat& Ap, const idv& id,
		const ll_hbonds& W, const double tol, const double f) noexcept;
template ll__::struct_report ll__::analyzeStructure(const fMat& Ap, const idv& id,
		const ll_hbondss& W, const double tol, const double f) noexcept;
