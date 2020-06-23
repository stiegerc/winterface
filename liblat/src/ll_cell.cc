// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "ll_cell.h"
#include "aux_sort.h"
#include "aux_io.h"
#include "ll_types.h"
#include "ll_lambda.h"
#include "ll_fn.h"
#include "ll_io.h"
#include "ll_bonds.h"
#include <functional>
#include <algorithm>
#include <cmath>
#include <regex>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <locale>
#include <ctime>
#include <cfloat>
#include <stdexcept>
#include <utility>
#include <numeric>
#include <cassert>


using namespace lm__;
using namespace ll__;
using namespace aux;

// empty string
const std::string ll_cell::e_("");

// contructors
ll_cell::ll_cell(fMat B, fMat Ap, aCv N, idv id) noexcept:
		mat_cb<fMat,fArray>(std::move(Ap)), B_(std::move(B)), id_(std::move(id)) {

	// check basis B
	assert(!B_.empty());
	assert(B_.square());
	assert(ops::nz(volnorm()));
	
	// check atomic positions
	assert(!mat_.empty());
	assert(B_.M()==mat_.M());
	assert(all(mat_.geq(0.0)&mat_.lt(1.0)));
	
	// check atomic count vector and ids
	assert(std::accumulate(N.cbegin(),N.cend(),size_t(0))==mat_.N());
	assert(id_.empty() || id_.size()==N.size());


	// force positions into positive sector
	mat_ %= 1.0;

	// construct T_
	T_.reserve(N.size()+1); T_.push_back(0);
	for (const auto i: N) T_.push_back(T_.back()+i);
	
	// sort Ap inside types
	sortAp_(types());

	// index ids
	indexId_();
}
ll_cell::ll_cell(fMat B, fMat Ap, idv id) noexcept {
	assert(Ap.N()==id.size());

	// sort ids and Ap
	const auto I = aux::sorted_order(id.cbegin(),id.cend());
	aux::reorder(id.begin(),I); aux::reorder(Ap.cBegin(),I);

	// get unique ids and N vector
	auto uid = id; uid.resize(std::distance(uid.begin(),std::unique(uid.begin(),uid.end())));
	aCv N(uid.size()); std::transform(uid.cbegin(),uid.cend(),N.begin(),
				[&id](const std::string& s)->size_t
				{ return std::distance(std::lower_bound(id.cbegin(),id.cend(),s),
						       std::upper_bound(id.cbegin(),id.cend(),s)); });
	
	// create cell
	*this = ll_cell(std::move(B),std::move(Ap),N,std::move(uid));
}
ll_cell::ll_cell(const std::string& fileName, const size_t d) {

	switch(detectFileType(fileName)) {
	case 0: // POSCAR file
	{
		auto P = readPOSCAR(fileName,d);
		*this = P.id.empty() ? ll_cell(std::move(P.B),std::move(P.Ap),P.N):
				       ll_cell(std::move(P.B),std::move(P.Ap),P.N,std::move(P.id));
	}
	break;
	case 1: // wannier90 wout file
	{
		const auto B = readB(fileName);
		const auto p = readAp(fileName,false);
		const auto Ap = B.leftDivide(p.Ap())%1.0;
		*this = ll_cell(std::move(B),std::move(Ap),std::move(p.id()));
	}
	break;
	case 2: // OMEN lattice file
	{
		const auto olf = readOlf(fileName);
		*this = ll_cell(std::move(olf.B),olf.B.leftDivide(olf.Ap)%1.0,
				std::move(olf.id));
	}
	break;
	}
}


// information
lm__::c_fColItr ll_cell::find(const fArray& inp, const aT t) const noexcept {
	assert(inp.M()==dim() && inp.N()==1);
	return std::find_if(ccBegin(t),ccEnd(t),[&inp](const fCol& c)->bool{
			return std::equal(c.cbegin(),c.cend(),inp.cbegin(),
				[](const double a, const double b)->bool{
					// distance for numbers on integer grid < mtol
					const double tmp = std::abs(a-b);
					return std::abs(tmp-std::round(tmp))<mtol();
				});
		});
}
bool ll_cell::lVec(const fCol& inp) const noexcept {
	assert(inp.L()==dim());
	
	// check for bad vector lambda
	const auto bv = [this,&inp](const aT t) -> bool {
		for (auto i=ccBegin(t),e=ccEnd(t),e_=ccEnd(); i!=e; ++i)
			if (find(*i+inp,t)==e_) return true;
		return false;
	};

	for (const auto t: types())
		if (bv(t)) return false;
	return true;
}
bool ll_cell::lVec(const fMat& inp) const noexcept {
	if (inp.empty() || inp.M()!=dim())
		return false;

	for (auto i=inp.ccBegin(),e=inp.ccEnd(); i!=e; ++i)
		if (!lVec(*i)) return false;
	return true;
}
bool ll_cell::primitive() const noexcept {
	if (empty()) return true;

	const auto lft = leastFreqType();
	for (auto i=ccBegin(lft)+1, e=ccEnd(lft); i!=e; ++i)
		for (auto j=ccBegin(lft); j!=i; ++j) {
			const auto vec = *j-*i;
			if (lVec(vec)) return false;
		}

	return true;
}
bool ll_cell::validBasis(fMat rhs) const noexcept {
	if (rhs.M()!=dim() || !rhs.square()) return false;
	if (ll__::volnorm(rhs)<mtol()) return false;

	B().leftDivideEq(rhs);
	for (auto i=rhs.ccBegin(),e=rhs.ccEnd(); i!=e; ++i)
		if (!lVec(*i)) return false;
	return true;
}
ll__::rv ll_cell::r(const double f) const noexcept {
	if (empty()) return rv(dim(),false);

	// nn bond length
	const double bl = getBonds(-1.1,ll__::NN).radius();

	// find restricted dimensions by comparing the largest gap
	// in each dimension to the bond length using factor f
	rv res; res.reserve(dim());
	for (size_t d=0; d!=dim(); ++d) {
		auto tmp = Ap().rGet(d); std::sort(tmp.begin(),tmp.end());

		// find largest gap in tmp
		const double h = tmp.front()+1.0;
		for (auto j=tmp.begin()+1, e=tmp.end(); j<e; ++j)
			*(j-1) = *j - *(j-1);
		tmp.back() = h-tmp.back();
		const double lg = *std::max_element(tmp.cbegin(),tmp.cend());

		res.push_back(lg*norm(B().cAt(d))>f*bl);
	}

	return res;
}
double ll_cell::directTol(const double tol) const noexcept {
	if (empty()) return 0.0;
	
	const fMat IB = lm__::inv(B())*tol;
	std::vector<double> tmp(dim());
	std::transform(IB.crBegin(),IB.crEnd(),tmp.begin(),
		[](const auto& r)->double{return lm__::dot(r,r);});
	return std::sqrt(*std::max_element(tmp.cbegin(),tmp.cend()));
}


// valid types and indices
std::vector<aTv> ll_cell::equalTypes() const noexcept {
	
	std::vector<aTv> res;
	
	if (id().empty())
		for (const auto t: types())
			res.push_back({t});
	else {
		auto T = types();
		for (auto i=T.cbegin(); i!=T.cend();) {
			res.push_back({});
			const auto sid = stripId(id(*i));
			for (auto j=i; j!=T.cend();)
				if (sid==stripId(id(*j)))
					res.back().push_back(*j), T.erase(j);
				else ++j;
		}
	}
	
	return res;
}


// id information	
idv ll_cell::fundamentalIds() const noexcept {
	if (id().empty()) return {};

	const std::regex rgx("[(:_]|\\)-[0-9]+|\\)");

	// find all contenders
	idv res;
	for (const auto& s: id())
		for (std::sregex_token_iterator i(s.cbegin(),s.cend(),rgx,-1),e; i!=e; ++i)
			if (!i->str().empty()) res.push_back(*i);
	
	// sort and remove duplicates
	std::sort(res.begin(),res.end());
	res.resize(std::distance(res.begin(),std::unique(res.begin(),res.end())));
	res.shrink_to_fit();

	return res;
}
idv ll_cell::compositeIds() const {
	if (id().empty()) return {};

	// find all contenders
	std::vector<std::string> res;
	try {
		for (const auto& s: id()) {
			// look for parenthesis, ensure first is not ')'
			const size_t i = s.find_first_of("()");
			if (i==std::string::npos) continue;
			if (s[i]==')') throw(std::invalid_argument(s));

			// find all (nestes) ids
			for (auto i=s.cbegin(),e=s.cend(); i!=e; ++i)
				if (*i=='(') {
					auto j=i; int brcnt=0;
					do {
						if (j==e || brcnt<0) throw(std::invalid_argument(s));
						switch (*j++) {
							case '(': ++brcnt; break;
							case ')': --brcnt; break;
						}
					} while(brcnt);
					while (j!=e && (*j=='-' || std::isdigit(*j)))
						++j;
					res.push_back(std::string(i,j));
				}
		}
	} catch(const std::invalid_argument& e) {
		throw(std::invalid_argument("id '"+std::string(e.what())+"': bad format"));
	}
	
	// sort and remove duplicates
	std::sort(res.begin(),res.end());
	res.resize(std::distance(res.begin(),std::unique(res.begin(),res.end())));
	res.shrink_to_fit();

	return res;
}

// modification
ll_cell& ll_cell::swapDim(const size_t n1, const size_t n2) noexcept {
	assert(n1<dim());
	assert(n2<dim());
	
	swap(B_.cAt(n1),B_.cAt(n2));
	swap(mat_.rAt(n1),mat_.rAt(n2));
	sortAp_(types());
	return *this;
}
ll_cell& ll_cell::invDim(const size_t n) noexcept {
	assert(n<dim());
	
	B_.cAt(n) = -B_.cAt(n);
	mat_.rAt(n) = -mat_.rAt(n)%1.0;
	sortAp_(types());
	return *this;
}
ll_cell& ll_cell::orient(const double sgn, const rv& r) noexcept {
	if (std::signbit(sign())==std::signbit(sgn)) return *this;
	
	const auto nri = ninds(r);
	if (!nri.size()) return *this;
	if (nri.size()==1) return invDim(nri.front());
	return swapDim(nri.front(),nri.back());
}
ll_cell& ll_cell::permute(const fMat& P) noexcept {
	assert(size(B())==size(P));
	assert(P.permutation());
	
	B_ = B_.prod(P);
	mat_ = P.prod(mat_);
	sortAp_(types());
	return *this;
}
ll_cell& ll_cell::rotate(const fMat& R) noexcept {
	assert(size(B())==size(R));
	assert(R.onb());

	B_ = R.prod(B_);
	return *this;
}
ll_cell& ll_cell::scale(const double f) noexcept {
	assert(f);

	B_ *= f;
	return *this;
}
ll_cell& ll_cell::scale(const fArray& f) noexcept {
	assert(f.size()==dim());
	assert(std::none_of(f.cbegin(),f.cend(),
		[](const auto i)->bool{return !i;}));

	auto j=f.cbegin();
	for (auto i=B_.cBegin(),e=B_.cEnd(); i!=e; ++i,++j)
		*i *= *j;
	return *this;
}




ll_cell& ll_cell::changeBasis(const fMat& NB) noexcept {
	assert(validBasis(NB));

	// transformation matrix and limits
	const auto TRm = NB.leftDivide(B());
	const auto ulim = 2.0*round(nmax(ceil(abs(B().leftDivide(NB)))).mat);
	const auto llim = -ulim;

	// new containers for Ap and T
	fMat nAp(dim(),std::round(N()/ll__::vol(TRm)));
	std::vector<size_t> nT(T_.size());
	std::transform(T_.cbegin(),T_.cend(),nT.begin(),
		[N=N(),nN=nAp.N()](const aT t){return (t*nN)/N;});

	// counter
	size_t cnt=0; aCv tcnt(Nspecies(),0);

	// lambda to insert position into nAp
	auto ins_ = [&nAp,&nT,&tcnt,&cnt](const fMat pos, const aT t) -> bool {
		const auto b = nAp.cBegin()+nT[t];
		const auto e = b+tcnt[t];
		const auto i = std::lower_bound(b,e,pos);

		// return false if pos is already included
		if (i!=e && *i==pos) return false;

		// found new position, shift positions >=i
		memmove(i->data()+i->M(),i->data(),(e-i)*i->M()*sizeof(double));

		// insert new pos at i,increase counter and return
		*i=pos; ++tcnt[t]; ++cnt;
		return true;
	};

	// lambda to insert for all using riv
	auto ins__ = [&ins_,&tcnt,&nT,&TRm,this](const fMat riv) -> void {
		for (const auto t: types()) {
			const aC N = nT[t+1]-nT[t];
			if (tcnt[t]==N) continue;
			for (auto j=ccBegin(t),e=ccEnd(t); j!=e; ++j) {
				ins_(TRm.prod(*j+riv)%1.0,t);
				if (tcnt[t]==N) break;
			}
		}
	};
	
	// find new positions using random image vectors
	for (size_t cc=0; cc!=100000; ++cc) {
		for (const auto t: types())
		for (auto j=ccBegin(t),e=ccEnd(t); j!=e; ++j) {
			const auto riv = randi<fMat>(llim,ulim);
			if (ins_(TRm.prod(*j+riv)%1.0,t))
				ins__(riv);

			// all positions found, replace and return
			if (cnt==nAp.N()) {
				B_ = std::move(NB); mat_ = std::move(nAp); T_ = std::move(nT);
				return *this;
			}
		}
	}
	
	assert(false);
	return *this;
}
ll_cell& ll_cell::makePrimitive_(fMat C, const rv& r, const std::function<bool(const fMat&)>& p) noexcept {
	assert(r.size()==dim());
	if (empty() || N()==1 || std::all_of(r.begin(),r.end(),[](const bool i){return i;})) return *this;

	// recursive find smaller basis lambda
	std::function<size_t(fMat&, const size_t)> recur;
	recur = [this,&r,&p,&recur](fMat& C, const size_t id) -> size_t {
		// advance pos lambda
		auto adv = [this,&r](size_t& pos) -> void {
			while (r[pos] && pos!=dim()) {
				++pos;
			}
		};
		
		// column position in C, rank of C and success bool
		size_t R=rank(C), pos=0; adv(pos); bool found=false;

		// extend base vectors lambda
		const auto ebv = [&C,&pos,&adv,&R,&found,this](const fMat& vec) -> size_t {
			assert(pos<dim());
			C.cAt(pos) = vec;
			if (rank(C)!=R) ++R, adv(++pos), found=true;
			return R;
		};
	
		// find internal lattice vector lambda
		auto filv = [this,&p,&ebv,&found](const aT t) -> bool {
			assert(validType(t));
			for (auto i=ccBegin(t),e=ccEnd(t); i!=e; ++i)
				for (auto j=i+1; j!=e; ++j) {
					const auto vec = *i-*j;
					if (!p(vec) || !lVec(vec.cFront())) continue;
					if (ebv(vec)==dim()) return found;
				}
			return found;
		};
	
		// find new basis, return if no internal vector is found
		if (!filv(leastFreqType())) return id;

		// extend by unit vectors until complete or return if not possible
		for (size_t i=0; R!=dim() && i!=dim(); ++i)
			if (!r[i] && ebv(cId<fMat>(dim(),i))==dim()) break;
		if (R!=dim()) return id;

		// change basis
		changeBasis(B().prod(C));
		
		// create new C
		for (size_t i=0; i!=dim(); ++i)
			C.cAt(i) = r[i] ? cId<fMat>(dim(),i): zeros<fMat>(dim(),1);

		return recur(C,id+1);
	};

	// recursively find smaller basis
	const auto sgn = sign();
	recur(C,0);
	return orient(sgn,r);
}
ll_cell& ll_cell::makePrimitive(const rv& r) noexcept {
	fMat C(dim(),0); C.reserve(dim());
	for (size_t i=0; i!=r.size(); ++i)
		C.push_back(r[i] ? lm__::cId<fMat>(dim(),i): lm__::zeros<fMat>(dim(),1));
	
	return makePrimitive_(C,r,[](const fMat& vec){return true;});
}
ll_cell& ll_cell::makePrimitiveInSubspace(const rv& r) noexcept {
	assert(r.size()==dim());
#ifndef NDEBUG
	{ const auto S = B().get(ll__::inds(r),{});
	assert(T(S).prod(S)==diag(mnormsq(S))); }
#endif

	fMat C(dim(),0); C.reserve(dim());
	for (size_t i=0; i!=r.size(); ++i)
		C.push_back(r[i] ? lm__::cId<fMat>(dim(),i): lm__::zeros<fMat>(dim(),1));

	auto inSubspace = [ri=ll__::inds(r)](const auto& vec) -> bool {
		const auto i=vec.cbegin();
		for (const auto j: ri)
			if (ops::nz(*(i+j))) return false;
		return true;
	};

	return makePrimitive_(C,r,inSubspace);
}
ll_cell& ll_cell::shift(const fCol& sh) noexcept {
	assert(sh.L()==dim());
	cadd(mat_,sh)%=1.0; sortAp_(types());
	return *this;
}
ll_cell& ll_cell::autoShift(const rv& r) noexcept {
	assert(r.size()==dim());

	const auto nri = ll__::ninds(r);
	const auto ri =  ll__::inds(r);
	if (!nri.size()) return *this;

	// get positions with minimal distance to origin in restricted dimensions
	fMat pos;
	if (ri.size()) {
		fMat rAp = Ap();
		for (const auto i: nri)
			rAp.rAt(i) = .0;
		const auto n = mnorm(B().prod(rAp));
		pos = Ap().get({},lm__::find(n.eq(min(n))));
	} else  pos = Ap();

	// find shift resulting in minimal volume encapsulating frame
	fMat sh = -pos;
	for (const auto i: ri)
		sh.rAt(i) = .0;
	auto j = sh.ccBegin(); double V = DBL_MAX;
	for (auto i=sh.ccBegin(),ie=sh.ccEnd(); i!=ie; ++i) {
		auto spos = pos; cadd(spos,*i); spos %= 1.0;
		
		double cV = 0.0;
		for (auto s=spos.ccBegin(), se=spos.ccEnd(); s!=se; ++s) {
			const double ccV = prod(*s);
			if (ccV>cV) cV=ccV;
		}

		if (cV<V) V=cV, j=i;
	}

	return shift(*j);
}
ll_cell& ll_cell::diversify(const aTv& ts) noexcept {
	assert(std::all_of(ts.cbegin(),ts.cend(),[this](const size_t t){return validType(t);}));
	assert(std::is_sorted(ts.cbegin(),ts.cend()));
	if (ts.empty()) return *this;

	// get number of elements for ts
	size_t N=0;
	for (auto t: ts) N += Ntype(t);
	
	// get new id_
	if (!id().empty()) {
		idv nid; nid.reserve(id_.size()+N-ts.size());
		
		auto ct = ts.cbegin();
		for (auto t: types())
			if (t == *ct) {
				for (size_t j=T_[t]+1; j<=T_[t+1]; ++j)
					nid.push_back(id(t));
				++ct;
			} else nid.push_back(id(t));

		id_ = std::move(nid);
	}

	// get new T_
	{
		std::vector<size_t> nT_; nT_.reserve(T_.size()+N-ts.size());
		
		nT_.push_back(0);
		auto ct = ts.cbegin();
		for (auto t: types())
			if (t == *ct) {
				for (size_t j=T_[t]+1; j<=T_[t+1]; ++j)
					nT_.push_back(j);
				++ct;
			} else nT_.push_back(T_[t+1]);

		T_ = std::move(nT_);
	}

	indexId_();

	return *this;
}
ll_cell& ll_cell::collectivize(const aTv& ts) noexcept {
	assert(std::all_of(ts.cbegin(),ts.cend(),[this](const size_t t){return validType(t);}));
	assert(std::is_sorted(ts.cbegin(),ts.cend()));
	if (ts.size()<2) return *this;

	// get new id_
	if (!id().empty()) {
		// remove from the back of ts except the first
		for (auto i=ts.cend()-1,e=ts.cbegin(); i!=e; --i)
			id_.erase(id_.cbegin()+*i);
	}
	
	// get new Ap and T_
	{
		std::vector<size_t> nT_; nT_.reserve(T_.size()-ts.size());
		fMat nAp(dim(),0);  nAp.reserve(N());

		nT_.push_back(0);
		auto ct = ts.cbegin();
		for (const auto t: types()) {
			if (t == ts.front()) {
				size_t N=0;
				for (const auto t_: ts) {
					N+=Ntype(t_);
					std::for_each(ccBegin(t_),ccEnd(t_),
						[&nAp](const auto& c)->void{nAp.push_back(c);});
				}
				nT_.push_back(nT_.back()+N);
			}
			if (ct!=ts.cend() && *ct==t) ++ct;
			else {
				std::for_each(ccBegin(t),ccEnd(t),
					[&nAp](const auto& c)->void{nAp.push_back(c);});
				nT_.push_back(nT_.back()+Ntype(t));
			}
		
		}
			
		T_ = std::move(nT_);
		mat_ = std::move(nAp);
	}

	// sort Ap in collectivized types and reindex id
	std::sort(cBegin(ts.front()),cEnd(ts.front()));
	indexId_();

	return *this;
}
ll_cell& ll_cell::merge(const ll_cell& inp) noexcept {
	assert(B()==inp.B());
	if (inp.empty()) return *this;
	if (this->empty()) return *this=inp;
	assert(!id().empty() && !inp.id().empty());

	// create union of both ids
	idv nid(this->Nspecies()+inp.Nspecies());
	nid.resize(std::distance(nid.begin(),
		std::set_union(this->id().cbegin(),this->id().cend(),
				 inp.id().cbegin(),  inp.id().cend(),
			         nid.begin())));
	nid.shrink_to_fit();

	// create new Ap and T_ vector omitting duplicates
	std::vector<size_t> nT; nT.reserve(nid.size()+1); nT.push_back(0);
	fMat nAp(dim(),N()+inp.N());
	
	auto itr = nAp.cBegin();
	for (const auto& s: nid) {
		const aT t1 = this->type(s), t2 = inp.type(s);
		itr = (this->validType(t1) && this->validType(t2)) ?
			std::set_union(this->ccBegin(t1),this->ccEnd(t1),
				         inp.ccBegin(t2),  inp.ccEnd(t2),itr):
		      (this->validType(t1)) ?
			std::copy(     this->ccBegin(t1),this->ccEnd(t1),itr):
		      (  inp.validType(t2)) ?
			std::copy(       inp.ccBegin(t2),  inp.ccEnd(t2),itr):
			itr;
		nT.push_back((size_t)itr);
	}

	// shrink Ap to fit and replace containers
	nAp.resize((size_t)itr); nAp.shrink_to_fit();
	this->mat_ = std::move(nAp); this->id_ = std::move(nid); this->T_ = std::move(nT);

	// reindex
	indexId_();

	return *this;
}


// conversion
fMat ll_cell::getAp(const aTv& ts) const noexcept {
	assert(std::all_of(ts.begin(),ts.end(),[this](const aT t){return validType(t);}));

	const auto N = Ntype(ts);
	fMat res(dim(),0); res.reserve(std::accumulate(N.cbegin(),N.cend(),aC(0)));
	for (const auto t: ts)
		res.push_back(getAp(t));
	return res;
}
ll_cell ll_cell::getSubCell(const aTv& ts) const noexcept {
	assert(std::all_of(ts.cbegin(),ts.cend(),[this](const aT t){return validType(t);}));	
	assert(std::is_sorted(ts.cbegin(),ts.cend()));
	assert(std::set<aT>(ts.cbegin(),ts.cend()).size()==ts.size());

	// new N and Ap	
	const auto nN = Ntype(ts);
	fMat nAp(dim(),0); nAp.reserve(std::accumulate(nN.begin(),nN.end(),0));
	for (auto t: ts)
		for (auto i=ccBegin(t),e=ccEnd(t); i!=e; ++i)
			nAp.push_back(*i);

	if (id().empty())
		return ll_cell(B(),std::move(nAp),std::move(nN));
	
	// new id
	idv nid; nid.reserve(ts.size());
	for (const auto t: ts) nid.push_back(id(t));
	
	return ll_cell(B_,std::move(nAp),std::move(nN),std::move(nid));
}
ll_bonds<> ll_cell::getBonds(const fMat& NN,
		const std::function<bool(const fMat&,const i_i&)>& keep) const noexcept {
		return ll_bonds<>(*this,NN,keep);
}
ll_bonds<i_i_R> ll_cell::getBonds(const aTv& T1, const aTv& T2,
				  const double f, const fMat& NN) const noexcept {
	// get indices for T1,T2
	auto I1 = inds(T1); std::sort(I1.begin(),I1.end());
	auto I2 = inds(T2); std::sort(I2.begin(),I2.end());

	// get distance matrix
	auto D = genDmat(B(),Ap().get({},I1),Ap().get({},I2),f,NN);
	for (auto& d: D) if(lm__::ops::z(d)) d=DBL_MAX;
	
	// linearize D and sort, find cutoff length l
	double l;
	{
		auto d = D.R(); d.push_back(mnorm(B()));
		std::sort(d.begin(),d.end());
		l = d[0];

		auto i=d.cbegin()+1,e=d.cend();
		if (f<-1.0)	while ((i!=e) && (*i < -f*l)) l=*i++;
		else		while ((i!=e) && (*i < f)) l=*i++;
	}

	return ll_bonds<i_i_R>(*this,NN,l,I1,I2);
}
ll_bonds<i_i_R> ll_cell::getBonds(const double f, const fMat& NN) const noexcept {
	assert(f<-1.0 || f>=0.0);
	if (empty()) return ll_bonds<i_i_R>();

	// get distance matrix
	auto D = genDmat(B(),Ap(),f,NN);
	
	// linearize D and sort, find cutoff length l
	double l;
	{
		auto d = D.upper().R(); d.push_back(mnorm(B()));
		std::sort(d.begin(),d.end());
		l = d[0];

		auto i=d.cbegin()+1,e=d.cend();
		if (f<-1.0)	while ((i!=e) && (*i < -f*l)) l=*i++;
		else		while ((i!=e) && (*i < f)) l=*i++;
	}

	// return bonds
	return ll_bonds<i_i_R>(*this,NN,l);
}

// comparison
fMat ll_cell::getPmat(const ll_cell& inp) const noexcept {
	if (empty() || inp.empty() || dim()!=inp.dim()) return fMat(dim(),0);
	
	auto res = zeros<fMat>(dim());
	for (auto i=B().ccBegin(),e=B().ccEnd(); i!=e; ++i) {
		const auto j = std::find(inp.B().ccBegin(),inp.B().ccEnd(),*i);
		if (size_t(j.i())==dim()) return fMat(dim(),0);
		res(j.i(),i.i()) = 1.0;
	}
	return res;
}
fMat ll_cell::getAvec(const ll_cell& inp, const fMat& P) const noexcept {
	if (empty() || inp.empty()) return fMat(dim(),0);
	if (size(P)!=size(B())) return fMat(dim(),0);
	if (dim()!=inp.dim()) return fMat(dim(),0);
	if (Ntype()!=inp.Ntype()) return fMat(dim(),0);
	if (B().prod(P)!=inp.B()) return fMat(dim(),0);

	const aT lft = leastFreqType();
	const auto col = P.prod(Ap().cAt(T_[lft]));

	// check shift lambda
	auto ckt = [this,&P,&inp](const fMat& sh) -> bool {
		for (const auto t: types())
			for (auto i=ccBegin(t),e=ccEnd(t); i!=e; ++i)
				if (inp.find(P.prod(*i)+sh,t)==inp.ccEnd(t)) return false;
		return true;
	};

	// find Avec by checking possible shifts in lft
	for (auto i=inp.ccBegin(lft), e=inp.ccEnd(lft); i!=e; ++i) {
		const auto res = *i - col;
		if (ckt(res)) return T(P).prod(res);
	}
	
	return fMat(dim(),0);
}
bool ll_cell::sameLattice(const ll_cell& inp) const noexcept {
	if (empty() && inp.empty()) return true;
	if (dim()!=inp.dim()) return false;

	const double r = inp.vol()/vol();
	if (std::abs(inp.N()-N()*r)>mtol()) return false;
	if (!validBasis(inp.B()) || !inp.validBasis(B())) return false;
	
	const auto tmp = r>=1.0 ? ll_cell(inp).changeBasis(B()): ll_cell(*this).changeBasis(inp.B());
	return !getAvec(tmp,getPmat(tmp)).empty();
}

// printing
std::string ll_cell::print(const bool direct, const size_t prec) const noexcept {
	if(empty()) return "";
	
	std::stringstream outp;
	
	outp << "cell contents\n1.0\n"
	     << T(B()).print(prec);
	if (!id().empty())
		outp << "\n" << id();
	outp << "\n" << Ntype()
	     << "\n" << (direct ? "Direct": "Cartesian")
	     << "\n" << (direct ? T(Ap()).print(prec): T(B().prod(Ap())).print(prec));

	return outp.str();
}
void ll_cell::printToFile(const std::string& fileName) const { printPOSCAR(fileName,*this); }

// index non unique ids
void ll_cell::indexId_() noexcept {
	if (id().empty()) return;
	assert(Apsorted_(types()));

	// get sorted type order with ids as is
	// first by ids, then number of positions then by positions
	auto T = types();
	const auto I = aux::sorted_order(T.cbegin(),T.cend(),
		[this](const aT t1, const aT t2) -> bool {
			if (id(t1)<id(t2)) return true;
			if (id(t1)>id(t2)) return false;
			if (Ntype(t1)<Ntype(t2)) return true;
			if (Ntype(t1)>Ntype(t2)) return false;
			for (auto i1=ccBegin(t1),e1=ccEnd(t1),i2=ccBegin(t2); i1!=e1; ++i1,++i2) {
				if ((*i1)<(*i2)) return true;
				if ((*i1)>(*i2)) return false;
			}
			return false;
		});

	// reorder containers if needed
	if (!std::is_sorted(I.cbegin(),I.cend())) {

		// number of types in old order
		auto Nt = Ntype();
		
		// find reordered nAp from I and Nt
		fMat nAp(dim(),0); nAp.reserve(N());
		for (const auto i: I) {
			const size_t pos = std::accumulate(Nt.cbegin(),Nt.cbegin()+i,size_t(0));
			std::for_each(this->ccBegin()+pos,this->ccBegin()+pos+Nt[i],
				[&nAp](const auto& j)->void{nAp.push_back(j);});
		}
		mat_ = std::move(nAp);

		// reorder nN and reconstruct T_
		aux::reorder(Nt.begin(),I);
		T_.clear(); T_.push_back(0);
		for (const auto i: Nt) T_.push_back(T_.back()+i);

		// reorder id_
		aux::reorder(id_.begin(),I);
	}

	// index non unique ids
	for (auto i=id_.begin(),e=id_.end(); i!=e;) {
		
		// current stripped id
		const std::string sid = softstripId(*i);
		
		// find range of duplicate ids
		size_t N=0;
		for (auto j=i; j!=e && sid==softstripId(*j); ++j,++N);

		// set to stripped if no duplicates
		if (N==1) { *i++ = sid; continue; }
		
		// index range
		for (size_t ind=0; ind<N; ++ind,++i)
			*i = appendIndex_(sid,ind+1);
	}
}


// detect file type
size_t ll_cell::detectFileType(const std::string& fileName) const {
	
	// check for wannier90 wout file
	{
		auto file = aux::openFile<std::ifstream>(fileName);

		std::string line;
		std::getline(file,line);
		std::getline(file,line);
		std::getline(file,line);
		std::getline(file,line);
		if (line.find("WANNIER90")!=std::string::npos)
			return 1;

		file.close();
	}

	// check for OMEN lattice file
	{
		// check if line has 5 integers on it
		const auto cl5i = [](const std::string& line) -> bool {
			const std::regex rgx("[\\s]+");
			
			size_t cnt=0;
			for (std::sregex_token_iterator
				i(line.begin(),line.end(),rgx,-1), e; i!=e; ++i) {
			
				const std::string s = i->str();
				if (s.empty()) continue;
				if (std::any_of(s.cbegin(),s.cend(),[](const char c)->bool
					{return !std::isdigit(c);}))
					return false;
				++cnt;
			}

			return cnt==5;
		};
		
		auto file = aux::openFile<std::ifstream>(fileName);

		std::string line;
		std::getline(file,line);
		if (cl5i(line))
			return 2;
	}
	
	// default POSCAR
	return 0;
}

// streaming
std::ostream& operator<<(std::ostream& os, const ll_cell& inp) noexcept {
	return (os<<inp.print());
}
