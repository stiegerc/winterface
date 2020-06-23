// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "ll_hbondss.h"
#include "ll_fn.h"
#include "ll_io.h"
#include "ll_hio.h"

using namespace ll__;
using namespace lm__;


// constructor
ll_hbondss::ll_hbondss(const ll_hbondss_input& inp, std::ostream& os) {
	assert(!inp.wbh.empty());
	
	// read wbhs
	this->vec_.reserve(inp.wbh.size());
	for (const auto& s: inp.wbh) {
		vec_.push_back(ll_hbonds(s));
		if (inp.verbosity & PRINTBIT__)
			os << "\nwannier bonds file '" << s << "' loaded containing "
			   << BLUE__ << vec_.back().cell().N() << RESET__ << " positions and "
			   << BLUE__ << vec_.back().cardinality() << RESET__ << " bonds";
	}
	if (std::any_of(this->cbegin(),this->cend(),[this](const auto& w)->bool
			{return w.dim()!=this->vec_.front().dim();}))
		throw(std::invalid_argument("dimensions on wbh don't match!"));
	if (inp.verbosity & PRINTBIT__)
		os << "\n";
}


// information from structure
ll_hbondss::pb ll_hbondss::getPerimeterConnections(const fMat& Ap,
	const std::vector<size_t>& I, const ll_hbondss_input& inp, std::ostream& os) const {
	
	assert(Ap.M()==dim());
	assert(Ap.N()==I.size());

	if (!inp.perimeter_radius) return {};
	
	// cutoff level for which bonds to include
	const double Rcut = inp.perimeter_radius<.0 ? -inp.perimeter_radius:
		this->radius()*(inp.perimeter_radius>1.0 ? 1.0: inp.perimeter_radius);
	
	if (inp.verbosity & PRINTBIT__)
		os << "\nsearching for approximate matches across section perimeters...\n"
		   << "using Rcut: " << BLUE__ << Rcut << RESET__;

	// get indices for each section, check basic properties
	std::vector<std::vector<size_t>> SI(size());
	for (auto& i: SI)
		i.reserve(I.size());
	for (size_t i=0; i!=I.size(); ++i) {
		const size_t j = j_(I[i]);
		if (j>=size() || i_(I[i])>=(*this)[j].cell().N())
			throw(std::invalid_argument("bad type"));
		SI[j].push_back(i);
	}
	if (std::any_of(SI.cbegin(),SI.cend(),
		[](const auto& i)->bool{ return i.empty(); }))
		throw(std::invalid_argument("structure appears incomplete"));


	// result index container, reserve estimate space
	std::vector<ll_hbondss::pbrt> RT; RT.reserve(Ap.N()/2);
	fMat bnds12(dim(),0);		 bnds12.reserve(Ap.N()/2);

	// find bonds across perimeter
	for (size_t j1=0   ; j1!=SI.size(); ++j1)
	for (size_t j2=j1+1; j2!=SI.size(); ++j2)
		for (const size_t si1: SI[j1])
		for (const size_t si2: SI[j2]) {

			// indices inside sections			
			const size_t i1 = i_(I[si1]);
			const size_t i2 = i_(I[si2]);

			// bond 1 -> 2
			const fMat bnd12 = Ap.cAt(si2)-Ap.cAt(si1);
			if (norm(bnd12)>Rcut) continue;
			
			// get iterators to beginning of index range
			auto itr = std::lower_bound(RT.cbegin(),RT.cend(),
				std::array<size_t,4>{{j1,i1,j2,i2}},
				[](const auto& i, const auto& j) -> bool {
					return std::lexicographical_compare(
						i.dat.cbegin(),i.dat.cbegin()+4,
						j.cbegin(),j.cend());
				}
			);
			auto jtr = bnds12.ccBegin() + std::distance(RT.cbegin(),itr);
			
			// check if index, bond combo is included already
			while (itr!=RT.cend() && (itr->j1!=j1 || itr->i1!=i1 ||
						  itr->j2!=j2 || itr->i2!=i2 || *jtr!=bnd12))
				++itr,++jtr;
			if (itr!=RT.cend() && itr->j1==j1 && itr->i1==i1 &&
					      itr->j2==j2 && itr->i2==i2 && *jtr==bnd12)
				continue;

			// check if interaction is found, 1 -> 2
			size_t a2=0, e2=(*this)[j1].cell().N();
			for (; a2!=e2; ++a2)
				if ((*this)[j1].getInteraction({i1,a2},bnd12)!=eH())
					break;
			
			// check if interaction is found, 2 -> 1
			const fMat bnd21 = -bnd12;
			size_t a1=0, e1=(*this)[j2].cell().N();
			for (; a1!=e1; ++a1)
				if ((*this)[j2].getInteraction({i2,a1},bnd21)!=eH())
					break;

			// one or both interactions not found
			if (a2==e2 || a1==e1)  continue;

			// insert new entry
			RT.insert(itr,pbrt{{j1,i1,j2,i2,a2,a1}});
			bnds12.cInsert(jtr,bnd12);
		}


	// check section connectivity
	if (RT.front().j1)
		throw(std::invalid_argument("get perimeter connections: bad section connectivity"));
	std::vector<bool> conn(size(),false); conn.front() = true;
	for (const auto& rt: RT) {
		if (conn[rt.j1]) conn[rt.j2] = true;
		if (conn[rt.j2]) conn[rt.j1] = true;
	}
	if (!std::all_of(conn.cbegin(),conn.cend(),[](const bool i)->bool{return i;}))
		throw(std::invalid_argument("get perimeter connections: bad section connectivity"));


	// all good, print and return
	if (inp.verbosity & PRINTBIT__)
		os << ", found " << BLUE__ << RT.size() << RESET__ << " bond pairs\n";
	
	return {std::move(RT),std::move(bnds12)};
}
std::vector<std::vector<i_i>> ll_hbondss::getSubstitutes(const fMat& Ap,
			const std::vector<size_t>& I, const double f) const noexcept {
	
	std::vector<std::vector<i_i>> res;
	
	std::vector<size_t> cI(I.size());
	fMat cp(Ap.M(),I.size());
	for (size_t j=0; j!=this->size(); ++j) {
		cI.resize(0); cp.resize(0);

		auto p = Ap.ccBegin();
		for (auto i=I.cbegin(),e=I.cend(); i!=e; ++i,++p)
			if (j_(*i)==j)
				cI.push_back(i_(*i)),
				cp.push_back(*p);

		res.push_back((*this)[j].getSubstitutes(cp,cI,f));
	}

	return res;
}


// modification
void ll_hbondss::equalizeShifts(const ll_hbondss::pb& B,
			const ll_hbondss_input& inp, std::ostream& os) {
	assert(B.RT.size()==B.bnds12.N());
	
	if (B.RT.empty()) return;
	
	if (inp.verbosity & PRINTBIT__)
		os << "\nequalizing energy shifts...";

	std::vector<size_t> leveled = {0};
	while (leveled.size() < size()) {
		
		auto q=B.RT.cbegin(),qe=B.RT.cend(); auto b = B.bnds12.ccBegin();
		while (q!=qe) {
			
			// find range of this perimeter
			const auto cq=q; const auto cb=b;
			while (q!=qe && q->j1==cq->j1 && q->j2==cq->j2)
				++q,++b;

			// check if exactly one of the sections involved is leveled already
			const bool f1 = std::binary_search(leveled.cbegin(),leveled.cend(),cq->j1);
			const bool f2 = std::binary_search(leveled.cbegin(),leveled.cend(),cq->j2);
			if (f1+f2 != 1) continue;

			// get number of total orbitals
			size_t Norb = 0;
			for (auto ccq=cq; ccq!=q; ++ccq)
				Norb += (*this)[cq->j1].Norb(ccq->i1);

			// list of matching already included
			std::vector<i_i> iincl; iincl.reserve(std::distance(cq,q));

			// gather differences along diagonals
			fMat dd(1,0); dd.reserve(Norb);
			std::for_each(cq,q,[this,&dd,&iincl,&inp,&os](const auto& rt) -> void {
			{
				const i_i cp{rt.a2,rt.i2};
				const auto itr = std::lower_bound(iincl.cbegin(),
									   iincl.cend(),cp);
				if (itr == iincl.cend() || *itr!=cp) {
					iincl.insert(itr,cp);
					dd.push_back(
					real(diag(*(*this)[rt.j1].search({rt.a2,rt.a2})->center())).T()
				       -real(diag(*(*this)[rt.j2].search({rt.i2,rt.i2})->center())).T());
					
					if (inp.verbosity & DEBUGBIT__)
						os << "\n including {'"
					   	   << (*this).id(rt.j1,rt.a2) << "' - '"
					   	   << (*this).id(rt.j2,rt.i2) << "'}";
				}
			}
			{
				const i_i cp{rt.i1,rt.a1};
				const auto itr = std::lower_bound(iincl.cbegin(),
									   iincl.cend(),cp);
				if (itr == iincl.cend() || *itr!=cp) {
					iincl.insert(itr,cp);
					dd.push_back(
					real(diag(*(*this)[rt.j1].search({rt.i1,rt.i1})->center())).T()
				       -real(diag(*(*this)[rt.j2].search({rt.a1,rt.a1})->center())).T());
					
					if (inp.verbosity & DEBUGBIT__)
						os << "\n including {'"
					   	   << (*this).id(rt.j1,rt.i1) << "' - '"
					   	   << (*this).id(rt.j2,rt.a1) << "'}";
				}
			}
			});

			// shift section 2 to match section 1
			if (f1) {
				const double m = mean(dd); const double d = stdd(dd,m);
				if (inp.verbosity & PRINTBIT__)
					os << "\nshifting sections " << BLUE__ << cq->j1 << RESET__
					   << " < " << BLUE__ << cq->j2 << RESET__ << " by "
					   << BLUE__ << m << RESET__ << " ("
					   << BLUE__ << d << RESET__ << ")";
				this->vec_[cq->j2].shiftSpectrum(m);
				leveled.insert(std::lower_bound(leveled.cbegin(),leveled.cend(),cq->j2),cq->j2);
			} else {
				const double m = mean(dd); const double d = stdd(dd,m);
				if (inp.verbosity & PRINTBIT__)
					os << "\nshifting sections " << BLUE__ << cq->j2 << RESET__
					   << " < " << BLUE__ << cq->j1 << RESET__ << " by "
					   << BLUE__ << -m << RESET__ << " ("
					   << BLUE__ <<  d << RESET__ << ")";
				this->vec_[cq->j1].shiftSpectrum(-m);
				leveled.insert(std::lower_bound(leveled.cbegin(),leveled.cend(),cq->j1),cq->j1);
			}
		}
	}
	
	if (inp.verbosity & PRINTBIT__)
		os << "\n";
}
void ll_hbondss::equalizeSigns(const ll_hbondss::pb& B, const std::vector<std::vector<i_i>>& SI,
			const ll_hbondss_input& inp, std::ostream& os) {
	assert(B.RT.size()==B.bnds12.N());
	
	if (inp.verbosity & PRINTBIT__)
		os << "\nequalizing wannier signs...";

	std::vector<size_t> leveled = {0};
	while (leveled.size() < size()) {
		
	auto q=B.RT.cbegin(),qe=B.RT.cend(); auto b = B.bnds12.ccBegin();
	while (q!=qe) {
		
		// find range of this perimeter
		const auto cq=q; const auto cb=b;
		while (q!=qe && q->j1==cq->j1 && q->j2==cq->j2)
			++q,++b;

		// check if exactly one of the sections involved is leveled already
		const bool f1 = std::binary_search(leveled.cbegin(),leveled.cend(),cq->j1);
		const bool f2 = std::binary_search(leveled.cbegin(),leveled.cend(),cq->j2);
		if (f1+f2 != 1) continue;
		if (!f1 && inp.verbosity & PRINTBIT__)
			os << "\nadapting section " << BLUE__ << cq->j2
			   << RESET__ << " <- " << BLUE__ << cq->j1 << RESET__;
		if (!f2 && inp.verbosity & PRINTBIT__)
			os << "\nadapting section " << BLUE__ << cq->j1
			   << RESET__ << " <- " << BLUE__ << cq->j2 << RESET__;

	
		// get unique indices for each section
		std::vector<size_t> I1; I1.reserve(std::distance(cq,q));
		std::vector<size_t> I2; I2.reserve(std::distance(cq,q));
		std::for_each(cq,q,[&I1,&I2](const auto& rt) -> void {
			const auto itr1 = std::lower_bound(I1.cbegin(),I1.cend(),rt.i1);
			if (itr1==I1.cend() || *itr1!=rt.i1)
				I1.insert(itr1,rt.i1);
			const auto itr2 = std::lower_bound(I1.cbegin(),I1.cend(),rt.a2);
			if (itr2==I1.cend() || *itr2!=rt.a2)
				I1.insert(itr2,rt.a2);
			const auto itr3 = std::lower_bound(I2.cbegin(),I2.cend(),rt.i2);
			if (itr3==I1.cend() || *itr3!=rt.i2)
				I2.insert(itr3,rt.i2);
			const auto itr4 = std::lower_bound(I2.cbegin(),I2.cend(),rt.a1);
			if (itr4==I1.cend() || *itr4!=rt.a1)
				I2.insert(itr4,rt.a1);
		});

		// get indices and references to treat the larger and smaller
		// index sets I1,I2 as *b (larger) and *s (smaller)
		const size_t ib = I1.size()>I2.size() ? 0: 1;
		const size_t is = I1.size()<I2.size() ? 0: 1;
		const size_t jb = I1.size()>I2.size() ? cq->j1: cq->j2;
		const size_t js = I1.size()<I2.size() ? cq->j1: cq->j2;
		const auto& Ib  = I1.size()>I2.size() ? I1: I2;
		const auto& Is  = I1.size()<I2.size() ? I1: I2;

		// find mapping of the section with less unique indices
		// onto those of the other section
		std::vector<size_t> MI; MI.reserve(Ib.size());
		for (const auto i: Ib) {
			const auto& rt = *std::find_if(cq,q,[i,ib](const auto& rt) -> bool {
				return rt(ib,1) == i || rt(2,ib) == i;
			});
			if (rt(ib,1) == i) MI.push_back(rt(2,is));
			if (rt(2,ib) == i) MI.push_back(rt(is,1));
		}
		if (inp.verbosity & VERBOBIT__) {
			os << "\n type mapping of section " << BLUE__ << jb
			   << RESET__ << " <- " << BLUE__ << js << RESET__;
			for (size_t i=0; i!= Ib.size(); ++i)
				os << "\n  '" << (*this).id(jb,Ib[i]) << "' <- '"
				   << (*this).id(js,MI[i]) << "'";
		}
		
		// find the signs matrix from hamiltonian test matrices
		// for each section
		fMat S;
		{
			// atomic positions for section with more indices
			const fMat Ap = (*this)[jb].cell().getcAp(Ib);
		
			// find total number of orbitals for all atoms involved
			const size_t Norb = std::accumulate(Ib.cbegin(),Ib.cend(),size_t(0),
				[this,jb](const size_t acc, const size_t i) -> size_t {
					return acc + (*this)[jb].Norb(i);
			});

			fMat Hb = zeros<fMat>(Norb,Norb);
			hctor((*this)[jb],eye<fMat>(dim()),Ap,Ib,ll_writerR0<fMat>(Hb),DBL_MAX,false,1);	
			if (inp.verbosity & DEBUGBIT__) {
				const std::string fileName = inp.prefix+"Hb_"
					+std::to_string(jb)+"_"+std::to_string(js)+".mat";
				Hb.writeToFile(fileName);
				if (inp.verbosity & VERBOBIT__)
					os << "\nfile '" << fileName << "' written";
			}
			fMat Hs = zeros<fMat>(Norb,Norb);
			hctor((*this)[js],eye<fMat>(dim()),Ap,MI,ll_writerR0<fMat>(Hs),DBL_MAX,false,1);
			if (inp.verbosity & DEBUGBIT__) {
				const std::string fileName = inp.prefix+"Hs_"
					+std::to_string(jb)+"_"+std::to_string(js)+".mat";
				Hs.writeToFile(fileName);
				if (inp.verbosity & VERBOBIT__)
					os << "\nfile '" << fileName << "' written";
			}
		
			// get signs matrix using signs_tol
			set_mtol(inp.signs_tol);
			S = sign(Hb)*sign(Hs);
			reset_mtol();
		}
		if (inp.verbosity & DEBUGBIT__) {
			const std::string fileName = inp.prefix+"S_"
				+std::to_string(jb)+"_"+std::to_string(js)+".mat";
			S.writeToFile(fileName);
			if (inp.verbosity & VERBOBIT__)
				os << "\nfile '" << fileName << "' written";
		}

		// find sign flip vector P
		fMat P;
		{
			// extract relevant equations from S
			fMat A(S.M(),0); A.reserve(S.N()); size_t rk=0;
			for (size_t n=0; n!=S.N(); ++n)
			for (size_t m=n+1; m!=S.M(); ++m)
				if (S(m,n)) {
					fMat ne = zeros<fMat>(S.M(),1);
					ne[n] = 1.0, ne[m] = -S(m,n);
					A.push_back(ne);
					
					if (rank(A)>rk) ++rk;
					else A.pop_back();
				}

			// add identities until the rank of A is full
			fMat b = zeros<fMat>(S.M(),1);
			for (size_t m=0; m!=S.M() && rk<S.M(); ++m) {
				fMat ne = zeros<fMat>(S.M(),1); ne[m] = 1.0;
				A.push_back(ne);

				if (rank(A)>rk) b[rk]=1.0, ++rk;
				else A.pop_back();
			}

			P = T(A).leftDivide(b);
		}
		if (std::any_of(P.cbegin(),P.cend(),
		[](const double p)->bool{return !p;}))
			throw(std::invalid_argument("solution to signs equations invalid,"
						    " try increasing tolerance"));
		if (inp.verbosity & DEBUGBIT__) {
			const std::string fileName = inp.prefix+"P_"	
				+std::to_string(jb)+"_"+std::to_string(js)+".mat";
			diag(P).writeToFile(fileName);
			if (inp.verbosity & VERBOBIT__)
				os << "\nfile '" << fileName << "' written";
		}

		// decompose P into smaller vectors corresponding to the indices
		std::vector<std::vector<bool>> dP;
		{
			auto ip = P.cbegin();
			for (const size_t i: Ib) {
				const size_t Norb = (*this)[jb].Norb(i);
				dP.push_back(std::vector<bool>(Norb));
				std::transform(ip,ip+Norb,dP.back().begin(),
					[](const double p)->bool{return p<.0;});
				ip += Norb;
			}
		}

		// check consistency with duplicate indices from the section
		// with less unique indices
		for (const size_t i: Is) {
			// find indices where MI == i
			std::vector<size_t> J; J.reserve(MI.size());
			for (size_t j=0; j!=MI.size(); ++j)
				if (MI[j]==i) J.push_back(i);

			// check corresponding entries in DP are identical
			if (std::any_of(J.cbegin(),J.cend(),
					[&dP,&J](const size_t j) -> bool {
						return dP[J.front()] != dP[j];
					}))
				throw(std::invalid_argument("inconsistent "
				"wannier signs, try increasing tolerance"));
		}		

		// tracker for inds flipped, and recursive flip matched inds lambda
		std::vector<size_t> fI; fI.reserve(std::max(I1.size(),I2.size())
							   +2*SI[cq->j1].size());
		std::function<void(const size_t, const size_t,
				   const std::vector<bool>&)> flipper;
		flipper = [this,&fI,&SI,&flipper,&inp,&os](const size_t j,
				const size_t i, const std::vector<bool>& p) -> void {
			
			const auto itr = std::lower_bound(fI.cbegin(),fI.cend(),i);
			if (itr!=fI.cend() && *itr == i)
				return; // already flipped

			// new ind, insert into fI
			fI.insert(itr,i);

			// flip signs for this ind
			if (inp.verbosity & VERBOBIT__)
				os << "\n  " << std::setw(12)
				   << (+"'"+(*this).id(j,i)+"'")
				   << "' < " << p;
			this->vec_[j].flipWannierSigns(i,p);

			// call recursively for all matched inds to this ind
			for (auto& si: SI[j]) {
				if (si.i1() == i) flipper(j,si.i2(),p);
				if (si.i2() == i) flipper(j,si.i1(),p);
			}
		};


		// all good, flip wannier signs!
		if (inp.verbosity & VERBOBIT__)
			os << "\nflipping wannier signs...";

		
		if (!f1) {
			// leveling section 1
			if (I1.size()>I2.size()) {
				auto i = I1.cbegin();
				for (auto idp=dP.cbegin(),edp=dP.cend(); idp!=edp; ++idp,++i)
					flipper(cq->j1,*i,*idp);
			} else
				for (const size_t i: I1) {
					const size_t pos = std::distance(
					MI.cbegin(),std::find(MI.cbegin(),MI.cend(),i));
					flipper(cq->j1,i,dP[pos]);
				}
			
			leveled.insert(std::lower_bound(leveled.cbegin(),leveled.cend(),
				cq->j1),cq->j1);
		} else {
			// leveling section 2
			if (I2.size()>I1.size()) {
				auto i = I2.cbegin();
				for (auto idp=dP.cbegin(),edp=dP.cend(); idp!=edp; ++idp,++i)
					flipper(cq->j2,*i,*idp);
			} else
				for (const size_t i: I2) {
					const size_t pos = std::distance(
					MI.cbegin(),std::find(MI.cbegin(),MI.cend(),i));
					flipper(cq->j2,i,dP[pos]);
				}
			
			leveled.insert(std::lower_bound(leveled.cbegin(),leveled.cend(),
				cq->j2),cq->j2);
		}
	}
	}
	
	if (inp.verbosity & PRINTBIT__)
		os << "\n";
}


// approximate searching in multiple hbonds
ll_hbonds::am_ ll_hbondss::getApproximateInteraction(const i_i& j, const fArray& b) const noexcept {

	const auto I1 = decompInd_(j.i1()), I2 = decompInd_(j.i2());
	if (I1.j>=this->size() || I2.j>=this->size())
		return {eH(),NPOS__,NPOS__};

	cMat H12; size_t i2;
	for (i2=0; i2!=(*this)[I1.j].cell().N(); ++i2)
		if ((*this)[I1.j].Norb(i2) == (*this)[I2.j].Norb(I2.i)) {
			H12 = (*this)[I1.j].getInteraction({I1.i,i2},b);
			if (!H12.empty()) break;
		}

	cMat H21; size_t i1;
	for (i1=0; i1!=(*this)[I2.j].cell().N(); ++i1)
		if ((*this)[I2.j].Norb(i1) == (*this)[I1.j].Norb(I1.i)) {
			H21 = T((*this)[I2.j].getInteraction({I2.i,i1},-b));
			if (!H21.empty()) break;
		}

	if (!H12.empty()) {
		if (!H21.empty()) return {(H12+H21)*.5, compInd_(I1.j,i2), compInd_(I2.j,i1)};
		else		  return {H12*.5,       compInd_(I1.j,i2), NPOS__};
	}
	if (!H21.empty())	  return {H21*.5,       NPOS__,            compInd_(I2.j,i1)};
	else			  return {eH(),         NPOS__,            NPOS__};
}
