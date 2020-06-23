// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "ll_hbonds.h"
#include "lm_fn.h"
#include "ll_fn.h"
#include "ll_io.h"
#include "ll_lambda.h"
#include "aux_md5.h"
#include <fstream>
#include <limits>

using namespace lm__;
using namespace ll__;

// matching constructor
ll_hbonds::ll_hbonds(const ll_wmatching_input& inp, std::ostream& os) {
	
	auto w = readWp(inp.wout,inp.wbl);

	ll_cell cell; wi I;

	set_mtol(WTOL__); // set tolerance due to wannier input precision
	try {
	
	if (inp.verbosity & PRINTBIT__) {
		os << "generating material using matching mode '" << inp.mode << "'\n"
		   << "wannier data from files:\n"
		   << " '" << inp.wout << "'"
		   << (inp.verbosity & MD5BIT__ ? ", "+aux::md5(inp.wout): "") << "\n"
		   << " '" << inp.hrdat << "'"
		   << (inp.verbosity & MD5BIT__ ? ", "+aux::md5(inp.hrdat): "") << "\n"
		   << "using tolerance: " << BLUE__ << std::setprecision(4)
		   << std::scientific << inp.tol << RESET__ << "\n";
		if (inp.abl.size()) os << " abl = " << inp.abl;
		if (inp.wbl.size()) os << " wbl = " << inp.wbl;
	}

	// read initial restriction vector, basis B and Wp
	rv ir(hrDim(inp.hrdat));				// initial r
	fMat B = orthogonalize(readB(inp.wout),ir);		// get orth. B
	B.leftDivideEq(w.Wp());					// force Wp to basis B
	for (const auto i: ll__::inds(ir))			// force into positive section
		w.Wp().rAt(i) %= 1.0;				// in restricted directions

	// get cell and wannier matching
	switch (fnvHash(inp.mode.front().c_str())) {
		case "cluster"_h:
			I = clusterize(B,w.Wp(),inp.wradius,inp.minpts,NN);
			cell = ll_cell(B,genCenters(B,w.Wp(),w.s(),I),aCv(I.size(),1),"CC");
		break;
		case "wannier"_h:
			cell = ll_cell(B,w.Wp()%1.0,aCv(w.N(),1),"WW");
			I.reserve(w.N());
			for (size_t i=0; i!=w.N(); ++i)
				I.push_back({i});
		break;
		case "atomic"_h:
		{
			cell = ll_cell(inp.wout).diversify();
			I = matchCenters(cell,w.Wp(),NN);
		}
		break;
		case "atomic_bc"_h:
		{
			cell = ll_cell(inp.wout).diversify();
			cell.merge(cell.getBonds(inp.bond_factor,NN).getBondCenters());
			I = matchCenters(cell,w.Wp(),NN);
		}
		break;
		case "composite"_h:
		{
			if (inp.mode.size()==1)
				throw(std::invalid_argument("composite mode: "
							"additional entries expected"));
			// unit cell from wannier90 output
			const ll_cell mCell(inp.wout);
			
			// lambda for basic property check
			const auto bpc = [](const std::string& s) -> void {
				// check total weight of braces is zero
				if (std::count(s.cbegin(),s.cend(),'(') !=
				    std::count(s.cbegin(),s.cend(),')'))
					throw(std::invalid_argument(
					"composite mode: bad entry"));
					
				// check there is at most one central ':'
				size_t cnt=0;
				for (size_t pos=s.find(':'); pos!=std::string::npos; pos=s.find(':',pos+1)) {
					if (std::count(s.cbegin(),s.cbegin()+pos,'(') ==
					    std::count(s.cbegin(),s.cbegin()+pos,')'))
						++cnt;
					if (cnt>1)
						throw(std::invalid_argument(
						"composite mode: bad entry"));
				}
			};

			// recursive lambda to parse entries
			std::function<ll_cell(const std::string& s)> pe;
			pe = [&mCell,&pe,bf=inp.bond_factor](const std::string& s) -> ll_cell {

				if (s.empty())
					throw(std::invalid_argument(
					"composite mode: bad entry"));
				
				// look for central ':'
				size_t pos;
				for (pos=s.find(':'); pos!=std::string::npos &&
					(std::count(s.cbegin(),s.cbegin()+pos,'(') !=
					 std::count(s.cbegin(),s.cbegin()+pos,')'));
					pos=s.find(':',pos+1));

				// no central ':'
				if (pos == std::string::npos) {
					// braces around whole string
					if (s.front()=='(' && s.back()==')')
						return pe(s.substr(1,s.size()-2));
					
					// no braces around, no central ':' -> atomic type
					const aT t = mCell.type(s);
					if (!mCell.validType(t))
						throw(std::invalid_argument(
						"composite mode: bad entry"));
					return mCell.getSubCell({t}).diversify();
				} 
				
				// found central ':' -> bond type
				const std::string s1 = s.substr(0,pos);
				const std::string s2 = s.substr(pos+1);
				const ll_cell CC = pe(s1).merge(pe(s2));
				
				return CC.getBonds(CC.strippedType(s1),
						   CC.strippedType(s2),bf,NN).
					getBondCenters().diversify();
			};
			

			cell = ll_cell(mCell.B());
			for (auto i=inp.mode.cbegin()+1,e=inp.mode.cend(); i!=e; ++i)
				bpc(*i), cell.merge(pe(*i));
			I = matchCenters(cell,w.Wp(),NN);
		}
		break;
		default:
			cell = ll_cell(inp.mode.front()).diversify();
			I = matchCenters(cell.B(),cell.Ap(),w.Wp(),NN);
		break;
	}

	// simplify composite ids
	if (inp.simplify_composite || std::any_of(cell.id().cbegin(),cell.id().cend(),
			[](const std::string& s)->bool{return s.size()>256;}))
		cell.simplifyCompositeIds();

	// remove umatched
	if (std::any_of(I.cbegin(),I.cend(),[](const auto& i)->bool{return i.empty();})) {

		// get indices to remove
		aTv RMT; RMT.reserve(cell.N());
		if (inp.rm_unmatched & 1)
			for (const auto t: cell.compositeTypes())
				if (I[t].empty()) RMT.push_back(t);
		if (inp.rm_unmatched & 2)
			for (const auto t: cell.fundamentalTypes())
				if (I[t].empty()) RMT.push_back(t);
		std::sort(RMT.begin(),RMT.end());

		fMat Ap_ = cell.moveAp();
		idv id_ = cell.moveId();

		// remove from the back
		for (auto i=RMT.crbegin(),e=RMT.crend(); i!=e; ++i)
			Ap_.cRm(*i), id_.erase(id_.begin()+(*i)), I.erase(I.begin()+(*i));

		// reconstruct cell
		cell = ll_cell(cell.B(),std::move(Ap_),aCv(Ap_.N(),1),std::move(id_));
	}

	// balance matching numbers
	if (inp.balance_matches)
		balance(I,cell,w.Wp()%1.0,NN);

	// center positions in restricted dimensions
	{
		const auto sh = .5-com(cell.Ap());
		cell.shift(sh*fMat(ir));
		for (const auto i: ll__::inds(ir))
			w.Wp().rAt(i)+=sh[i], w.Wp().rAt(i)%=1.0;
	}

	if (inp.verbosity & PRINTBIT__) {
		// # of positions with wannier centers matched to them
		const size_t Nm = std::accumulate(I.cbegin(),I.cend(),size_t(0),
			[](const size_t s, const auto& i)->size_t{return s+(bool)i.size();});
		
		os << "\nfound " << w.N() << " wannier centers, matched to "
		   << BLUE__ << Nm << "/" << cell.N() << RESET__ << " positions\n";
	
		// distance matrix positions vs wannier centers, bond length
		fMat D(1,0); D.reserve(w.N());
		for (auto i=cell.ccBegin(),e=cell.ccEnd(); i!=e; ++i)
			if (!I[(size_t)i].empty()) D.push_back(
				genDmat(cell.B(),*i,w.Wp().get({},I[(size_t)i])%1.0,NN));
		const double bl = cell.getSubCell(cell.fundamentalTypes()).
			getBonds(inp.bond_factor,NN).radius();

		const double mean_d = mean(D), max_d = max(D);
		os << " mean: " << std::setprecision(3) << std::fixed
		   << (mean_d<bl ? GREEN__: RED__) << mean_d << RESET__ << " [Angstrom]\n"
		   << " max:  " << std::setprecision(3) << std::fixed
		   << (max_d<bl ? GREEN__: RED__) << max_d << RESET__ << " [Angstrom]\n"
		   << " tot: " << std::setprecision(3) << std::fixed
		   << BLUE__ << sum(D) << RESET__ << " [Angstrom]\n";
	}
	if (inp.verbosity & WRITEBIT__) {
		// write inital wannier positions
		printPOSCAR(inp.prefix+"iwp.psc",cell.B(),w.Wp(),wiToT(I),{},1.0,true);
		if (inp.verbosity & VERBOBIT__) os << "wannier centers POSCAR file '"
			<< inp.prefix << "iwp.psc' written"
			<< (inp.verbosity & MD5BIT__ ? ", "+md5(inp.prefix+"iwp.psc"): "") << "\n";

		// write initial wannier positions, one file for each matched center
		if (inp.verbosity & DEBUGBIT__)
		for (size_t i=0; i!=I.size(); ++i) {
			if (I[i].empty()) continue;
			
			// positions, type vector and stripped id
			fMat Ap(DIM__,0);
			Ap.push_back(cell.Ap().cAt(i));
			Ap.push_back(w.Wp().get({},I[i]));
			aTv T(I[i].size()+1,1); T[0]=0;

			// distance matrix
			const auto D = genDmat(cell.B(),Ap.cFront(),w.Wp().get({},I[i])%1.0,NN);

			std::stringstream sstr;
			sstr << inp.prefix << "iwp_" << std::setfill('0') << std::setw(4) << (i+1)
				<< "_[" << ll_cell::stripId(cell.id(cell.type(i)))
				<< "]_" << std::setfill('0') << std::setw(2) << I[i].size()
				<< "_(" << std::fixed << std::setw(4) << std::setprecision(3)
				<< std::setfill('0') << mean(D)
				<< "," << std::fixed << std::setw(4) << std::setprecision(3)
				<< std::setfill('0') << max(D) << ").psc";

			printPOSCAR(sstr.str(),cell.B(),Ap,T,idv{"U","H"});
		}
		
		// write initial cell
		printPOSCAR(inp.prefix+"icell.psc",cell);
		if (inp.verbosity & VERBOBIT__) os << "cell POSCAR file '"
			<< inp.prefix << "icell.psc' written"
			<< (inp.verbosity & MD5BIT__ ? ", "+md5(inp.prefix+"icell.psc"): "") << "\n";

		// write initial cell, one file for each equal stripped types
		if (inp.verbosity & DEBUGBIT__)
		for (const auto& i: cell.equalTypes()) {
			const auto tmp = cell.getSubCell(i).collectivize();
			const std::string sid = tmp.stripId(tmp.id(0));
			printPOSCAR(inp.prefix+"icell_["+sid+"].psc",
					tmp.B(),tmp.Ap(),aTv(tmp.N(),0),idv{sid});
		}
		
		printWiToFile(inp.prefix+"wi",I);
		if (inp.verbosity & VERBOBIT__) os << "wannier matching file '"
			<< inp.prefix << "wi' written"
			<< (inp.verbosity & MD5BIT__ ? ", "+md5(inp.prefix+"wi"): "") << "\n";
	}

	} catch(const std::invalid_argument& e) {
		reset_mtol();
		throw(e);
	}

	// keep lambda, determines data to discard
	const auto keep = [tol=inp.tol](const cMat& inp)->bool
		{ return std::any_of(inp.cbegin(),inp.cend(),
			[tol](const hel& h)->bool{return cmph(h,tol);});};
	
	*this = ll_hbonds(cell,w.Wp(),I,readHr(inp.hrdat,inp.wbl,keep),keep);
	
	if (inp.verbosity & PRINTBIT__)
		os << "\ngenerated wannier bonds hamiltonian containing "
		   << BLUE__ << this->cardinality() << RESET__ << " bonds total\n";
}


// from cell constructor
ll_hbonds::ll_hbonds(ll_cell cell, const fMat& Wp, const wi& I,
			const R_H& hr, const std::function<cMat(cMat&&)>& transf) noexcept:
		ll_bonds(std::move(cell)) {
	if (hr.empty()) return;
	assert(!this->cell_.empty());
	assert(this->cell_.N()==this->cell_.Nspecies());
	assert(!Wp.empty());
	assert(this->cell_.dim()==Wp.M());
	assert(I.size()==this->cell_.Nspecies());
	assert(std::all_of(I.cbegin(),I.cend(),[&Wp](const wTv& i)->bool{
		return std::all_of(i.cbegin(),i.cend(),[&Wp](const wT j)->bool{
		return j<Wp.N();});}));
	assert(hr.dim()==Wp.M());
	assert(hr.front().M()==Wp.N());

	set_mtol(WTOL__);

	// generate Norb_
	Norb_.reserve(I.size());
	for (const auto& i: I)
		Norb_.push_back(i.size());

	// generate "wannier bonds"
	fMat Wb(Wp.M(),Wp.N());
	for (size_t i=0; i!=I.size(); ++i)
		for (const auto j: I[i]) Wb.cAt(j) = distb(this->cell().B(),this->cell().cAt(i),
							Wp.cAt(j)%1.0,NN).bonds.cFront();

	// wannier functions can be periodic in vacuum directions
	// ==> find vacuum in according to hr.R, i.e. row with all 0.0
	const auto VFIX = nany(hr.R().neq(0.0));

	// lambda to find hamiltonian block for an atomic bond
	const auto findBlock = [&Wp,&Wb,&hr,&I,&VFIX](const auto& p1, const aT t1,
				    const auto& p2, const aT t2) -> cMat {
		assert(t1<I.size());
		assert(t2<I.size());
		
		auto res = zeros<cMat>(I[t1].size(),I[t2].size());
		for (size_t i=0; i!=res.M(); ++i) {
			const wT wt1 = I[t1][i];
			const fMat R1 = p1+Wb.cAt(wt1)-Wp.cAt(wt1);
			
			for (size_t j=0; j!=res.N(); ++j) {
				const wT wt2 = I[t2][j];
				
				// get dR and force vacuum direction to 0.0
				const fMat dR = (p2+Wb.cAt(wt2)-Wp.cAt(wt2) - R1)*VFIX;
				
				const auto itr = std::lower_bound(
					hr.ccBegin(),hr.ccEnd(),dR,vcmp);
				if (itr!=hr.ccEnd() && *itr==dR)
					res(i,j) = hr[(size_t)itr](wt1,wt2);
			}
		}
		
		return res;
	};

	// find bonds, hamiltonians, indices, Nbonds_
	vec_.reserve(cell.N()*cell.N());
	for (auto p1=this->cell().ccBegin(),pe=this->cell().ccEnd(); p1!=pe; ++p1)
		for (auto p2=this->cell().ccBegin(); p2!=pe; ++p2) {
			vec_.push_back(i_i_R_H(size_t(p1),size_t(p2),this->cell_.dim()));

			vec_.back().reserve(hr.N());
			for (auto r=hr.ccBegin(),re=hr.ccEnd(); r!=re; ++r) {
				const auto H = transf(findBlock(*p1,vec_.back().i1(),
						   *p2+*r,vec_.back().i2()));
				if (!H.empty())
					vec_.back().push_back(*r,std::move(H));
			}

			if (!vec_.back().empty())
				vec_.back().shrink_to_fit();
			else    vec_.pop_back();
		}
	vec_.shrink_to_fit();

	reset_mtol();
}
ll_hbonds::ll_hbonds(const std::string& fileName): ll_bonds() {
	
	// open file
	auto file = aux::openFile<std::ifstream>(fileName, std::ios::binary);
	
	// read header
	{
		char head[5];
		file.read((char*) &head, 5*sizeof(char));
		if (strncmp(head,"wad90",5)) throw(std::runtime_error("bad header in file '"+fileName+"'"));
	}

	// read cell
	uint32_t D, N;
	file.read((char*) &D, sizeof(uint32_t));
	file.read((char*) &N, sizeof(uint32_t));
	{
		fMat B(D,D), Ap(D,N);
		file.read((char*)  B.data(), D*D*sizeof(double));
		file.read((char*) Ap.data(), D*N*sizeof(double));
		
		idv id; id.reserve(N);
		uint8_t l; char buff[std::numeric_limits<uint8_t>::max()+1];
		while (id.size()<id.capacity()) {
			file.read((char*) &l, sizeof(uint8_t));
			file.read((char*) buff, l*sizeof(char));
			id.push_back(std::string(buff));
		}
		
		set_mtol(WTOL__);
		this->cell_ = ll_cell(std::move(B),std::move(Ap),std::move(id));
		reset_mtol();
	}

	// read Norb
	Norb_.reserve(N);
	while (Norb_.size()<Norb_.capacity()) {
		file.read((char*) &N, sizeof(uint32_t));
		Norb_.push_back(N);
	}

	// read #pairs
	{
		uint32_t Np;
		file.read((char*) &Np, sizeof(uint32_t));
		vec_.reserve(Np);
	}

	// read pairs, R and H
	while (vec_.size()<vec_.capacity()) {
		uint32_t tmp[3]; // i1,i2,size
		file.read((char*) tmp, 3*sizeof(uint32_t));

		fMat R(D,tmp[2]);
		file.read((char*) R.data(), R.size()*sizeof(double));

		std::vector<cMat> H; H.reserve(tmp[2]);
		while (H.size()<H.capacity()) {
			cMat cH(Norb_[tmp[0]],Norb_[tmp[1]]);
			file.read((char*) cH.data(), cH.size()*sizeof(std::complex<double>));
			H.push_back(cH);
		}

		vec_.push_back(i_i_R_H(tmp[0],tmp[1],std::move(R),std::move(H)));
	}

	// close file
	file.close();
}

// information
std::vector<i_i> ll_hbonds::getSubstitutes(const fMat& Ap, const std::vector<size_t>& I,
		const double f) const noexcept {

	// cutoff radius	
	const double Rcut = (f>.0) ? radius()*f: -f;
	
	// find substitute matched positions
	std::vector<i_i> res; res.reserve(I.size());
	{
		// insert new entry lambda
		const auto ine = [&res](i_i ne) -> void {
			if ((ne.i1() == ne.i2()) ||
			    (ne.i1() == NPOS__)  ||
			    (ne.i2() == NPOS__)) return;
			const auto itr = std::lower_bound(res.cbegin(),res.cend(),ne);
			if (itr!=res.cend() && *itr == ne)
				return;
			res.insert(itr,ne);
		};

		auto i1 = I.cbegin();
		for (auto p1=Ap.ccBegin(),pe=Ap.ccEnd(); p1!=pe; ++p1,++i1) {
			auto i2 = I.cbegin();
			for (auto p2=Ap.ccBegin(); p2!=pe; ++p2,++i2) {
				const auto b = *p2-*p1;
				const i_i J{*i1,*i2};
				if (norm(b)>Rcut) continue;
				
				const auto H = getApproximateInteraction(J,b);
				if (!H.H.empty())
					ine(J.i1()<H.mi ? i_i{J.i1(),H.mi}: i_i{H.mi,J.i1()}),
					ine(J.i2()<H.pi ? i_i{J.i2(),H.pi}: i_i{H.pi,J.i2()});
			}
		}
	}

	return res;
}

// modification
size_t ll_hbonds::filter(const std::function<bool(
				const cMat& H, const fMat& B,
				const size_t i1, const size_t i2,
				const fMat& p1, const fMat& p2,
				const fCol& R)>& eval) noexcept {
	
	std::vector<size_t> rm; rm.reserve(this->size());
	size_t res=0;
	for (auto i=begin(); i<end(); ++i) {
		const fMat p1 = cell().Ap().cAt(i->i1()),
		      	   p2 = cell().Ap().cAt(i->i2());
		auto jh = i->cbegin();
		auto jr = i->ccBegin();
		
		// get indices to be removed
		std::vector<size_t> rm_; rm_.reserve(i->size());
		for (const auto e=i->cend(); jh!=e; ++jh,++jr)
			if (eval(*jh,cell().B(),i->i1(),i->i2(),p1,p2,*jr))
				rm_.push_back((size_t)jr);

		// remove in reverse order
		std::for_each(rm_.crbegin(),rm_.crend(),
			[i](const size_t ind)->void{i->erase(ind);});
		res+=rm_.size();

		if (i->empty()) rm.push_back(std::distance(begin(),i));
	}

	// remove empty index pairs
	std::for_each(rm.crbegin(),rm.crend(),
			[this](const size_t ind)->void{this->erase(ind);});

	// enforce symmetry
	if (!this->symmetric())
		res += filter([this](const cMat& H, const fMat& B,
			     const size_t i1, const size_t i2,
			     const fMat& p1, const fMat& p2,
			     const fCol& R) -> bool {
			const auto jtr = std::lower_bound(this->cbegin(),this->cend(),i_i(i2,i1));
			if (jtr==this->cend() || jtr->i1()!=i2 || jtr->i2()!=i1)
				return true;
			return !std::binary_search(jtr->ccBegin(),jtr->ccEnd(),-R);
		});
	
	return res;
}
void ll_hbonds::flipWannierSigns(const size_t ind, const std::vector<bool>& F) noexcept {
	assert(ind<cell().N());
	assert(F.size()==Norb(ind));

	for (auto& I: *this) {
		
		// flip from the left
		if (I.i1() == ind)
			for (auto& h: I.H()) {
				auto jf = F.cbegin();
				for (auto jh=h.rBegin(),eh=h.rEnd(); jh!=eh; ++jh,++jf)
					if (*jf) *jh = -(*jh);
			}
		
		// flip from the right
		if (I.i2() == ind)
			for (auto& h: I.H()) {
				auto jf = F.cbegin();
				for (auto jh=h.cBegin(),eh=h.cEnd(); jh!=eh; ++jh,++jf)
					if (*jf) *jh = -(*jh);
			}
	}
}


// io
void ll_hbonds::writeToFile(const std::string& fileName, const double Ef) const {
	// open file
	std::ofstream file;
	file.open(fileName, std::ios::binary);
	if (!file.good()) throw(std::invalid_argument("failed to open '"+fileName+"'"));

	// write header
	{
		char head[] = "wad90";
		file.write((char*) &head, 5*sizeof(char));
	}

	// write cell basis, positions and ids
	{
		// B, Ap
		uint32_t D = cell().dim(), N = cell().N();
		file.write((char*) &D, sizeof(uint32_t));
		file.write((char*) &N, sizeof(uint32_t));
		file.write((char*) cell().B().data(), D*D*sizeof(double));
		file.write((char*) cell().Ap().data(), D*N*sizeof(double));
		
		// ids
		uint8_t l;
		for (const auto& i: cell().id()) {
			l=i.size()+1;
			file.write((char*) &l, sizeof(uint8_t));
			file.write((char*) i.c_str(), l*sizeof(char));
		}
	}

	// write Norb
	{
		uint32_t tmp;
		for (const auto& i: Norb_) {
			tmp = i;
			file.write((char*) &tmp, sizeof(uint32_t));
		}
	}

	// write #pairs, pairs, #bonds, R and H
	{
		uint32_t tmp[3];
		
		// write #pairs
		tmp[0] = this->size();
		file.write((char*) tmp, sizeof(uint32_t));

		// write pairs, #bonds, R and H
		for (const auto& i: *this) {
			tmp[0]=i.i1(), tmp[1]=i.i2(), tmp[2]=i.size();
			file.write((char*) tmp, 3*sizeof(uint32_t));
	
			file.write((char*) i.R().data(), i.R().size()*sizeof(double));
			for (const auto& j: i)
				file.write((char*) j.data(), j.size()*sizeof(std::complex<double>));
		}
	}

	// write Ef
	file.write((char*) &Ef, sizeof(double));

	// close file
	file.close();
}

// searching
ll_hbonds::am_ ll_hbonds::getApproximateInteraction(const i_i& j, const fArray& b) const noexcept {
	
	// find interactions both ways
	cMat H12; size_t i2;
	for (i2=0; i2!=cell().N(); ++i2)
		if (Norb(i2)==Norb(j.i2())) {
			H12 = getInteraction({j.i1(),i2},b);
			if (!H12.empty()) break;
		}
	
	cMat H21; size_t i1;
	for (i1=0; i1!=cell().N(); ++i1)
		if (Norb(i1)==Norb(j.i1())) {
			H21 = T(getInteraction({j.i2(),i1},-b));
			if (!H21.empty()) break;
		}

	// return mean of interactions and actual target indices
	if (!H12.empty()) {
		if (!H21.empty()) return {(H12+H21)*.5,i2,i1};
		else		  return {H12*.5,i2,NPOS__};
	}
	if (!H21.empty())	  return {H21*.5,NPOS__,i1};
	else			  return {eH(),NPOS__,NPOS__};
}

// printing
std::ostream& ll_hbonds::print(std::ostream& os, const size_t mode) const noexcept {
	
	ll_bonds<ll__::i_i_R_H>::print(os,mode);

	switch (mode) {
	case 2:
	{
		// get envelope
		fMat lENV(cell().dim(),0); lENV.reserve(cardinality());
		fMat uENV(cell().dim(),0); uENV.reserve(cardinality());
		for (const auto& i: *this)
			lENV.push_back(nmin(i.R()).mat),
			uENV.push_back(nmax(i.R()).mat);
		
		// get widths
		const size_t lwidth =
			std::max_element(lENV.ccBegin(),lENV.ccEnd(),
			[](const auto& i, const auto& j)->bool{
				return T(i).print(0).size()<T(j).print(0).size();
				})->print(0).size()-1;
		const size_t uwidth =
			std::max_element(uENV.ccBegin(),uENV.ccEnd(),
			[](const auto& i, const auto& j)->bool{
				return T(i).print(0).size()<T(j).print(0).size();
				})->print(0).size()-1;

		// print cell
		os << cell() << "\n\n";

		// print bonds information
		for (const auto& i: *this)
			os << std::setw(11) << std::left << (ll__::i_i)i
			   << std::setw(21) << std::left <<
				(" size(H) = ["+std::to_string(Norb(i.i1()))+","
						+std::to_string(Norb(i.i2()))+"]")
			   << std::setw(13) << std::left <<
				(" NR = "+std::to_string(i.R().N()))
			   << "[" << std::setw(lwidth) << std::left <<
				nmin(i.R()).mat.T().print(0)
			   << "; " << std::setw(uwidth) << std::left <<
				nmax(i.R()).mat.T().print(0) << "]\n";
	}
	break;
	}

	return os;
}

// empty matrix
const cMat ll_hbonds::eH_;
