// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "ll_omen.h"
#include "ll_hbonds.h"
#include "ll_fn.h"
#include "ll_io.h"
#include "ll_hio.h"
#include "aux_io.h"
#include "aux_md5.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>


using namespace ll__;

// generate data for OMEN input
ll__::omen::prepper::prepper(const ll_wf_input& inp, std::ostream& os) {
	
	
	// write to disk helper lambda
	const auto wtd = [&inp,&os](const std::string& fileName, const auto& F,
				    const std::string& discr = "file") -> void {
		if (~inp.verbosity & PRINTBIT__) return;
		F.printToFile(fileName);
		if (inp.verbosity & PRINTBIT__)
			os << (discr.empty() ? "": discr+" ") << "'" << fileName << "' written"
		   	   << (inp.verbosity & MD5BIT__ ? ", "+aux::md5(fileName): "") << "\n";
	};


	// generate wbh
	W = ll_hbonds(inp,os);
	if (W.empty())
		throw(std::invalid_argument("generated empty wannier bonds hamiltonian"
				       	", check input parameters"));

	// find smallest orthorhombic cell
	{
		// check integrity of C, xyz
		if (!inp.C.empty()) {
			// we will use C, check integrity
			if (msize(inp.C)!=msize(W.cell().B()) || inp.C!=round(inp.C))
			throw(std::invalid_argument("C must be a square matrix of integers of size("
				+std::to_string(DIM__)+","+std::to_string(DIM__)+")"));
			
			// set C to input C
			C = inp.C;
		} else {
			// we will use xyz, check integrity
			if (msize(inp.xyz)!=msize(W.cell().B()) || !inp.xyz.ob())
			throw(std::invalid_argument("xyz must be a square orthogonal matrix of size("
				+std::to_string(DIM__)+","+std::to_string(DIM__)+")"));
			
			// find C using xyz as a template
			C = findBasisExpansion(W.cell().B(),inp.xyz,1000.0,inp.itol);
		}
		
		// find r for expansion C
		{
			const fMat ir = T(W.rmat());
			const auto rC = abs(ir.prod(C));
			if (!rC.logical() || sum(round(rC))!=sum(ir))
				throw(std::invalid_argument("bad expansion to orthorhombic cell, "
					"non vacuum conserving"));
			r.resize(ir.size());
			std::transform(rC.cbegin(),rC.cend(),r.begin(),
				[](const double i)->bool{return i;});
		}

		// enforce positive orientation
		auto P = eye<fMat>(msize(C));
		if (det(W.cell().B().prod(C))<.0) {
			const auto i = inds(r);
			if (i.size()) P(i.back(),i.back()) = -1.0;
			else P.back() = -1.0;
		}
		C = C.prod(P);

		fMat S = eye<fMat>(W.dim(),W.dim());
		// automatic stress tensor from xyz
		if (inp.C.empty() && inp.S.empty()) {
			
			// vanilla expanded basis
			const fMat B = W.cell().B().prod(C);
			
			// get Bhat, i.e. directions from xyz with length and sign of columns in B
			fMat Bhat(B.M(),0); Bhat.reserve(B.N());
			auto s = P.cdbegin();
			for (auto i=B.ccBegin(),e=B.ccEnd(), j=inp.xyz.ccBegin(); i!=e; ++i,++j,++s)
				Bhat.push_back((*j) * (*s)*(norm(*i)/norm(*j)));

			// scale non restricted subspace to conserve volume
			const auto Jn = ninds(r);
			const auto Jm = find(nany(inp.xyz.get({},Jn)));
			
			const auto Bhat_s = Bhat.get(Jm,Jn), B_s = B.get(Jm,Jn);
			const double f = std::pow(det(B_s)/det(Bhat_s),1.0/B_s.M());
			Bhat.set(Bhat_s*f,Jm,Jn);

			// stress matrix is S s.t. S*B = Bhat
			S = Bhat.prod(inv(B));
		}
		// user defined stress tensor
		if (!inp.S.empty()) {
			if (msize(inp.S)!=msize(W.cell().B()) || det(inp.S)<=.0)
				throw(std::invalid_argument("S must be a square matrix with det(S)>0"));
			S = inp.S;
		}

		// stress W
		W.stress(S);
		if (inp.verbosity & PRINTBIT__) {
			const double sn = norm(eye<fMat>(DIM__)-S);
			os << "strain applied, frobenius(1-S): "
			   << (sn<.05 ? GREEN__: RED__) << sn << RESET__ << "\n";
		}
		
		// stress wbh and get cell in new orthorhombic basis
		set_mtol(WTOL__);
		ORcell = W.cell().copy().expand(C);
		reset_mtol();

		// check resulting basis is orthorhombic
		if (!ORcell.B().ob())
			throw(std::invalid_argument("B*C not orthorhombic"));
	}
	if (inp.verbosity & PRINTBIT__)
		os << "\nfound minimal orthorhombic cell of " << BLUE__
		   << ORcell.N() << RESET__ << " positions\n" << "basis expansion coefficients:\n"
		   << T(C).print(0,1) << "\n\nrestricted dimensions are: " << r << "\n";
	if (inp.verbosity & WRITEBIT__)
		wtd(inp.prefix+"cell.psc",ORcell,"orthorhombic cell POSCAR file");


	// align ORcell and wbh with cartesian axis
	{
		const auto R = diag(mnorm(ORcell.B())).prod(inv(ORcell.B()));
		ORcell.rotate(R); W.rotate(R);
	}

	// find expansion to nn interaction only, filter wbh and expand ORcell
	size_t fcnt=0;
	{
		auto l = nmax(abs(W.range(ORcell.B()))).mat;
		for (size_t d=0; d!=DIM__; ++d)
			if (inp.l[d]) l[d]=inp.l[d];
			else if (l[d]==0.0) l[d]=1.0;
		NNE = diag(l);
		
		if (std::any_of(inp.filter_wbh.cbegin(),inp.filter_wbh.cend(),
				[](const bool i)->bool{return i;})) {
		fcnt = W.filter([fwb = inp.filter_wbh.cbegin(), TRM = inv(ORcell.B().prod(NNE))]
			(const size_t i1, const size_t i2,
			 const fMat& p1, const fMat& p2) -> bool {
				assert(i1||true||i2||p1.size());
				auto j = fwb; const fMat p2_ = TRM.prod(p2);
				for (const double i: p2_.lt(-1.0) | p2_.geq(2.0))
					if (*j++ && i) return true;
				return false;
			});
		}

		set_mtol(WTOL__);
		ORcell.expand(NNE);
		reset_mtol();

		set_mtol(1e-10);
		ORcell.autoShift(r);
		reset_mtol();
	}
	if (inp.verbosity & PRINTBIT__) {
		const auto Norb = W.Norb(W.cell().inds());
		const size_t Ntot = std::accumulate(Norb.cbegin(),Norb.cend(),0)
					* ORcell.N()/W.cell().N();
		
		os << "expansion to nn only is ("
		   << (fcnt ? RED__: GREEN__) << fcnt << RESET__
		   << "):\n" << T(NNE).print(0,1) << "\n"
		   << "\nresulting in " << BLUE__ << ORcell.N()
		   << RESET__ << " positions and "
		   << BLUE__ << Ntot << RESET__ " orbitals total\n"
		   << BLUE__ << fcnt << RESET__ << " bonds filtered\n";
	}
}


// write structural input
void ll__::omen::writeInput(const ll_wf_input& inp, std::ostream& os) {

	// generate the data
	const prepper P(inp,os);

	
	// get band edges
	const double Ef = std::isnan(inp.Ef) ? readEf(inp.outcar): inp.Ef;
	const auto E = findBandEdges(Ef,readEig(inp.weig));
	{
		std::ofstream file; file.open(inp.prefix+"Ef.dat");
		if (!file.good())
			throw(std::invalid_argument("open file \'"+inp.prefix+"Ef.dat"+"\' failed"));
		file << Ef;
	}
	if (inp.verbosity & MD5BIT__)
		os << "band edge information derived from files:\n"
		   << " '" << inp.weig << "'"
		   << (inp.verbosity & MD5BIT__ ? ", "+aux::md5(inp.weig): "") << "\n"
		   << (!std::isnan(inp.Ef) ? "": " '"+inp.outcar+"'"+
			(inp.verbosity & MD5BIT__ ? ", "+aux::md5(inp.outcar): "")+"\n");


	// write wbh to file
	P.W.writeToFile(WBH__,Ef);
	if (inp.verbosity & PRINTBIT__)
		os << "\nfile '" << WBH__ << "' written"
		   << (inp.verbosity & MD5BIT__ ? ", "+aux::md5(WBH__): "") << "\n";

	// run mesh BS tests
	{
		rv r; fMat EXP;
		switch(inp.expand_mesh) {
			case 0: EXP = eye<fMat>(P.W.dim(),P.W.dim()), r = P.W.r(); break;
			case 1: EXP = P.C, r = P.r; break;
			default: EXP = P.C.prod(P.NNE), r = P.r; break;
		}
		meshBStest(P.W,EXP,r,E,inp,os);
	}
	{
		rv r; fMat EXP;
		switch(inp.expand_trace) {
			case 0: EXP = eye<fMat>(P.W.dim(),P.W.dim()), r = P.W.r(); break;
			case 1: EXP = P.C, r = P.r; break;
			default: EXP = P.C.prod(P.NNE), r = P.r; break;
		}
		traceBStest(P.W,EXP,r,inp,os);
	}
	
	// get final structure and do cosmetics
	auto B_ = P.ORcell.B();
	auto Ap_ = P.ORcell.getcAp();
	auto id_ = P.ORcell.id(P.ORcell.type());
	{
		// shift restricted directions to origin
		for (const auto i: inds(P.r))
			Ap_.rAt(i) -= min(Ap_.rAt(i));

		// even out small rounding errors in cAp
		for (auto r=Ap_.rBegin(),re=Ap_.rEnd(); r!=re; ++r) {
			for (const auto i: *r) {
				// find average
				double s = i; size_t cnt=1;
				for (const auto j: *r)
					if (std::abs(i-j)<=WTOL__)
						s+=j, ++cnt;

				// set all equal coordinates to average
				s/=double(cnt);
				for (auto& j: *r)
					if (std::abs(i-j)<=WTOL__)
						j=s;
			}
		}
		
		// fix small deviance from 0.0 to 0.0
		for (auto r=Ap_.rBegin(),re=Ap_.rEnd(); r!=re; ++r) {
			const double m = min(*r);
			if (std::abs(m)<WTOL__) *r-=m;
		}
	}
	
	// get structural parameters
	double bl; size_t Nnn;
	{
		const auto b = P.W.cell().getSubCell(P.W.cell().fundamentalTypes()).
						getBonds(inp.bond_factor,NN);
		// bond length
		bl = b.radius()+WTOL__;
	
		// # next neigbours
		const auto Nb = b.Nindex(b.inds());
		Nnn = *std::max_element(Nb.cbegin(),Nb.cend())/2;
	}

	// adjust vacuum in B_
	for (const auto i: inds(P.r))
		B_.cAt(i) *= (max(Ap_.rAt(i))-min(Ap_.rAt(i))
			    +(inp.vac>=bl ? inp.vac: bl))
					/norm(B_.cAt(i));

	// sort Ap_, id_
	{
		const auto SO = aux::sorted_order(Ap_.cBegin(),Ap_.cEnd(),vcmp);
		aux::reorder(Ap_.cBegin(),SO);
		aux::reorder(id_.begin(),SO);
	}
	
	// get unique ids in order of occurrence
	idv uid_; uid_.reserve(P.ORcell.Nspecies());
	for (const auto& s: id_) {
		if (uid_.size()==uid_.capacity()) break;
		if (std::find(uid_.cbegin(),uid_.cend(),s)==uid_.cend())
			uid_.push_back(s);
	}


	// write lattice file
	printOlf(OLF__,B_,Ap_,id_,Nnn,bl);
	if (inp.verbosity & PRINTBIT__)
		os << "file '" << OLF__ << "' written"
		   << (inp.verbosity & MD5BIT__ ? ", "+aux::md5(OLF__): "") << "\n";
	
	// write material file
	printOmf(OMF__,
		 E,
		 P.W.Norb(uid_),
		 P.W.cell().mass(uid_));
	if (inp.verbosity & PRINTBIT__)
		os << "file '" << OMF__ << "' written"
		   << (inp.verbosity & MD5BIT__ ? ", "+aux::md5(OMF__): "") << "\n";

	// write input stump
	if (inp.device_length>.0) {
	switch (Nfalse(P.r)) {
		case 1:
			writeStump1D(B_,inp.device_length);
		break;
		case 2:
			writeStump2D(B_,inp.device_length);
		break;
		case 3:
			writeStump3D(B_,inp.device_length);
		break;
		default:
			throw(std::runtime_error("bad dimension for OMEN stump"));
		break;

		if (inp.verbosity & PRINTBIT__)
			os << "file '" << inp.prefix+OIS__ << "' written"
			   << (inp.verbosity & MD5BIT__ ? ", "+aux::md5(OIS__): "") << "\n";
	}
	}
}


// write hamiltonian matrices
void ll__::omen::hctor(const ll_wf_input& inp, const lm__::fMat& LM,
		const idv& id, const lm__::fMat& L, std::ostream& os) {
	assert(LM.M()==L.size());
	assert(LM.N()==id.size());

	// read wbhs, set tolerance
	ll_hbondss W(inp,os); W.fixDiagonals(1e-4); W.setQueryTolCartesian(inp.sptol);
	
	// convert ids to type indices
	const aTv T = W.ind(id);
	
	// adapt W to structure
	ll__::adaptWBH(W,inp,LM,T,os);

	// get interaction radius from inp or derived from wbh
	const double IR = std::isnan(inp.IR) ? W.radius(): inp.IR;

	// get R vector according to spacial restrictions
	const fMat R = rToR({true,L[1]>1e6,L[2]>1e6});

	// matrix block size
	const size_t Nw = std::accumulate(T.cbegin(),T.cend(),size_t(0),
		[&W](const size_t s, const aT t) -> size_t
			{ return s + W.Norb(t); });
	
	if (inp.verbosity & PRINTBIT__) {
		os << "\nlayer matrix of " << LM.N() << " positions loaded\n"
		   << "total number of orbitals in layer matrix: " << Nw << "\n"
		   << "number of R vectors is " << R.N() << "\n"
		   << "using interaction cutoff radius of " << IR << "\n\n";
	}

	// call hctor
	ll__::omen::writer writer(R,Nw);
	ll__::hctor(W,lm__::diag(L),LM,T,writer,IR,!inp.strict_matching,
			inp.Nthreads,inp.verbosity,os);

	// report md5's
	if (inp.verbosity & MD5BIT__) {
		os << "\nmd5 sums:";
		for (auto i=writer.R().ccBegin(),e=writer.R().ccEnd(); i!=e; ++i)
			os << "\n \'" << writer.fileName(*i) << "\', "
			   << aux::md5(writer.fileName(*i));
	}

	// check consistency
	if (!inp.force_check) return;

	os << "\n\ncheck forced!\n";
	for (auto i=writer.R().ccBegin(),e=writer.R().ccBegin()+R.N()/2+1; i!=e; ++i) {
		const std::string filep = writer.fileName( *i);
		const std::string filem = writer.fileName(-*i);
		
		if (inp.verbosity & PRINTBIT__)
			os << " reading file '" << filep << "' ...";
		const ll_sparse Hp(filep);

		ll_sparse Hm;
		if (filep!=filem) {
			if (inp.verbosity & PRINTBIT__)
				os << ", '" << filem << "' ...";
			Hm = ll_sparse(filem).T();
		} else  Hm = ll_sparse(Hp).T();

		if (inp.verbosity & PRINTBIT__)
			os << (Hp==Hm ? GREEN__ "   ok!": RED__ "   bad!") << RESET__ << "\n";
	}
}


// call material generator
void ll__::omen::callInputGenerator() {
	

	#ifndef NLOGFILE_
//	auto fs = aux::openFile<std::ofstream>(aux::timeStamp()+".log");
	auto fs = aux::openFile<std::ofstream>("winterface.log",std::ios::app);
	fs << "\n\n" << aux::timeStamp() << "\n\n";
	aux_tee os(std::cout,fs);
	#else
	std::ostream& os = std::cout;
	#endif


	// WINTERFACE
	if (std::ifstream("winput").good()) {
	
		// run interface
		hw_winterFace(os);
		try {
			writeInput(aux::parseFile<ll_wf_input>("winput",os),os);
		} catch(const std::exception& e) {
			os << "\n" << ll_wf_input() << "\n\n\n"
			   << RED__ << "ERROR: " << RESET__ << e.what() << "\n";
		}
		gb_winterFace(os);
	}
	

	// PHONON INTERFACE
	if (std::ifstream("phinput").good()) {
	
		// run interface
		try {
			hw_phinterFace(os);
			const auto inp = aux::parseFile<ll_ph_input>("phinput",os);
			
			/*
			 * your stuff here
			 */
			
			gb_phinterFace(os);
		} catch(const std::exception& e) {
			os << "\n" << ll_ph_input() << "\n\n\n"
			   << RED__ << "ERROR: " << RESET__ << e.what() << "\n";
		}
	}
}


// call scaler
void ll__::omen::callHamiltonianConstructor(double* const Layer_Matrix, const idv& id,
					const double Lx, const double Ly, const double Lz) {
	
#ifndef NLOGFILE_
//	auto fs = aux::openFile<std::ofstream>(aux::timeStamp()+".log");
	auto fs = aux::openFile<std::ofstream>("winterface.log",std::ios::app);
	fs << aux::timeStamp() << "\n\n";
	aux_tee os(std::cout,fs);
#else
	std::ostream& os = std::cout;
#endif
	
	try {
		// check if input file exists
		{
			std::ifstream file;
			file.open("winput");
			if (!file.good()) return;
		}

		// check timestamps
		{
			// lattice_dat
			struct stat st_olf;
			if(stat(OLF__, &st_olf)) return;
			
			// mat_par
			struct stat st_omf;
			if(stat(OMF__, &st_omf)) return;
			
			// H_4.bin
			struct stat st_h4;
			int ierr = stat("H_4.bin", &st_h4);

			// do nothing if H_4.bin is younger than
			// lattice_dat and mat_par
			if (!ierr &&
			st_h4.st_mtime>st_olf.st_mtime &&
			st_h4.st_mtime>st_omf.st_mtime) return;
		}
		
		// run interface
		{
			hw_winterFace(os);

			hctor(aux::parseFile<ll_wf_input>("winput",os),
				fMat(Layer_Matrix,3,id.size()),id,{Lx,Ly,Lz},os);

			gb_winterFace(os);
		}
	} catch(const std::exception& e) {
		os << "\n" << ll_wf_input() << "\n\n\n"
		   << RED__ << "ERROR: " << RESET__ << e.what() << "\n";
	}
}


// stumps
void ll__::omen::writeStump1D(const lm__::fMat& B, const double L, const double ftc) {
	// add code here
	std::cout << "\n\n1D stump not implemented, sorry\n\n";
}
void ll__::omen::writeStump2D(const lm__::fMat& B, const double L, const double ftc) {
	using namespace lm__;
	
	assert(B.square());
	assert(B.M()==DIM__);
	assert(std::abs(det(B))>mtol());
	assert(B==diag(diag(B)));

	// open file
	std::ofstream file;
	file.open(OIS__);
	if (!file.good())
		throw(std::invalid_argument("failed to open '"+std::string(OIS__)+"'"));
	
	// write upto structure section
	file << "/*Parameters*/\n";
	file << "mat_name = " << OMF__ << ";\n";
	file << "lattice_type = " << OLF__ << ";\n";
	file << "a0 = 0.1;\n";
	file << "first_atom = anion;\n";
	file << "poisson_solver	= 1;\n";
	file << "poisson_iteration = 15;\n";
	file << "charge_transfer = 0;\n";
	file << "read_hamiltonian = 1;\n";
	file << "NzFold	= 1;\n";
	file << "NDim = 2;\n";
	file << "Nkz = 1;\n";
	file << "tb = 10;\n";
	file << "dsp3 = 0;\n";
	file << "Temp = 300;\n";
	file << "n_of_modes = 8;\n";
	file << "Nk = 251;\n";
	file << "bs_solver = full;\n";
	file << "max_bond_def = 0.1;\n";
	file << "last_first = 0;\n";
	file << "x = [1 0 0];\n";
	file << "y = [0 1 0];\n";
	file << "z = [0 0 1];\n";
	file << "eta_res = 0;\n";
	file << "eta = 0;\n";
	file << "Elimit = 50e-3;\n";
	file << "Emin_tail = 5.0e-3;\n";
	file << "EOffset = 10*UT;\n";
	file << "dE_in = 5.0e-5;\n";
	file << "dE_f = 5.0e-5;\n";
	file << "dE_sep	= 2.0e-5;\n";
	file << "NEmax = 250;\n";
	file << "CPU_per_kz_point = all;\n";
	file << "CPU_ppoint = 1;\n";
	file << "Eps_wire = 6.0;\n";
	file << "Eps_ox = 20.0;\n";
	file << "Xi_wire = 4.05;\n";
	file << "phi_m = 4.05;\n";
	file << "NVG = [16 2];\n";
	file << "Vgmin = [-0.35 0.35];\n";
	file << "Vgmax = [0.15 0.35001];\n";
	file << "NVS = 1;\n";
	file << "Vsmin = 0.0;\n";
	file << "Vsmax = 0.0;\n";
	file << "NVD = 1;\n";
	file << "Vdmin = 0.4;\n";
	file << "Vdmax = 0.4;\n";
	file << "update_atoms = 0;\n";
	file << "atom_file = atom_pos_dat;\n";

	file << "\n/*Structure*/\n";
	file << "grid_accuracy = 1;\n";
	file << "no_mat	= 2;\n";
	file << "no_channel_mat = 1;\n";
	file << "no_oxide_mat = 1;\n";

	// write structure info from B and L
	const size_t f = std::ceil(L/B(0,0));
	const double l = (double)f*B(0,0)/30.0;
	file << "\n// number of unit cells: " << f << "\n";
	file << "Lc = " << std::setprecision(12) << l << ";\n";
	file << "Ls = " << std::setprecision(12) << l << ";\n";
	file << "Ld = " << std::setprecision(12) << l << ";\n";
	file << "tc = " << std::setprecision(12) << (B(1,1)*ftc*0.1) << ";\n\n";

	// write tail
	file << "tox = 2.9;\n";
	file << "hox = 3.0;\n";
	file << "mat_type(1) = square;\n";
	file << "mat_id(1) = 1;\n";
	file << "mat_cs(1) = yes;\n";
	file << "mat_coord(1,1)	= [0.0 0.0];\n";
	file << "mat_coord(1,2) = [Ls+Lc+Ld 0.0];\n";
	file << "mat_coord(1,3) = [Ls+Lc+Ld tc];\n";
	file << "mat_coord(1,4)	= [0.0 tc];\n";
	file << "ox_type(1) = square;\n";
	file << "ox_id(1) = 1;\n";
	file << "ox_cs(1) = no;\n";
	file << "ox_coord(1,1) = [0.0 0.0-hox];\n";
	file << "ox_coord(1,2) = [Ls+Lc+Ld 0.0-hox];\n";
	file << "ox_coord(1,3) = [Ls+Lc+Ld tc+tox];\n";
	file << "ox_coord(1,4) = [0.0 tc+tox];\n";
	file << "no_gate = 2;\n";
	file << "gate_type(1) = square;\n";
	file << "gate_coord(1,1) = [0.0 tc+tox];\n";
	file << "gate_coord(1,2) = [Ls+Lc+Ld tc+tox];\n";
	file << "gate_type(2) = square;\n";
	file << "gate_coord(2,1) = [0.0 0.0-hox];\n";
	file << "gate_coord(2,2) = [Ls+Lc+Ld 0.0-hox];\n";

	file << "\n/*Commands*/\n";
	file << "command(1) = Write_Layer_Matrix;\n";
	file << "command(2) = Write_Fermi_Level;\n";
	file << "command(3) = CB_Bandstructure;\n";
	file << "command(4) = VB_Bandstructure;\n";
	file << "command(5) = Write_Grid_Matrix;\n";

	file.close();
}
void ll__::omen::writeStump3D(const lm__::fMat& B, const double L, const double ftc) {
	// add code here
	std::cout << "\n\n3D stump not implemented, sorry\n\n";
}
