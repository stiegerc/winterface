// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "b_switch.h"
#include "ll_io.h"
#include "ll_fn.h"
#include "ll_sparse.h"
#include "ll_omen.h"

using namespace lm__;
using namespace ll__;

void b_switch(const b_input& inp, std::ostream& os) {
try {
	// check input properties
	if (inp.C!=round(inp.C) || std::abs(det(inp.C))<.5)
		throw(std::invalid_argument("parameter C invalid"));
	if (inp.R!=round(inp.R) || inp.R.M()!=DIM__)
		throw(std::invalid_argument("parameter R invalid"));
	if (inp.Nk.size()!=1 && (inp.Nk.size()!=inp.kpts.N()-1))
		throw(std::invalid_argument("bad Nk"));
	if (inp.rho_k<0.0)
		throw(std::invalid_argument("bad rho_k"));
	if (!inp.k.empty() && inp.k.M()!=DIM__)
		throw(std::invalid_argument("bad k"));
	

	// lambdas to write unit cell as POSCAR file
	const auto wPSC = [&inp,&os](const std::string& fileName,
			const fMat& B, const fMat& Ap, const idv& id,
			const bool direct) -> void {
		if (~inp.verbosity & WRITEBIT__) return;

		printPOSCAR(fileName,B,Ap,id,1.0,direct);
		if (inp.verbosity & PRINTBIT__)
			os << "file '" << fileName << "' written\n";
	};
	const auto wPSCc = [&wPSC](const std::string& fileName, const ll_cell& cell) -> void {
		wPSC(fileName,cell.B(),cell.Ap(),cell.id(cell.type()),true);
	};

	// lambda to generate reciprocal basis
	const auto getRB = [](const fMat& B)->fMat{return 2.0*M_PI*T(inv(B));};

	// lambda to generate path using basis B
	const auto gnwPath = [&inp,&os,&getRB](const fMat& B) -> p_p {
		p_p res;
		if (inp.Nk.size()==1) {
			const fMat RB = getRB(B);
			res = genPath(inp.kpts,inp.Nk.front(),RB);
			if (inp.verbosity & PRINTBIT__)
				os << "\ngenerated path of " << BLUE__ << res.N()
				   << RESET__ << " kpoints using\n"
				   << "B:\n" << T(B).print(9,1)
				   << "\nRB:\n" << T(RB).print(9,1) << "\n\n";
		} else {
			res = genPath(inp.kpts,inp.Nk);
			if (inp.verbosity & PRINTBIT__)
				os << "\ngenerated path of " << BLUE__ << res.N()
				   << RESET__ << "kpoints using Nk: [" << inp.Nk << "]\n";
		}

		if (inp.verbosity & WRITEBIT__) {
			res.pos().writeToFile(inp.prefix+"pos.bin");
			if (inp.verbosity & PRINTBIT__)
				os << "file '" << inp.prefix << "pos.bin' written\n";
			res.path().writeToFile(inp.prefix+"path.bin");
			if (inp.verbosity & PRINTBIT__)
				os << "file '" << inp.prefix << "path.bin' written\n";
		}
		
		return res;
	};
	// lambda to generate mesh using basis B and restriction vector r
	const auto gnwMesh = [&inp,&os,&getRB](const fMat& B, const rv& r) -> ll_mesh<> {
		// check bzbounds
		if (inp.bzbounds.size()!=2)
			throw(std::invalid_argument("bad bzbounds"));
		if (inp.bzbounds[0]>=inp.bzbounds[1])
			throw(std::invalid_argument("bad bzbounds"));

		std::vector<size_t> maj;
		switch (fnvHash(inp.maj_style.c_str())) {
			case "inorder"_h: maj = maj_default(r.size()); break;
			case "matlab"_h: maj = maj_MATLAB(r); break;
			case "MATLAB"_h: maj = maj_MATLAB(r); break;
			default:
				throw(std::invalid_argument("maj style '"+
					inp.maj_style+"'not recognized"));
		}

		const fMat RB = getRB(B);
		const auto res = genMesh_cell(inp.rho_k,RB,inp.bzbounds[0],inp.bzbounds[1],maj,r);

		if (inp.verbosity & PRINTBIT__) {

			// fractions of column in RB vs number of mesh points
			std::vector<double> l; l.reserve(Nfalse(r));
			for (const auto i: ninds(r))
				l.push_back(norm(RB.cAt(i))/double(res.D(i)));

			os << "\ngenerated mesh of dimensions: [" << BLUE__ << res.D()
			   << RESET__ << "], majority: [" << BLUE__ << res.maj()
			   << RESET__ << "]\nrho_k: "
			   << BLUE__ << inp.rho_k << RESET__
			   << ", r: [" << BLUE__ << r << RESET__ << "]"
			   << ", l: [" << BLUE__ << std::fixed << std::setprecision(5)
			   << l << RESET__ << "]"
			   << "\nB:\n" << T(B).print(9,1)
			   << "\nRB:\n" << T(RB).print(9,1) << "\n\n";
		}

		if (inp.verbosity & WRITEBIT__) {
			res.writeToFile(inp.prefix+"mesh.bin");
			if (inp.verbosity & PRINTBIT__)
				os << "file '" << inp.prefix << "mesh.bin' written\n";
		}

		return res;
	};


	// lambda to write EGC to disk
	const auto wEGC = [inp,&os](const auto& EGC, const std::string& idstr) -> void {
		if (!(inp.verbosity & WRITEBIT__)) return;

		{
			const std::string fileName = inp.prefix+"E"+idstr+".bin";
			EGC.E.writeToFile(fileName);
			if (inp.verbosity & PRINTBIT__)
				os << "file '" << fileName << "' written\n";
		}
		{
			const std::string fileName = inp.prefix+"G"+idstr+".bin";
			EGC.G.writeToFile(fileName);
			if (inp.verbosity & PRINTBIT__)
				os << "file '" << fileName << "' written\n";
		}
		{
			const std::string fileName = inp.prefix+"C"+idstr+".bin";
			EGC.C.writeToFile(fileName);
			if (inp.verbosity & PRINTBIT__)
				os << "file '" << fileName << "' written\n";
		}
	};


	// mode switch
	switch (fnvHash(inp.mode.c_str())) {
	case "fold"_h:
	{
		if (inp.verbosity & PRINTBIT__)
			os << "computing bandstructure using '"
			   << GREEN__ << "folding" << RESET__ << "' mode\n";

		// read data from files
		const auto B = readB(inp.wout);
		const auto hr = readHr(inp.hrdat);
		if (inp.verbosity & PRINTBIT__)
			os << "input data from files '" << inp.wout
			   << "' and '" << inp.hrdat << "'\n";

		const auto NB = B.prod(inp.C);
		if (inp.verbosity & PRINTBIT__)
			os << "expanding using matrix\n" << inp.C.print(0,1) << "\n";
		if (inp.verbosity & WRITEBIT__)
			NB.printToFile(inp.prefix+"NB.mat");

		// trace
		if (inp.Nk.front())
			wEGC(calcFoldedBS_gc(hr,
				gnwPath(NB).path(),
				NB,B,inp.Nthreads),"fold");
		// mesh
		if (inp.rho_k) {
			// get correct r for matrix C
			const auto r_ = T(fMat(inp.r.empty() ? hr.r(): inp.r)).prod(inp.C);
			rv r(r_.size());
			std::transform(r_.cbegin(),r_.cend(),r.begin(),
				[](const double i)->bool{return i;});

			wEGC(calcFoldedBS_gc(hr,
				gnwMesh(NB,r),
				NB,B,inp.Nthreads),"fold_mesh");
		}
		// k spec
		if (!inp.k.empty())
			wEGC(calcFoldedBS_gc(hr,inp.k,NB,B,inp.Nthreads),"fold_spec");
	}
	break;
	case "scale"_h:
	{
		if (inp.verbosity & PRINTBIT__)
			os << "computing bandstructure using '"
			   << GREEN__ << "scaling" << RESET__ << "' mode\n";
		
		// read wbh from file
		const ll_hbondss W(inp,os); W.setQueryTolCartesian(inp.sptol);
		
		// read structure
		ll__::olf OLF;
		{
			if (inp.lattice_dat.empty()) {
				if (W.size()>1)
					throw(std::invalid_argument("found more than one whb"
						", need lattice_dat file"));

				// use cell from wbh and expand with C matrix
				if (inp.verbosity & PRINTBIT__)
					os << "using cell from file wannier bonds file\n"
					   << "expanding using matrix\n" << inp.C.print(0,1) << "\n";
				const auto cell = W.front().cell().copy().expand(inp.C);
				if (inp.verbosity & WRITEBIT__)
					cell.B().printToFile(inp.prefix+"NB.mat");
				
				OLF.B = cell.B();
				OLF.Ap = cell.getcAp();
				OLF.T = cell.type();
				OLF.id = cell.id(OLF.T);
				OLF.nn = 0;
				OLF.bl = .0;

			} else {
				if (inp.verbosity & PRINTBIT__)
					os << "using structure from file '" << inp.lattice_dat << "'\n";
				OLF = readOlf(inp.lattice_dat);
			}
		}

		// restriction vector
		const rv r = inp.r.empty() ? W.r(): inp.r;

		if (inp.re) {
			// generate hamiltonians
			set_mtol(inp.sptol);
			const auto HH = !inp.R.empty() ?
				genHam<fMat>(OLF.B,OLF.Ap,OLF.id,W,inp.R,inp.strict_matching):
				genHam<fMat>(OLF.B,OLF.Ap,OLF.id,W.r(),W,inp.strict_matching);
			reset_mtol();
			if (inp.dump_hamiltonians) HH.writeToFile(inp.prefix+"scal");

			// trace
			if (inp.Nk.front())
				wEGC(calcBS_gc(HH,gnwPath(OLF.B).path(),
						OLF.B,inp.Nthreads),"scal");
			// mesh
			if (inp.rho_k)
				wEGC(calcBS_gc(HH,gnwMesh(OLF.B,r),
						OLF.B,inp.Nthreads),"scal_mesh");
			// spec
			if (!inp.k.empty()) 
				wEGC(calcBS_gc(HH,inp.k,OLF.B,inp.Nthreads),"scal_spec");
		} else {
			// generate hamiltonians
			const auto HH = !inp.R.empty() ?
				genHam<cMat>(OLF.B,OLF.Ap,OLF.id,W,inp.R,inp.strict_matching):
				genHam<cMat>(OLF.B,OLF.Ap,OLF.id,r,W,inp.strict_matching);
			if (inp.dump_hamiltonians) HH.writeToFile(inp.prefix+"scal");

			// trace
			if (inp.Nk.front())
				wEGC(calcBS_gc(HH,gnwPath(OLF.B).path(),
						OLF.B,inp.Nthreads),"scal");
			// mesh
			if (inp.rho_k)
				wEGC(calcBS_gc(HH,gnwMesh(OLF.B,r),
						OLF.B,inp.Nthreads),"scal_mesh");
			// spec
			if (!inp.k.empty()) 
				wEGC(calcBS_gc(HH,inp.k,OLF.B,inp.Nthreads),"scal_spec");
		}
	}
	break;
	case "legacy"_h:
	{
		if (inp.verbosity & PRINTBIT__)
			os << "computing bandstructure using '"
			   << GREEN__ << "legacy" << RESET__ << "' mode\n";
		
		// read data from files
		const auto OMF = readOmf(inp.inprefix+inp.mat_par);
		const auto OLF = readOlf(inp.inprefix+inp.lattice_dat);
		if (inp.verbosity & PRINTBIT__)
			os << "input data from files '" << inp.mat_par
			   << "' and '" << inp.lattice_dat << "'\n";

		// Norb for structure in L
		const size_t Norb = std::accumulate(OLF.T.cbegin(),OLF.T.cend(),size_t(0),
			[&OMF](const size_t s, const aT t)->size_t{return s+OMF.Norb[t];});

		// get hr from subblocks of sparse device hamiltonians
		R_H<cMat> hr;
		{
			const auto dh = omen::readOMENhr(inp.inprefix,inp.r.empty() ?
				rv{true,true,false}: rv{true,inp.r[1],inp.r[2]});

			fMat R(DIM__,0); R.reserve(dh.size()*DIM__);
			std::vector<cMat> H; H.reserve(dh.size()*DIM__);

			auto r = dh.ccBegin();
			for (auto i=dh.cbegin(),e=dh.cend(); i!=e; ++i,++r) {
				
				R.push_back(*r); R.cBack()[0] = -1.0;
				H.push_back(i->get(Norb,0,Norb,Norb).full());
				
				R.push_back(*r); R.cBack()[0] =  0.0;
				H.push_back(i->get(0,0,Norb,Norb).full());
				
				R.push_back(*r); R.cBack()[0] =  1.0;
				H.push_back(i->get(0,Norb,Norb,Norb).full());
			}

			hr = {std::move(R),std::move(H)};
		}
		if (inp.dump_hamiltonians) hr.writeToFile(inp.prefix+"scal");
	
		// trace
		if (inp.Nk.front())
			wEGC(calcBS_gc(hr,
				gnwPath(OLF.B).path(),
				OLF.B,inp.Nthreads),"scal");

		// mesh
		if (inp.rho_k)
			wEGC(calcBS_gc(hr,
				gnwMesh(OLF.B,inp.r.empty() ? rv{false,true,false}: inp.r),
				OLF.B,inp.Nthreads),"scal_mesh");

		// spec
		if(!inp.k.empty())
			wEGC(calcBS_gc(hr,inp.k,OLF.B,inp.Nthreads),"scal_spec");
	}
	break;
	case "hr"_h:
	{
		if (inp.verbosity & PRINTBIT__)
			os << "computing bandstructure using '"
			   << GREEN__ << "hr" << RESET__ << "' mode\n";
		
		// read layer matrix
		const auto OLF = readOlf(inp.layer_matrix);

		if (!inp.re) { // keep complex parts
			// read hr and convert to full
			R_H<cMat> hr;
			{
				const auto hr_ = readHrSparse(inp.hr);
				std::vector<cMat> H; H.reserve(hr.size());
				for (auto i=hr_.cbegin(),e=hr_.cend(); i!=e; ++i)
					H.push_back(i->full());
				hr = {hr_.R(),std::move(H)};
			}


			// trace
			if (inp.Nk.front())
				wEGC(calcBS_gc(hr,
					gnwPath(OLF.B).path(),
					OLF.B,inp.Nthreads),"scal");

			// mesh
			if (inp.rho_k)
				wEGC(calcBS_gc(hr,
					gnwMesh(OLF.B,inp.r.empty() ? rv{false,true,false}: inp.r),
					OLF.B,inp.Nthreads),"scal_mesh");

			// spec
			if(!inp.k.empty())
				wEGC(calcBS_gc(hr,inp.k,OLF.B,inp.Nthreads),"scal_spec");
		
		} else { // ignore complex parts
			// read hr and convert to full
			R_H<fMat> hr;
			{
				const auto hr_ = readHrSparse(inp.hr);
				std::vector<fMat> H; H.reserve(hr.size());
				for (auto i=hr_.cbegin(),e=hr_.cend(); i!=e; ++i)
					H.push_back(real(i->full()));
				hr = {hr_.R(),std::move(H)};
			}


			// trace
			if (inp.Nk.front())
				wEGC(calcBS_gc(hr,
					gnwPath(OLF.B).path(),
					OLF.B,inp.Nthreads),"scal");

			// mesh
			if (inp.rho_k)
				wEGC(calcBS_gc(hr,
					gnwMesh(OLF.B,inp.r.empty() ? rv{false,true,false}: inp.r),
					OLF.B,inp.Nthreads),"scal_mesh");

			// spec
			if(!inp.k.empty())
				wEGC(calcBS_gc(hr,inp.k,OLF.B,inp.Nthreads),"scal_spec");
		}
	}
	break;
	case "local"_h:
	{
		// read wbh, get sym cell
		const ll_hbonds W(inp.wbh.front()); W.setQueryTolCartesian(inp.sptol);
		const ll_cell scell = W.symCell();

		const rv r = inp.r.empty() ? W.r(): inp.r;
		if (inp.verbosity & PRINTBIT__)
			os << "found similarized cell of " << BLUE__
			   << scell.N() << RESET__ << " positions, "
			   << BLUE__ << scell.Nspecies() << RESET__ << " types\n";
		if (inp.verbosity & DEBUGBIT__)
			wPSCc(inp.prefix+"scell.psc",scell);

		// find primitive version of sym cell
		auto pcell = scell;
		if (!inp.LB.empty()) {
			// user defined basis
			set_mtol(pcell.directTol(inp.sptol));
			if (std::abs(det(inp.LB))>pcell.vol() || !pcell.validBasis(inp.LB)) {
				reset_mtol();
				throw(std::invalid_argument("specified basis invalid or bad tolerance level"));
			}
			pcell.changeBasis(inp.LB);
			reset_mtol();
		} else {
			// find primitive cell automatically
			set_mtol(pcell.directTol(inp.sptol));
			pcell.makePrimitive(r);
			reset_mtol();
		}

		if (inp.verbosity & PRINTBIT__)
			os << "found primitive similarized cell of " << BLUE__ << pcell.N()
			   << RESET__ << " positions\n";
		if (inp.verbosity & DEBUGBIT__)
			wPSCc(inp.prefix+"pcell.psc",pcell);

		// find expansion R vectors and sort them
		fMat R;
		{
			ll_cell rcell(pcell.B(),zeros<fMat>(scell.dim(),1),{"A"});
			set_mtol(rcell.directTol(inp.sptol));
			rcell.changeBasis(scell.B());
			reset_mtol();

			R = round(pcell.B().leftDivide(rcell.getcAp()));
			std::sort(R.cBegin(),R.cEnd());
		}
		if (inp.verbosity & DEBUGBIT__) {
			const std::string fileName = inp.prefix+"Rvecs.mat";
			R.printToFile(fileName,0);
			if (inp.verbosity & PRINTBIT__)
				os << "file '" << fileName << "' written\n";
		}

		// transformation matrix pcell -> scell, possible shiftvectors
		const fMat TM = scell.B().leftDivide(pcell.B());

		// generate mesh and trace
		const ll_mesh<> mesh = inp.rho_k ? gnwMesh(pcell.B(),inp.r.empty() ? W.r(): inp.r):
						   ll_mesh<>();
		const auto PP = inp.Nk.front() ? gnwPath(pcell.B()): ll__::p_p();
			
		// find cells for each R vector
		for (auto r=R.ccBegin(),re=R.ccEnd(); r!=re; ++r) {

			// get shifted positions of pcell in the basis of scell
			fMat shm = pcell.Ap(); cadd(shm,*r); shm = (TM.prod(shm)%1.0);

			// get id strings and positions in cartesian
			set_mtol(inp.sptol);
			idv id; id.reserve(shm.N());
			fMat cAp(W.cell().dim(),0); cAp.reserve(shm.N());
			for (auto s=shm.ccBegin(),se=shm.ccEnd(); s!=se; ++s) {
				const auto itr = std::find(W.cell().ccBegin(),W.cell().ccEnd(),*s);

				if (itr==W.cell().ccEnd()) {
					reset_mtol();
					throw(std::invalid_argument("position not found! Increase tolerance"));
				}

				id.push_back(W.cell().id((size_t)itr));
				cAp.push_back(W.cell().B().prod(*itr));
			}
			reset_mtol();
			if (inp.verbosity & PRINTBIT__)
				os << "\nR: " << lm__::T(*r).print(0) << ",  ID: " << id << "\n";

			// type vector
			aTv T(id.size()); std::iota(T.begin(),T.end(),0);

			// idstring for filenames
			std::string idstr;
			{
				std::stringstream sstr;
				sstr << "local_" << std::setw(std::ceil(std::log10(R.N()))+1)
				     << std::setfill('0') << ((size_t)r+1);
				idstr = sstr.str();
			}

			// print POSCAR file
			if (inp.verbosity & WRITEBIT__) {
				const std::string fileName = inp.prefix+idstr+".psc";
				printPOSCAR(fileName,pcell.B(),cAp,id,1.0,false);
				if (inp.verbosity & PRINTBIT__)
					os << "file '" << fileName << "' written\n";
			}

			// generate Hamiltonian and compute bandstructure
			if (inp.re) {
				const auto HH = genHam<fMat>(pcell.B(),cAp,id,W.r(),W,false);
				if (inp.dump_hamiltonians) HH.writeToFile(inp.prefix+idstr);
				if (inp.Nk.front())
					wEGC(calcBS_gc(HH,PP.path(),pcell.B(),inp.Nthreads),"_trace_"+idstr);
				if (inp.rho_k)
					wEGC(calcBS_gc(HH,mesh,pcell.B(),inp.Nthreads),"_mesh_"+idstr);
			} else {
				const auto HH = genHam<cMat>(pcell.B(),cAp,id,W.r(),W,false);
				if (inp.dump_hamiltonians) HH.writeToFile(inp.prefix+idstr);
				if (inp.Nk.front())
					wEGC(calcBS_gc(HH,PP.path(),pcell.B(),inp.Nthreads),"_trace_"+idstr);
				if (inp.rho_k)
					wEGC(calcBS_gc(HH,mesh,pcell.B(),inp.Nthreads),"_mesh_"+idstr);
			}
		}
	}
	break;
	default:
		throw(std::invalid_argument("mode '"+inp.mode+"' unknown"));
	}
} catch(const std::exception& e) {
	os << b_input() << "\n\n" << RED__
	   << "ERROR: " << e.what() << "\n" << RESET__;
}
}
