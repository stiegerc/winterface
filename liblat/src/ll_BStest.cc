// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "ll_BStest.h"
#include "ll_fn.h"
#include "ll_io.h"
#include "ll_cell.h"
#include "aux_md5.h"


// bandstructure tests
void ll__::meshBStest(const ll_hbonds& W, const fMat& EXP, const rv& r,
		const vb_cb& E, const ll_BStest_input& inp, std::ostream& os) {
	if (!(inp.rho_k && inp.verbosity)) return;
	if (inp.rho_k<.0) throw(std::invalid_argument("negative rho_k not allowed"));
	if (inp.rho_k>1e6) throw(std::invalid_argument("rho_k > 1e6, don't be silly!"));
	if (EXP.empty() || std::abs(det(EXP))<mtol() || !(EXP == round(EXP)))
		throw(std::invalid_argument("bad expansion matrix"));


	// write to disk helper lambda
	const auto wtd = [&inp,&os](const std::string& fileName, const auto& F,
				    const std::string& discr = "file") -> void {
		if (~inp.verbosity & PRINTBIT__) return;
		F.writeToFile(fileName);
		if (inp.verbosity & PRINTBIT__)
			os << (discr.empty() ? "": discr+" ") << "'" << fileName << "' written"
		   	   << (inp.verbosity & MD5BIT__ ? ", "+aux::md5(fileName): "") << "\n";
	};

	// generate test cell
	set_mtol(WTOL__); // set tolerance due to wannier input precision
	const auto tCell = W.cell().copy().expand(EXP);
	const auto RB = tCell.getRB();
	reset_mtol();
	
	if (inp.verbosity & PRINTBIT__)
		os << "\n" << CYAN__ << "MESH BANDSTRUCTURE TEST" << RESET__
		   << "\nexpansion C:\n" << EXP.print(0,1)
		   << "\nreciprocal basis RB:\n"
		   << RB.print(PPREC__,1) << "\n\n";
	if (inp.verbosity & DEBUGBIT__)
		wtd(inp.prefix+"mesh_cell.psc",tCell);

	// read wannier90 hamiltonian
	const auto hr = readHr(inp.hrdat);

	// generate kmesh
	std::vector<size_t> maj;
	switch(aux::fnvHash(inp.maj_style.c_str())) {
		case "inorder"_h: maj = maj_default(r.size()); break;
		case "matlab"_h: maj = maj_MATLAB(r); break;
		case "MATLAB"_h: maj = maj_MATLAB(r); break;
		default:
			throw(std::invalid_argument(
			"maj style '"+inp.maj_style+"' not recognized"));
	}
	const auto kmesh = genMesh_cell(inp.rho_k,RB,
				inp.bzbounds[0],inp.bzbounds[1],maj,r);
	
	if (inp.verbosity & PRINTBIT__)
		os << "\nrunning test on k-mesh of dimensions ["
		   << BLUE__ << kmesh.D() << RESET__ << "], majority: ["
		   << BLUE__ << kmesh.maj() << RESET__ << "]"
		   << ", rho_k: " << BLUE__ << inp.rho_k << RESET__ << " ...\n\n";

	// get scaled and folded energies
	const auto EGCscal = inp.re ?
		calcBS_gc(genHam<fMat>(tCell,r,W,true),kmesh,tCell.B(),inp.Nthreads):
		calcBS_gc(genHam<cMat>(tCell,r,W,true),kmesh,tCell.B(),inp.Nthreads);
	set_mtol(WTOL__);
	const auto EGCfold = calcFoldedBS_gc(hr,kmesh,tCell.B(),W.cell().B(),inp.Nthreads);
	reset_mtol();

	// first conduction band index, bounds
	const size_t cbi = distance(EGCfold.E.base_.cFront().cbegin(),
		   std::lower_bound(EGCfold.E.base_.cFront().cbegin(),
				    EGCfold.E.base_.cFront().cend(),.5*(E.vb+E.cb)));
	const bool have_vbands=cbi>0, have_cbands=cbi<EGCfold.Nb();
	const double vb_bnd=E.vb-inp.Lvb, cb_bnd=E.cb+inp.Lcb;

	// write files
	wtd(inp.prefix+"mesh.bin",kmesh);
	wtd(inp.prefix+"Efold_mesh.bin",EGCfold.E);
	wtd(inp.prefix+"Escal_mesh.bin",EGCscal.E);
	wtd(inp.prefix+"Gfold_mesh.bin",EGCfold.G);
	wtd(inp.prefix+"Gscal_mesh.bin",EGCscal.G);
	wtd(inp.prefix+"Cfold_mesh.bin",EGCfold.C);
	wtd(inp.prefix+"Cscal_mesh.bin",EGCscal.C);
	if (inp.verbosity & WRITEBIT__) {
		// write gnuplot difference file
		writeGnuPlot(EGCfold.E.base_,1e3*(EGCfold.E.base_-EGCscal.E.base_),
			"Ediff",inp.prefix,
			{vb_bnd-.5,cb_bnd+.5,-inp.toldev-5.0,inp.toldev+5.0},
			{vb_bnd,E.vb,E.cb,cb_bnd});
		if (inp.verbosity & VERBOBIT__) os << "gnuplot files using seed '"
			  << inp.prefix << "Ediff' written\n";
	}

	// report results
	if (inp.verbosity & PRINTBIT__)
		os << "\n\ntest yields:\n"
		   << "#bands: " << BLUE__ << EGCfold.Nb() << RESET__ << "\n"
		   << "#kpoints: " << BLUE__ << EGCfold.Nk() << RESET__ << "\n";


	// check valence bands
	double Emax_vb_fold=0.0, Emax_vb_scal=0.0;
	size_t Emax_vb_ind_fold=0, Emax_vb_ind_scal=0;
	if (have_vbands) {
	
		double Edev_vb = 0.0;
		fMat Gmax_vb_fold(EGCfold.dim(),1), Cmax_vb_fold(EGCfold.dim()),
		     Gmax_vb_scal(EGCfold.dim(),1), Cmax_vb_scal(EGCfold.dim());

		// find band maxima and indices
		{
			const auto r = EGCfold.E.base_.crBegin()+cbi-1;
			const auto j = std::max_element(r->begin(),r->cend());
			Emax_vb_fold = *j; Emax_vb_ind_fold = std::distance(r->cbegin(),j);
		}
		{
			const auto r = EGCscal.E.base_.crBegin()+cbi-1;
			const auto j = std::max_element(r->begin(),r->cend());
			Emax_vb_scal = *j; Emax_vb_ind_scal = std::distance(r->cbegin(),j);
		}

		const auto pos_fold = EGCfold.E.iTopos(Emax_vb_ind_fold);
		const auto pos_scal = EGCscal.E.iTopos(Emax_vb_ind_scal);

		// get gradients at energy maxima
		for (size_t m=0; m!=Gmax_vb_fold.size(); ++m)
			Gmax_vb_fold[m] = EGCfold.G.cAt(cat(m,pos_fold))[cbi-1],
			Gmax_vb_scal[m] = EGCscal.G.cAt(cat(m,pos_scal))[cbi-1];

		// get curvature at energy maxima
		for (size_t n=0; n!=Cmax_vb_fold.M(); ++n)
		for (size_t m=0; m!=Cmax_vb_fold.N(); ++m)
			Cmax_vb_fold(m,n) = EGCfold.C.cAt(cat(m,n,pos_fold))[cbi-1],
			Cmax_vb_scal(m,n) = EGCscal.C.cAt(cat(m,n,pos_scal))[cbi-1];

		// find deviations
		for (size_t i=0; i!=EGCfold.Nk(); ++i) {

			// starting index for relevant window
			const size_t jstart = std::distance(EGCfold.E.base_.cAt(i).cbegin(),
				std::lower_bound(EGCfold.E.base_.cAt(i).cbegin(),
				EGCfold.E.base_.cAt(i).cbegin()+cbi,E.vb-inp.Lvb));

			// energy deviations
			for (auto jf=EGCfold.E.base_.cAt(i).cbegin()+jstart,
				  js=EGCscal.E.base_.cAt(i).cbegin()+jstart,
				  jfe=EGCfold.E.base_.cAt(i).cbegin()+cbi; jf!=jfe; ++jf,++js) {
				const double dev = *jf - *js;
				if (std::abs(dev)>std::abs(Edev_vb)) Edev_vb=dev;
			}
		}

		// report results
		if (inp.verbosity & PRINTBIT__)
		os << CYAN__ "\n\nVALENCE BANDS " << RESET__
		   << "[" << 1 << ":" << cbi << "]\n"
		   << std::setprecision(6)
		   << "results at band egde:"
		   << "\n k_fold = " << lm__::T(kmesh.base_.cAt(Emax_vb_ind_fold))
		   << "\n k_scal = " << lm__::T(kmesh.base_.cAt(Emax_vb_ind_scal))
		   << "\n E_fold [eV]  = " << Emax_vb_fold
		   << "\n E_scal [eV]  = " << Emax_vb_scal
		   << "\n E_diff [meV] = " << BLUE__ << 1e3*(Emax_vb_fold-Emax_vb_scal) << RESET__
		   << "\n G_fold [eV*A] = " << lm__::T(Gmax_vb_fold).print(6)
		   << "\n G_scal [eV*A] = " << lm__::T(Gmax_vb_scal).print(6)
		   << "\n diag(C_fold) [eV*A^2] = " << T(diag(Cmax_vb_fold))
		   << "\n diag(C_scal) [eV*A^2] = " << T(diag(Cmax_vb_scal))
		   << "\n diag(C_scal/C_fold)   = "
		   << BLUE__ << T(diag(Cmax_vb_scal)/diag(Cmax_vb_fold)) << RESET__
		   << "\nmax deviations in window [" << std::setprecision(3)
		   << (E.vb-inp.Lvb) << " " << Emax_vb_fold  << "]eV: "
		   << std::setprecision(6) << (Edev_vb*1e-3<inp.toldev ? GREEN__: RED__)
		   << (Edev_vb*1e3) << RESET__ << "meV\n\n";
	}

	
	// check conduction bands
	double Emin_cb_fold=0.0, Emin_cb_scal=0.0;
	size_t Emin_cb_ind_fold=0, Emin_cb_ind_scal=0;
	if (have_cbands) {
	
		double Edev_cb = 0.0;
		fMat Gmin_cb_fold(EGCfold.dim(),1), Cmin_cb_fold(EGCfold.dim()),
		     Gmin_cb_scal(EGCfold.dim(),1), Cmin_cb_scal(EGCfold.dim());

		// find band maxima and indices
		{
			const auto r = EGCfold.E.base_.crBegin()+cbi;
			const auto j = std::min_element(r->begin(),r->cend());
			Emin_cb_fold = *j; Emin_cb_ind_fold = std::distance(r->cbegin(),j);
		}
		{
			const auto r = EGCscal.E.base_.crBegin()+cbi;
			const auto j = std::min_element(r->begin(),r->cend());
			Emin_cb_scal = *j; Emin_cb_ind_scal = std::distance(r->cbegin(),j);
		}
		
		const auto pos_fold = EGCfold.E.iTopos(Emin_cb_ind_fold);
		const auto pos_scal = EGCscal.E.iTopos(Emin_cb_ind_scal);

		// get gradients at energy maxima 
		for (size_t m=0; m!=Gmin_cb_fold.size(); ++m)
			Gmin_cb_fold[m] = EGCfold.G.cAt(cat(m,pos_fold))[cbi],
			Gmin_cb_scal[m] = EGCscal.G.cAt(cat(m,pos_scal))[cbi];

		// get curvature at energy maxima
		for (size_t n=0; n!=Cmin_cb_fold.M(); ++n)
		for (size_t m=0; m!=Cmin_cb_fold.N(); ++m)
			Cmin_cb_fold(m,n) = EGCfold.C.cAt(cat(m,n,pos_fold))[cbi],
			Cmin_cb_scal(m,n) = EGCscal.C.cAt(cat(m,n,pos_scal))[cbi];

		// find deviations
		for (size_t i=0; i!=EGCfold.Nk(); ++i) {

			// ending index for relevant window
			const size_t jend = std::distance(EGCfold.E.base_.cAt(i).cbegin(),
				std::lower_bound(EGCfold.E.base_.cAt(i).cbegin()+cbi,
				EGCfold.E.base_.cAt(i).cend(),E.cb+inp.Lcb));
			
			// energy deviations
			for (auto jf=EGCfold.E.base_.cAt(i).cbegin()+cbi,
				  js=EGCscal.E.base_.cAt(i).cbegin()+cbi,
				  jfe=EGCfold.E.base_.cAt(i).cbegin()+jend; jf!=jfe; ++jf,++js) {
				const double dev = *jf - *js;
				if (std::abs(dev)>std::abs(Edev_cb)) Edev_cb=dev;
			}
		}

		// report results
		if (inp.verbosity & PRINTBIT__)
		os << CYAN__ "\nCONDUCTION BANDS " << RESET__
		   << "[" << cbi+1 << ":" << EGCfold.Nb() << "]\n"
		   << std::setprecision(6)
		   << "results at band egde:"
		   << "\n k_fold = " << lm__::T(kmesh.base_.cAt(Emin_cb_ind_fold))
		   << "\n k_scal = " << lm__::T(kmesh.base_.cAt(Emin_cb_ind_scal))
		   << "\n E_fold [eV]  = " << Emin_cb_fold
		   << "\n E_scal [eV]  = " << Emin_cb_scal
		   << "\n E_diff [meV] = " << BLUE__ << 1e3*(Emin_cb_fold-Emin_cb_scal) << RESET__
		   << "\n G_fold [eV*A] = " << lm__::T(Gmin_cb_fold).print(6)
		   << "\n G_scal [eV*A] = " << lm__::T(Gmin_cb_scal).print(6)
		   << "\n diag(C_fold) [eV*A^2] = " << T(diag(Cmin_cb_fold))
		   << "\n diag(C_scal) [eV*A^2] = " << T(diag(Cmin_cb_scal))
		   << "\n diag(C_scal/C_fold)   = "
		   << BLUE__ << T(diag(Cmin_cb_scal)/diag(Cmin_cb_fold)) << RESET__
		   << "\nmax deviations in window [" << std::setprecision(3)
		   << Emin_cb_fold << " " << (E.cb+inp.Lcb) << "]eV: "
		   << std::setprecision(6) << (Edev_cb*1e-3<inp.toldev ? GREEN__: RED__)
		   << (Edev_cb*1e3) << RESET__ << "meV\n\n";
	}



	if (have_vbands && have_cbands && (inp.verbosity & PRINTBIT__))
		os << std::setprecision(6)
		   << "\nEg folded:" << BLUE__ << (Emin_cb_fold-Emax_vb_fold) << RESET__
		   << "\nEg scaled:" << BLUE__ << (Emin_cb_scal-Emax_vb_scal) << RESET__
		   << "\n\n";
}
void ll__::traceBStest(const ll_hbonds& W, const fMat& EXP, const rv& r,
			 const ll_BStest_input& inp, std::ostream& os) {
	if (!(inp.Nk && inp.verbosity)) return;
	if (inp.Nk>1e6) throw(std::invalid_argument("Nk > 1e6, don't be silly!"));
	if (EXP.empty() || std::abs(det(EXP))<mtol() || !(EXP == round(EXP)))
		throw(std::invalid_argument("bad expansion matrix"));
	if (inp.kpts.empty() || inp.kpts.M()!=W.dim())
		throw(std::invalid_argument("kpts may not be empty and must match dim of cell"));
	
	// write to disk helper lambda
	const auto wtd = [&inp,&os](const std::string& fileName, const auto& F,
				    const std::string& discr = "file") -> void {
		if (~inp.verbosity & PRINTBIT__) return;
		F.writeToFile(fileName);
		if (inp.verbosity & PRINTBIT__)
			os << (discr.empty() ? "": discr+" ") << "'" << fileName << "' written"
		   	   << (inp.verbosity & MD5BIT__ ? ", "+aux::md5(fileName): "") << "\n";
	};

	// generate test cell
	set_mtol(WTOL__); // set tolerance due to wannier input precision
	const auto tCell = W.cell().copy().expand(EXP);
	const auto RB = tCell.getRB();
	reset_mtol();
	
	if (inp.verbosity & PRINTBIT__)
		os << "\n" << CYAN__ << "TRACE BANDSTRUCTURE TEST" << RESET__
		   << "\nexpansion C:\n" << EXP.print(0,1)
		   << "\nreciprocal basis RB:\n"
		   << RB.print(PPREC__,1) << "\n\n";
	if (inp.verbosity & DEBUGBIT__)
		wtd(inp.prefix+"trace_cell.psc",tCell);

	
	const auto hr = readHr(inp.hrdat);
	const auto PP = genPath(inp.kpts,inp.Nk,tCell.getRB());

	// get scaled and folded energies,gradients,curvatures
	const auto EGCscal = inp.re ?
		calcBS_gc(genHam<fMat>(tCell,r,W,true),PP.path(),tCell.B(),inp.Nthreads):
		calcBS_gc(genHam<cMat>(tCell,r,W,true),PP.path(),tCell.B(),inp.Nthreads);
	set_mtol(WTOL__);
	const auto EGCfold = calcFoldedBS_gc(hr,PP.path(),tCell.B(),W.cell().B(),inp.Nthreads);
	reset_mtol();

	// write files
	wtd(inp.prefix+"path.bin",PP.path());
	wtd(inp.prefix+"pos.bin",PP.pos());
	wtd(inp.prefix+"Efold.bin",EGCfold.E);
	wtd(inp.prefix+"Escal.bin",EGCscal.E);
	wtd(inp.prefix+"Gfold.bin",EGCfold.G);
	wtd(inp.prefix+"Gscal.bin",EGCscal.G);
	wtd(inp.prefix+"Cfold.bin",EGCfold.C);
	wtd(inp.prefix+"Cscal.bin",EGCscal.C);
}
