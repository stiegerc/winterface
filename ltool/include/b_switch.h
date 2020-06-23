// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup ltool
 * @{
 */

#ifndef _B_SWITCH_
#define _B_SWITCH_

#include "aux_parser.h"
#include "ll_defs.h"
#include "ll_types.h"
#include "ll_hbondss.h"


/** input for bandstructure calculations.
 * As a child of aux_parser, it redefines its own versions of
 * parseKey_ and printHelp_.
 */
struct b_input: public virtual aux_parser,
       		public virtual ll_hbondss_input	{
public:
	/** @name filenames
	 */
	std::string wout = WOUT__;			//!< wannier90.wout file
	std::string hrdat = HR__;			//!< wannier90_hr.dat file
	std::string lattice_dat = "";			//!< OMEN style lattice file
	std::string mat_par = "ph_mat_par";		//!< OMEN style material file
	std::string inprefix = "./";			//!< prefix for input files
	std::string hr = "wbh.hr";			//!< hr file
	std::string layer_matrix = "Layer_Matrix.dat";	//!< layer matrix file


	/** @name parameters
	 */
	std::string mode = "fold";			//!< computation mode
	ll__::rv r = {};				//!< restriction vector
	lm__::fMat C = lm__::eye<lm__::fMat>(DIM__);	//!< expansion matrix
	lm__::fMat R = lm__::fMat(DIM__,0);		//!< R vectors
	lm__::fMat kpts = lm__::fMat({-.5,.0,.0,
			         .5,.0,.0},DIM__,2);	//!< kpoints aling which a trace is generated
	lm__::fMat k = {};				//!< kpoints manual specification
	bool strict_matching = true;			//!< strict matching switch
	std::vector<size_t> Nk = {100};			//!< number of k-points along trace
	double rho_k = 1000.0;				//!< kpoint density for mesh generation
	std::string maj_style = "MATLAB";		//!< majority style in the mesh
	std::vector<double> bzbounds = {-.5,.5};	//!< BZ bounds
	bool re = false;				//!< use only real parts of Hamiltonians
	size_t Nthreads = 0;				//!< number of threads for parallel sections
	bool dump_hamiltonians = false;			//!< switch to dump generated Hamiltonian matrices


	/** @name local bandstructure
	 */
	lm__::fMat LB = {};				//!< local basis for local bandstructure


public:
	//! parseKey_ redefinition
	struct parseKey_: public virtual aux_parser::parseKey_,
			  public virtual ll_hbondss_input::parseKey_ {

	/** constructor defining keys
	 */
	inline explicit parseKey_(b_input& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator& i, std::sregex_token_iterator& e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e),
			ll_hbondss_input::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			// filenames
			case "wout"_h: PARSE__(p.wout); return;
			case "hrdat"_h: PARSE__(p.hrdat); return;
			case "lattice_dat"_h: PARSE__(p.lattice_dat); return;
			case "mat_par"_h: PARSE__(p.mat_par); return;
			case "inprefix"_h: PARSE__(p.inprefix); return;
			case "hr"_h: PARSE__(p.hr); return;
			case "layer_matrix"_h: PARSE__(p.layer_matrix); return;

			// parameters
			case "mode"_h: PARSE__(p.mode); return;
			case "r"_h: PARSE__(p.r); return;
			case "C"_h: PARSE__(p.C); return;
			case "R"_h: PARSE__(p.R); return;
			case "kpts"_h: PARSE__(p.kpts); return;
			case "k"_h: PARSE__(p.k); return;
			case "strict_matching"_h: PARSE__(p.strict_matching); return;
			case "Nk"_h: PARSE__(p.Nk); return;
			case "rho_k"_h: PARSE__(p.rho_k); return;
			case "maj_style"_h: PARSE__(p.maj_style); return;
			case "bzbounds"_h: PARSE__(p.bzbounds,2); return;
			case "re"_h: PARSE__(p.re); return;
			case "Nthreads"_h: PARSE__(p.Nthreads); return;
			case "dump_hamiltonians"_h: PARSE__(p.dump_hamiltonians); return;

			// local bandstructure
			case "LB"_h: PARSE__(p.LB); return;
		}
	}
	};


public:
	//! printHelp_ redefinition
	struct printHelp_: public virtual aux_parser::printHelp_,
		           public virtual ll_hbondss_input::printHelp_ {

	/** constructor defining help messages
	 */
	inline explicit printHelp_(const b_input& p, std::ostream& os):
			aux_parser::printHelp_(p,os),
	       		ll_hbondss_input::printHelp_(p,os) {

	printHelpTuple_(os,std::make_tuple(
		TOPIC__("-b: FILENAMES")
		,HELP__("wout",p.wout,
			"wannier90 wout filename")
		,HELP__("hrdat",p.hrdat,
			"wannier90 hrdat filename")
		,HELP__("lattice_dat",p.lattice_dat,
			"OMEN lattice file")
		,HELP__("mat_par",p.mat_par,
			"OMEN material file")
		,HELP__("inprefix",p.inprefix,
			"prefix for OMEN material, lattice and H_* files in legacy mode")
		,HELP__("hr",p.hr,
			"input file for hr mode")
		,HELP__("layer_matrix",p.layer_matrix,
			"layer matrix for which the hr was generated")
		
		,TOPIC__("-b: PARAMETERS")
		,HELP__("mode",p.mode,
			"mode switch:"
			" - fold: use wout/hrdat\n"
			" - scale: use lattice_dat/wbh\n"
			" - legacy: compute bandstructure from OMEN input files\n"
			" - local: compute local approximate bandstructures")
		,HELP__("r",p.r,
			"restriction vector")
		,HELP__("C",p.C,
			"expansion matrix for mode 1")
		,HELP__("R",p.R,
			"R vectors for mode 2")
		,HELP__("kpts",p.kpts,
			"list of k vectors in direct coordinates")
		,HELP__("k",p.k,
			"directly specified k points")
		,HELP__("strict_matching",p.strict_matching,
			"allow strict type matching only")
		,HELP__("Nk",p.Nk,
			"number of points along k trace")
		,HELP__("rho_k",p.rho_k,
			"density of the mesh in the BZ")
		,HELP__("maj_style",p.maj_style,
			"majority style:\n"
			" - inorder: in order, 0,1,2,...\n"
			" - MATLAB: MATLAB style, 1,0,2,...")
		,HELP__("bzbounds",p.bzbounds,
			"bounds for the mesh in reciprocal basis")
		,HELP__("re",p.re,
			"skip imaginary parts")
		,HELP__("Nthreads",p.Nthreads,
			"number of threads to use for bandstructure calculations and\n"
			"hamiltonian matrix generation\n"
			"setting of 0 will set it to number of cores")
		,HELP__("dump_hamiltonians",p.dump_hamiltonians,
			"switch whether to write the Hamiltonian matrices used in\n"
			"bandstructure calculations, not applicable when using\n"
			"the folding mode")
		
		,TOPIC__("-b: LOCAL BANDSTRUCTURE")
		,HELP__("LB",p.LB,
			"basis to use for local bandstructures\n"
			"if this is unspecified the code will attempt to find\n"
			"the primitive basis on its own using sptol")
	));
	}
	};
};


/** calculate bandstructures
 * @param inp user input
 * @param os stream to print into
 */
void b_switch(const b_input& inp, std::ostream& os=std::cout);

#endif // _B_SWITCH_

/** @}
 */
