// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_OMEN_
#define _LL_OMEN_

#include "ll_types.h"
#include "ll_BStest.h"
#include "ll_hbonds.h"
#include "ll_hbondss.h"
#include "ll_hio.h"


/** input for running winterface
 * As a child of aux_parser, it redefines its own versions of
 * parseKey_ and printHelp_.
 */
struct ll_wf_input: public virtual ll_wmatching_input,
		      public virtual ll_BStest_input,
		      public virtual ll_hbondss_input {
public:
	//! resolve ambiguities
	using ll_wmatching_input::hrdat;


	/** @name input for basic electronic properties
	 */
	std::string weig = WEIG__;			//!< wannier90.eig file
	std::string outcar = OUTCAR__;			//!< VASP OUTCAR file
	double Ef = nan("");				//!< Fermi energy (if no OUTCAR is present)


	/** @name orthorhombic cell (ORC) detection
	 * and next neighbour (NN) expansion
	 */
	double itol = 5e-2;				//!< tolerance used in fraction detection
	fMat xyz = lm__::eye<fMat>(DIM__,DIM__);	//!< xyz direction matrix
	fMat S = {};					//!< stress tensor
	fMat C = {};					//!< expansion to ORC manually
	std::vector<size_t> l =				//!< expansion to NN manually
		std::vector<size_t>(DIM__,0);
	double vac = 0.0;				//!< vacuum modification length
	ll__::rv filter_wbh = ll__::rv(DIM__,true);	//!< filter bonds in wbh (in each direction)
							//!< to match expansion in nn-only ORC


	/** @name BS test parameters specific to OMEN
	 */
	size_t expand_mesh = 1;				//!< level of cell expansion for mesh mode
	size_t expand_trace = 1;			//!< level of cell expansion for trace mode


	/** @name hamiltonian scaling
	 */
	double IR = nan("");				//!< interaction radius cutoff manually
	bool strict_matching = true;			//!< allow only exact type matches
	double device_length = 400.0;			//!< OMEN device length [ANGSTROM] for input stump
	bool force_check = false;			//!< force consistency check


public:
	//! parseKey_ redefinition
	struct parseKey_: public virtual ll_wmatching_input::parseKey_,
			  public virtual ll_BStest_input::parseKey_,
			  public virtual ll_hbondss_input::parseKey_ {
	/** constructor defining keys
	 */
	inline explicit parseKey_(ll_wf_input& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator& i, std::sregex_token_iterator& e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e),
			ll_wmatching_input::parseKey_(p,key,file,lcnt,i,e),
			ll_BStest_input::parseKey_(p,key,file,lcnt,i,e),
			ll_hbondss_input::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			// input for basis electronic properties
			case "weig"_h: PARSE__(p.weig); return;
			case "outcar"_h: PARSE__(p.outcar); return;
			case "Ef"_h: PARSE__(p.Ef); return;

			// ORC and NN
			case "itol"_h: PARSE__(p.itol); return;
			case "xyz"_h: PARSE__(p.xyz); return;
			case "S"_h: PARSE__(p.S); return;
			case "C"_h: PARSE__(p.C); return;
			case "l"_h: PARSE__(p.l,false,DIM__); return;
			case "vac"_h: PARSE__(p.vac); return;
			case "filter_wbh"_h: PARSE__(p.filter_wbh,DIM__); return;

			// BS parameters
			case "expand_mesh"_h: PARSE__(p.expand_mesh); return;
			case "expand_trace"_h: PARSE__(p.expand_trace); return;

			// hamiltonian scaling
			case "IR"_h: PARSE__(p.IR); return;
			case "strict_matching"_h: PARSE__(p.strict_matching); return;
			case "device_length"_h: PARSE__(p.device_length); return;
			case "force_check"_h: PARSE__(p.force_check); return;
		}
	}
	};

public:
	/** printHelp_ redinition
	 */
	struct printHelp_: public virtual ll_wmatching_input::printHelp_,
			   public virtual ll_BStest_input::printHelp_,
			   public virtual ll_hbondss_input::printHelp_ {
	/** constructor defining help messages
	 */
	inline explicit printHelp_(const ll_wf_input& p, std::ostream& os):
			aux_parser::printHelp_(p,os),
			ll_wmatching_input::printHelp_(p,os),
			ll_BStest_input::printHelp_(p,os),
			ll_hbondss_input::printHelp_(p,os) {
		printHelpTuple_(os,std::make_tuple(
		
		TOPIC__("OMEN: ELECTRONIC PROPERTIES")
		,HELP__("weig",p.weig,
			"name of the wannier90 *.weig file")
		,HELP__("outcar",p.outcar,
			"name of the VASP OUTCAR file")
		,HELP__("Ef",p.Ef,
			"manual specification of the Fermi energy")

		,TOPIC__("OMEN: BASIS EXPANSION")
		,HELP__("itol",p.itol,
			"cutoff tolerance used in fraction detection, higher values may\n"
			"help with automatic basis detection")
		,HELP__("xyz",p.xyz,DIM__,DIM__,
			"matrix to specify where x,y and z directions should point to, in cartesian.\n"
			"this is used to find the expansion from the primitive cell to the smallest\n"
			"orthorhombic cell, rows must be perpendicular")
		,HELP__("S",p.S,DIM__,DIM__,
			"stress tensor to be applied to primitive basis before anything else\n"
			"this is useful in case no exact orthorhombic supercell exists\n"
			"beware! there are no enforced boundaries, do not use this for anything else!\n"
			"as a warning frobenius(eye-S) is printed, if this is not small you are doing\n"
			"something fishy")
		,HELP__("C",p.C,DIM__,DIM__,
			"expansion coefficient matrix overwrite for expansion from primitive to\n"
			"minimal orthorhombic, expansion for each new vector in terms of the old\n"
			"ones must be given in the rows")
		,HELP__("l",p.l,DIM__,
			"expansion coefficients overwrite for minimal orthorhombic basis to NN only")
		,HELP__("vac",p.vac,
			"vacuum modification length, replaces large value used in DFT\n"
			"if this value is lower than a bond length, the bond length is used")
		,HELP__("filter_wbh",p.filter_wbh,
			"switch to turn on/off filtering the bonds in the wbh in each direction\n"
			"so as to match the expansion to nn-only ORC")

		,TOPIC__("OMEN: BSTEST CELL EXPANSION")
		,HELP__("expand_mesh",p.expand_mesh,
			"cell expansion level for the BS test in mesh mode:\n"
			" 0: no expansion, test on primitive cell\n"
			" 1: expansion to smallest orthorhombic cell\n"
			" 2: expansion to NN only cell")
		,HELP__("expand_trace",p.expand_trace,
			"cell expansion level for the BS test in trace mode:\n"
			" 0: no expansion, test on primitive cell\n"
			" 1: expansion to smallest orthorhombic cell\n"
			" 2: expansion to NN only cell")
		
		,TOPIC__("OMEN: HAMILTONIAN SCALING")
		,HELP__("IR",p.IR,
			"interaction radius cutoff, any bonds of length above this value\n"
			"will be discarded when scaling hamiltonians. If this is not set manually\n"
			"it will be set such that all bonds are included")
		,HELP__("strict_matching",p.strict_matching,
			" - true: allow only exact type {i1,i2} matches\n"
			" - false: allow {i1,?} matches")
		,HELP__("device_length,",p.device_length,
			"total device length [ANGSTROM] to use in the OMEN stump cmd file")
		,HELP__("force_check",p.force_check,
			"switch to force the consistency check of hamiltonian matrices\n"
			"after writing them to disk")
		));
	}
	};
};



//! struct containing all the relevant information for phinterface
struct ll_ph_input: public virtual aux_parser {
public:
	/** @name filenames
	 */
	std::string myfile = "this_n_that";	//!< my filename

	// ...

public:
	//! parseKey_ redefinition
	struct parseKey_: public virtual aux_parser::parseKey_ {
	/** constructor defining keys
	 */
	inline explicit parseKey_(ll_ph_input& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator& i, std::sregex_token_iterator& e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			// filenames
			case "myfile"_h: PARSE__(p.myfile); return;
		}
	}
	};

public:
	//! printHelp_ redefinition
	struct printHelp_: public virtual aux_parser::printHelp_ {
	/** constructor defining help messages
	 */
	inline explicit printHelp_(const ll_ph_input& p, std::ostream& os):
			aux_parser::printHelp_(p,os) {
		printHelpTuple_(os,std::make_tuple(
		
		TOPIC__("OMEN: FILENAMES")
		,HELP__("myfile",p.myfile,
			"myfile n stuff")
		));
	}
	};
};


namespace ll__ {
namespace omen {

	/* helper functions
	 */
	//! OMEN style Hamiltonian filename form R vector
	inline std::string H_fileName(const lm__::fArray& R) noexcept {
		assert(R.size()==DIM__);
		return "H_"+std::to_string((DIM__*DIM__)/2
			   + int(R[0]*DIM__*DIM__ + R[1]*DIM__ + R[2]))
			   + ".bin";
	}
	//! OMEN style R vectors according to spacial restrictions
	inline fMat rToR(const rv& r) noexcept {
		assert(r.size()==DIM__);
	
		// get bounds and total number of shift vectors
		fMat bnds(DIM__,1);
		for (size_t d=0; d!=DIM__; ++d)
			bnds[d] = r[d] ? 0.0: 1.0;
		
		// result, allocate space
		fMat res(DIM__,0); res.reserve(size_t(std::pow(DIM__,sum(bnds))));
		
		// generate R vecs
		fMat vec(DIM__,1);
		for (vec[0] = -bnds[0]; vec[0]<=bnds[0]; ++vec[0])
		for (vec[1] = -bnds[1]; vec[1]<=bnds[1]; ++vec[1])
		for (vec[2] = -bnds[2]; vec[2]<=bnds[2]; ++vec[2])
			res.push_back(vec);

		return res;
	}


	//! OMEN style binary sparse hamiltonian writer
	class writer final: public ll_writer {
	public:
		/** @name constructors
		 */
		/** constructor from R vectors and number of Wannier functions
		 * @param R R vectors
		 * @param Nw number of Wannier functions
		 * @param prefix prefix to put in front of filenames
		 */
		inline explicit writer(lm__::fMat R, const size_t Nw, std::string prefix="./"):
				ll_writer(std::move(R)), Nw_(Nw), prefix_(std::move(prefix)) {
			assert(dim()==DIM__);
		}


		/** @name general information
		 */
		//! number of Wannier functions
		inline size_t Nw() const noexcept { return Nw_; }
		//! filename for an R vector
		inline std::string fileName(const lm__::fArray& R) const noexcept {
			return prefix()+H_fileName(R);
		}
		//! prefix to put in front of filenames
		inline const std::string& prefix() const noexcept { return prefix_; }


		/** current block information
		 */
		//! the current filename
		inline std::string c_fileName() const noexcept { return fileName(c_R()); }
		//! description of the current block
		inline std::string c_descr() const noexcept { return c_fileName(); }


		/** @name initiate and finalize
		 */
		//! new Hamiltonian block initialization
		inline void newBlock() {
			// open file
			file_.open(c_fileName(), std::ios::binary | std::ios::out);
			if (!file_.good())
				throw(std::invalid_argument("failed to open '"+c_fileName()+"'"));
			
			// write fake header
			double head[3] = {0.0,0.0,0.0};
			file_.write((char*) &head, 3*sizeof(double));
		}

	protected:
		/** internals
		 */
		//! the number of bytes for each Hamiltonian entry
		inline size_t bytes_() const noexcept { return 2*sizeof(double)+sizeof(ll__::hel); }
		/** write data blob function
		 * @param inp vector of Hamiltonian entries
		 */
		inline void insert_(const hdat& inp) {
			for (const auto& i: inp) {
				double buff[] = {double(i.m),double(i.n),std::real(i.h),std::imag(i.h)};
				file_.write((char*) &buff, 4*sizeof(double));
			}
		}
		//! flush current Hamiltonian block function
		inline void flush_() {
			// close and reopen file
			file_.close();
			file_.open(c_fileName(), std::ios::binary | std::ios::out | std::ios::in);

			// write real header and close
			double head[3] = {double(Nw()),double(c_nnz()),0.0};
			file_.write((char*) &head, 3*sizeof(double));
			file_.close();
		}
	

	protected:
		/** member variables
		 */
		std::fstream file_;		//!< file stream
		size_t Nw_;			//!< number of Wannier functions
		std::string prefix_;		//!< prefix to put in front of filenames
	};

	
	
	/* running winterface
	 */
	//! struct containing all the data needed for OMEN input
	struct prepper {
	public:
		ll_hbonds W;		//!< wbh
		ll_cell ORcell;		//!< ORC cell
		rv r;			//!< restriction vector
		fMat C;			//!< expansion to ORC
		fMat NNE;		//!< next neighbor expansion matrix
		
		/** @name constructors
		 */
		/** constructor from user input
		 * @param inp user input
		 * @param os stream to print into
		 */
		prepper(const ll_wf_input& inp, std::ostream& os);
	};

	/* io
	 */
	/** read OMEN style sparse matrices into H_R
	 * @param prefix prefix to add in front of filenames
	 * @param r restriction vector
	 */
	inline R_H<ll_sparse> readOMENhr(const std::string& prefix, const rv& r) {
		
		const fMat R = rToR(r);
		std::vector<ll_sparse> H; H.reserve(R.N());
		for (auto i=R.ccBegin(),e=R.ccEnd(); i!=e; ++i)
			H.push_back(ll_sparse(prefix+H_fileName(*i)));
		
		return {std::move(R),std::move(H)};
	}
	/** write structural input
	 * @param inp user input
	 * @param os stream to print into
	 */
	void writeInput(const ll_wf_input& inp, std::ostream& os);
	/** write hamiltonian matrices
	 * @param inp user input
	 * @param LM layer matrix
	 * @param id id strings corresponding to LM
	 * @param L matrix containing the periodity lengths
	 * @param os stream to print into
	 */
	void hctor(const ll_wf_input& inp, const fMat& LM, const idv& id,
			const fMat& L, std::ostream& os);


	/* OMEN interface
	 */
	//! OMEN input generator, to be used from within OMEN
	void callInputGenerator();
	//! OMEN Hamiltonian constructor, to be used from within OMEN
	void callHamiltonianConstructor(double* const Layer_Matrix, const idv& id,
					const double Lx, const double Ly, const double Lz);
	//! write 1D minimal OMEN input script
	void writeStump1D(const lm__::fMat& B, const double L, const double ftc=.9);
	//! write 2D minimal OMEN input script
	void writeStump2D(const lm__::fMat& B, const double L, const double ftc=.9);
	//! write 3D minimal OMEN input script
	void writeStump3D(const lm__::fMat& B, const double L, const double ftc=.9);


	/* print fancy titles
	 */
	//! winterface title
	inline void hw_winterFace(std::ostream& os) noexcept {
		os
		<< CYAN__ << "\n"
		<< " █     █░ ██▓ ███▄    █ ▄▄▄█████▓▓█████  ██▀███    █████▒▄▄▄       ▄████▄  ▓█████ \n"
		<< "▓█░ █ ░█░▓██▒ ██ ▀█   █ ▓  ██▒ ▓▒▓█   ▀ ▓██ ▒ ██▒▓██   ▒▒████▄    ▒██▀ ▀█  ▓█   ▀ \n"
		<< "▒█░ █ ░█ ▒██▒▓██  ▀█ ██▒▒ ▓██░ ▒░▒███   ▓██ ░▄█ ▒▒████ ░▒██  ▀█▄  ▒▓█    ▄ ▒███   \n"
		<< "░█░ █ ░█ ░██░▓██▒  ▐▌██▒░ ▓██▓ ░ ▒▓█  ▄ ▒██▀▀█▄  ░▓█▒  ░░██▄▄▄▄██ ▒▓▓▄ ▄██▒▒▓█  ▄ \n"
		<< "░░██▒██▓ ░██░▒██░   ▓██░  ▒██▒ ░ ░▒████▒░██▓ ▒██▒░▒█░    ▓█   ▓██▒▒ ▓███▀ ░░▒████▒\n"
		<< "░ ▓░▒ ▒  ░▓  ░ ▒░   ▒ ▒   ▒ ░░   ░░ ▒░ ░░ ▒▓ ░▒▓░ ▒ ░    ▒▒   ▓▒█░░ ░▒ ▒  ░░░ ▒░ ░\n"
		<< "  ▒ ░ ░   ▒ ░░ ░░   ░ ▒░    ░     ░ ░  ░  ░▒ ░ ▒░ ░       ▒   ▒▒ ░  ░  ▒    ░ ░  ░\n"
		<< "  ░   ░   ▒ ░   ░   ░ ░   ░         ░     ░░   ░  ░ ░     ░   ▒   ░           ░   \n"
		<< "    ░     ░           ░             ░  ░   ░                  ░  ░░ ░         ░  ░\n"
		<< "                                                                  ░               \n"
		<< RESET__ << "compiled on: " << __TIME__ << ", " << __DATE__ << "\n";
	}
	//! winterface goodbye
	inline void gb_winterFace(std::ostream& os) noexcept {
		os
		<< CYAN__ << "\n"
		<< " ▄████  ▒█████   ▒█████  ▓█████▄     ▄▄▄▄ ▓██   ██▓▓█████  \n"
		<< " ██▒ ▀█▒▒██▒  ██▒▒██▒  ██▒▒██▀ ██▌   ▓█████▄▒██  ██▒▓█   ▀ \n"
		<< "▒██░▄▄▄░▒██░  ██▒▒██░  ██▒░██   █▌   ▒██▒ ▄██▒██ ██░▒███   \n"
		<< "░▓█  ██▓▒██   ██░▒██   ██░░▓█▄   ▌   ▒██░█▀  ░ ▐██▓░▒▓█  ▄ \n"
		<< "░▒▓███▀▒░ ████▓▒░░ ████▓▒░░▒████▓    ░▓█  ▀█▓░ ██▒▓░░▒████▒\n"
		<< " ░▒   ▒ ░ ▒░▒░▒░ ░ ▒░▒░▒░  ▒▒▓  ▒    ░▒▓███▀▒ ██▒▒▒ ░░ ▒░ ░\n"
		<< "  ░   ░   ░ ▒ ▒░   ░ ▒ ▒░  ░ ▒  ▒    ▒░▒   ░▓██ ░▒░  ░ ░  ░\n"
		<< "░ ░   ░ ░ ░ ░ ▒  ░ ░ ░ ▒   ░ ░  ░     ░    ░▒ ▒ ░░     ░   \n"
		<< "      ░     ░ ░      ░ ░     ░        ░     ░ ░        ░  ░\n"
		<< "                           ░               ░░ ░            \n"
		<< RESET__;
	}
	//! phinterface title
	inline void hw_phinterFace(std::ostream& os) noexcept {
		os
		<< CYAN__ << "\n"
		<< "██████╗ ██╗  ██╗██╗███╗   ██╗████████╗███████╗██████╗ ███████╗ █████╗  ██████╗███████╗  \n"
		<< "██╔══██╗██║  ██║██║████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔════╝██╔══██╗██╔════╝██╔════╝  \n"
		<< "██████╔╝███████║██║██╔██╗ ██║   ██║   █████╗  ██████╔╝█████╗  ███████║██║     █████╗    \n"
		<< "██╔═══╝ ██╔══██║██║██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██╔══╝  ██╔══██║██║     ██╔══╝    \n"
		<< "██║     ██║  ██║██║██║ ╚████║   ██║   ███████╗██║  ██║██║     ██║  ██║╚██████╗███████╗  \n"
		<< "╚═╝     ╚═╝  ╚═╝╚═╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝ ╚═════╝╚══════╝  \n"
		<< RESET__ << "compiled on: " << __TIME__ << ", " << __DATE__ << "\n";
	}
	//! phinterface goodbye
	inline void gb_phinterFace(std::ostream& os) noexcept {
		os
		<< CYAN__ << "\n"
		<< " ██████╗  ██████╗  ██████╗ ██████╗     ██████╗ ██╗   ██╗███████╗ \n"
		<< "██╔════╝ ██╔═══██╗██╔═══██╗██╔══██╗    ██╔══██╗╚██╗ ██╔╝██╔════╝ \n"
		<< "██║  ███╗██║   ██║██║   ██║██║  ██║    ██████╔╝ ╚████╔╝ █████╗   \n"
		<< "██║   ██║██║   ██║██║   ██║██║  ██║    ██╔══██╗  ╚██╔╝  ██╔══╝   \n"
		<< "╚██████╔╝╚██████╔╝╚██████╔╝██████╔╝    ██████╔╝   ██║   ███████╗ \n"
		<< " ╚═════╝  ╚═════╝  ╚═════╝ ╚═════╝     ╚═════╝    ╚═╝   ╚══════╝ \n"
		<< RESET__;
	}
}
}

#endif // _LL_OMEN_

/** @}
 */
