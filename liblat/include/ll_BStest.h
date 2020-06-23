// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_BSTEST_
#define _LL_BSTEST_

#include "aux_parser.h"
#include "ll_defs.h"
#include "ll_types.h"
#include "ll_fn.h"

class ll_hbonds;

/** input for bandstructure tests.
 * As a child of aux_parser, it redefines its own versions of
 * parseKey_ and printHelp_.
 */
struct ll_BStest_input: public virtual aux_parser {
public:
	/** @name types
	 */
	typedef lm__::fMat fMat;			//!< real matrix


public:
	/** @name mesh specific
	 */
	std::vector<double> bzbounds{.0,.5};		//!< bounds inside the BZ used in the mesh
	double rho_k = 1000.0;				//!< kpoint density for check in BZ
	std::string maj_style = "MATLAB";		//!< majority ordering style

	/** @name trace specific
	 */
	size_t Nk = 101;				//!< number of points along kpts
	fMat kpts = lm__::ncat(				// kpoints
		-.5*lm__::cId<fMat>(DIM__,0),		// left edge of Brillouin zone, x axis
		 .5*lm__::cId<fMat>(DIM__,0)		// right edge of Brillouin zone, x axis
	);						//!< kpoints for the trace
	
	/** @name common options
	 */
	std::string hrdat = HR__;			//!< hrdat to use for folding algorithm
	bool re = true;					//!< take only the real part in scaled blocks
	double Lvb = 2.0;				//!< region inside valence band to check
	double Lcb = 2.0;				//!< region inside conduction band to check
	double toldev = 30.0;				//!< tolerable deviation scaled vs folded in meV
	size_t Nthreads = 4;				//!< number of threads to use


public:
	/** parseKey_ redefinition
	 */
	struct parseKey_: public virtual aux_parser::parseKey_ {
	/** constructor defining keys
	 */
	inline explicit parseKey_(ll_BStest_input& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator i, std::sregex_token_iterator e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;

		switch(key) {
			// mesh specific
			case "bzbounds"_h: PARSE__(p.bzbounds,2); return;
			case "rho_k"_h: PARSE__(p.rho_k); return;
			case "maj_style"_h: PARSE__(p.maj_style); return;

			// trace specific
			case "Nk"_h: PARSE__(p.Nk); return;
			case "kpts"_h: PARSE__(p.kpts); return;

			// common options
			case "hrdat"_h: PARSE__(p.hrdat); return;
			case "re"_h: PARSE__(p.re); return;
			case "Lvb"_h: PARSE__(p.Lvb); return;
			case "Lcb"_h: PARSE__(p.Lcb); return;
			case "toldev"_h: PARSE__(p.toldev); return;
			case "Nthreads"_h: PARSE__(p.Nthreads); return;
		}
	}
	};

public:
	/** printHelp_ redefinition
	 */
	struct printHelp_: public virtual aux_parser::printHelp_ {
	/** constructor defining help messages
	 */
	inline explicit printHelp_(const ll_BStest_input& p, std::ostream& os):
			aux_parser::printHelp_(p,os) {
		printHelpTuple_(os,std::make_tuple(
		
		TOPIC__("BSTEST: MESH SPECIFIC")
		,HELP__("bzbounds",p.bzbounds,
			"bounds to use when generating the mesh, in reciprocal basis")
		,HELP__("rho_k",p.rho_k,
			"density of kpoints points in the Brillouin Zone of the NN only orthorhmobic cell\n"
			"a meshgrid approximately matching this density is used, in units A^-dim\n"
			"setting this to .0 will turn off bandstructure check")
		,HELP__("maj_style",p.maj_style,
			"setting this to MATLAB for a mesh like that from the meshgrid function\n"
			" - inorder: 0,1,2,...\n"
			" - MATLAB: 1,0,2,...")
		
		,TOPIC__("BSTEST: TRACE SPECIFIC")
		,HELP__("Nk",p.Nk,
			"number of points along the additional ktrace bandstructure calculation\n"
			"setting this to 0 will turn it off")
		,HELP__("kpts",p.kpts,DIM__,NPOS__,
			"kpoints along which to generate the ktrace")
		
		,TOPIC__("BSTEST: COMMON OPTIONS")
		,HELP__("hrdat",p.hrdat,
			"wannier90 hamiltonian data file to use for folding alorithms")
		,HELP__("re",p.re,
			"switch to use the real part only of the wannier90 hamiltonian data\n"
			"durch bandstructure calculations")
		,HELP__("Lvb",p.Lvb,
			"critical depth in among the valence bands, this is used in warning output\n"
			"as well as some generated files, such which bandstructure data to include\n"
			"in [eV]")
		,HELP__("Lcb",p.Lcb,
			"analog to Lvb for conduction bands")
		,HELP__("toldev",p.toldev,
			"tolerable mismatch between folding and scaling algorithms, affects warning\n"
			"message only, in [meV]")
		,HELP__("Nthreads",p.Nthreads,
			"number of threads to use when copmputing bandstructure")
		));
	}
	};
};


namespace ll__ {
	/** bandstructure test on a mesh in the Brillouin zone.
	 * @param W wbh to test
	 * @param EXP expansion matrix to test using a supercell, allows for
	 * 		testing the upscaling algorithm.
	 * @param r spacial restriction vector
	 * @param E valence and conduction band edges
	 * @param inp user input struct
	 * @param os stream to print into
	 */
	void meshBStest(const ll_hbonds& W, const fMat& EXP, const rv& r,
			const vb_cb& E, const ll_BStest_input& inp, std::ostream& os);
	/** bandstructure test on a trace in k-space. No automatic checking done here.
	 * @param W wbh to test
	 * @param EXP expansion matrix to test using a supercell, allows for
	 * 		testing the upscaling algorithm.
	 * @param r spacial restriction vector
	 * @param inp user input struct
	 * @param os stream to print into
	 */
	void traceBStest(const ll_hbonds& W, const fMat& EXP, const rv& r,
			const ll_BStest_input& inp, std::ostream& os);
}

#endif // _LL_BSTEST_

/** @}
 */
