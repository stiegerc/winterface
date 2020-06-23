// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup ltool
 * @{
 */

#ifndef _H_SWITCH_
#define _H_SWITCH_

#include <iostream>
#include "ll_omen.h"
#include "ll_defs.h"
#include "ll_types.h"
#include <cfloat>


/** input to generate Hamiltonian matrices
 * As a child of aux_parser, it redefines its own versions of
 * parseKey_ and printHelp_.
 */
struct h_input: public virtual ll_wf_input {
public:
	/** @name filenames
	 */
	std::string layer_matrix = "Layer_Matrix.dat";		//!< layer matrix file


	/** @name parameters
	 */
	ll__::rv r = ll__::rv(DIM__,false);			//!< restriction vector
	fMat device_definition = {};				//!< specifies cuboid regions, empty -> all space
	std::string hr_out = "wbh.hr";				//!< fileName for the output Hamiltonian
	std::string hr_format = "omen";				//!< output format
	size_t pprec = 6;					//!< floating point precision
	fMat R = {};						//!< manually specify R grid


public:
	//! parseKey_ redefinition
	struct parseKey_: public virtual ll_wf_input::parseKey_ {

	/** constructor defining keys
	 */
	inline explicit parseKey_(h_input& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator& i, std::sregex_token_iterator& e):
	       		aux_parser::parseKey_(p,key,file,lcnt,i,e),
	       		ll_wmatching_input::parseKey_(p,key,file,lcnt,i,e),
	       		ll_BStest_input::parseKey_(p,key,file,lcnt,i,e),
	       		ll_hbondss_input::parseKey_(p,key,file,lcnt,i,e),
			ll_wf_input::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			// filenames
			case "layer_matrix"_h: PARSE__(p.layer_matrix); return;

			// parameters
			case "r"_h: PARSE__(p.r,DIM__); return;
			case "device_definition"_h: PARSE__(p.device_definition); return;
			case "hr_out"_h: PARSE__(p.hr_out); return;
			case "hr_format"_h: PARSE__(p.hr_format); return;
			case "pprec"_h: PARSE__(p.pprec); return;
			case "R"_h: PARSE__(p.R); return;
		}
	}
	};


public:
	//! printHelp_ redefinition
	struct printHelp_: public virtual ll_wf_input::printHelp_ {

	/** constructor defining help messages
	 */
	inline explicit printHelp_(const h_input& p, std::ostream& os):
			aux_parser::printHelp_(p,os),
			ll_wmatching_input::printHelp_(p,os),
			ll_BStest_input::printHelp_(p,os),
			ll_hbondss_input::printHelp_(p,os),
			ll_wf_input::printHelp_(p,os) {

	printHelpTuple_(os,std::make_tuple(
		TOPIC__("-h: FILENAMES")
		,HELP__("layer_matrix",p.layer_matrix,
			"layer_matrix file")
		
		,TOPIC__("-h: PARAMETERS")
		,HELP__("r",p.r,
			"restricted directions vector")
		,HELP__("device_definition",p.device_definition,
			"specify cuboids whose insides are included in the\n"
			"[xmin1 ymin1 zmin1]\n"
		        "[xmax1 ymax1 zmax1]\n"
			"[xmin2 ymin2 zmin2]\n"
			"[xmax2 ymax2 zmax2]\n"
			"...")
		,HELP__("hr_format",p.hr_format,
			"H(R) format for 'generic' hctor mode\n"
			" - omen: OMEN style format\n"
			" - wannier90: sparse 'fake' wannier90 hr output\n"
			" - hr32r: sparse binary format, real parts only, single prec.\n"
			" - hr32c: sparse binary format, real and complex, single prec.\n"
			" - hr64r: sparse binary format, real parts only, double prec.\n"
			" - hr64c: sparse binary format, real and complex, double prec.\n"
			"for ideal MLWF data, real only single prec. ('bin32re') is sufficient\n"
			"since wannier90 outputs only 6 digits for Hamiltonian elements\n"
			"and complex parts should be negligible")
		,HELP__("hr_out",p.hr_out,
			"filename for the output Hamiltonian file")
		,HELP__("pprec",p.pprec,
			"print precision for hamiltonian elements")
		,HELP__("R",p.R,
			"manual specification of R vectors")

	));
	}
	};
};


/** generate Hamiltonian matrices
 * @param inp user input
 * @param os stream to print into
 */
void h_switch(const h_input& inp, std::ostream& os=std::cout);
/** generic Hamiltonian constructor
 * @param inp user input
 * @param B basis
 * @param Ap atomic positions in cartesian
 * @param id range of id strings
 * @param os stream to print into
 */
void h_hctor(const h_input& inp, const lm__::fMat& B, const lm__::fMat& Ap,
		const ll__::idv& id, std::ostream& os);

#endif // _H_SWITCH_

/** @}
 */
