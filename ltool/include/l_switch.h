// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup ltool
 * @{
 */

#ifndef _L_SWITCH_
#define _L_SWITCH_

#include "aux_parser.h"
#include "ll_defs.h"
#include "ll_types.h"


/** input for lattice manipulations
 * As a child of aux_parser, it redefines its own versions of
 * parseKey_ and printHelp_.
 */
struct l_input: public virtual aux_parser {
public:
	/** @name filenames
	 */
	std::string pscin = POSCAR__;			//!< input POSCAR file
	std::string wbh = "";				//!< wbh to extract unit cell from
	std::string pscout = "out.psc";			//!< output POSCAR file


	/** @name parameters
	 */
	ll__::rv r = ll__::rv(DIM__,false);		//!< restriction vector
	std::vector<double> vac = {};			//!< replacement vacuum length
	lm__::fMat C = lm__::eye<lm__::fMat>(DIM__);	//!< expansion matrix
	lm__::fMat R = lm__::zeros<lm__::fMat>(DIM__,1);//!< grid of R vectors
	double bond_factor = 0.0;			//!< bond factor for adding bond centers
	lm__::fMat ROT = lm__::eye<lm__::fMat>(DIM__);	//!< rotation matrix
	double phi_x = 0.0;				//!< rotation angle around x axis
	double phi_y = 0.0;				//!< rotation angle around y axis
	double phi_z = 0.0;				//!< rotation angle around z axis


	/** @name switches
	 */
	bool direct = true;			//!< switch for direct coordinates
	bool strip = false;			//!< switch for stripping indices


public:
	//! parseKey_ redefinition
	struct parseKey_: public virtual aux_parser::parseKey_ {

	/** constructor defining keys
	 */
	inline explicit parseKey_(l_input& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator& i, std::sregex_token_iterator& e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			// filenames
			case "pscin"_h: PARSE__(p.pscin); return;
			case "wbh"_h: PARSE__(p.wbh); return;
			case "pscout"_h: PARSE__(p.pscout); return;

			// parameters
			case "r"_h: PARSE__(p.r); return;
			case "vac"_h: PARSE__(p.vac); return;
			case "C"_h: PARSE__(p.C); return;
			case "R"_h: PARSE__(p.R); return;
			case "bond_factor"_h: PARSE__(p.bond_factor); return;
			case "ROT"_h: PARSE__(p.ROT); return;
			case "phi_x"_h: PARSE__(p.phi_x); return;
			case "phi_y"_h: PARSE__(p.phi_y); return;
			case "phi_z"_h: PARSE__(p.phi_z); return;

			// switches
			case "direct"_h: PARSE__(p.direct); return;
			case "strip"_h: PARSE__(p.strip); return;
		}
	}
	};


public:
	//! printHelp_ redefinition
	struct printHelp_: public virtual aux_parser::printHelp_ {

	/** constructor defining help messages
	 */
	inline explicit printHelp_(const l_input& p, std::ostream& os):
			aux_parser::printHelp_(p,os) {

	printHelpTuple_(os,std::make_tuple(
		TOPIC__("-l: FILENAMES")
		,HELP__("pscin",p.pscin,
			"POSCAR input file")
		,HELP__("wbh",p.wbh,
			"alternative to pscin, allows to use the cell stored in the wbh")
		,HELP__("pscout",p.pscout,
			"POSCAR output file")
		
		,TOPIC__("-l: PARAMETERS")
		,HELP__("r",p.r,
			"restriction vector as to which dimensions should be\n"
			"considered frozen")
		,HELP__("vac",p.vac,
			"vacuum vector, size 1 or cell.dim()")
		,HELP__("C",p.C,
			"expansion matrix")
		,HELP__("R",p.R,
			"R vectors as to which images to include")
		,HELP__("bond_factor",p.bond_factor,
			"bond factor to use when adding bond centers")
		,HELP__("ROT",p.ROT,
			"direct specification of rotation matrix")
		,HELP__("phi_x",p.phi_x,
			"angle around x-axis in degrees, applied after ROT")
		,HELP__("phi_y",p.phi_y,
			"angle around y-axis in degrees, applied after ROT")
		,HELP__("phi_z",p.phi_z,
			"angle around z-axis in degrees, applied after ROT")
		
		,TOPIC__("-l: SWITCHES")
		,HELP__("direct",p.direct,
			"switch to output direct coordinates or cartesian")
		,HELP__("strip",p.strip,
			"switch to strip indices off id strings")
	));
	}
	};
};


/** manipulate lattices
 * @param inp user input
 * @param os stream to print into
 */
void l_switch(const l_input& inp, std::ostream& os=std::cout);

#endif // _L_SWITCH_

/** @}
 */
