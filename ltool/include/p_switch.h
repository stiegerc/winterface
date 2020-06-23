// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup ltool
 * @{
 */

#ifndef _P_SWITCH_
#define _P_SWITCH_

#include "aux_parser.h"
#include "ll_defs.h"
#include "ll_types.h"


/** input struct for finding primitive cells
 * As a child of aux_parser, it redefines its own versions of
 * parseKey_ and printHelp_.
 */
struct p_input: public virtual aux_parser {
public:
	/** @name filenames
	 */
	std::string pscin = POSCAR__;		//!< input POSCAR
	std::string pscout = POSCAR__"_mod";	//!< output POSCAR
	

	/** @name parameters
	 */
	bool direct = true;			//!< switch for direct coordinates
	bool strip = false;			//!< switch for stripping indices
	double tol = WTOL__;			//!< tolerance level
	ll__::rv r = ll__::rv(DIM__,false);	//!< restriction vector


public:
	//! parseKey_ redefinition
	struct parseKey_: public virtual aux_parser::parseKey_ {

	/** constructor defining keys
	 */
	inline explicit parseKey_(p_input& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator& i, std::sregex_token_iterator& e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			// filenames
			case "pscin"_h: PARSE__(p.pscin); return;
			case "pscout"_h: PARSE__(p.pscout); return;

			// parameters
			case "r"_h: PARSE__(p.r); return;
			case "tol"_h: PARSE__(p.tol); return;

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
	inline explicit printHelp_(const p_input& p, std::ostream& os):
			aux_parser::printHelp_(p,os) {

	printHelpTuple_(os,std::make_tuple(
		TOPIC__("-p: FILENAMES")
		,HELP__("pscin",p.pscin,
			"POSCAR input file")
		,HELP__("pscout",p.pscout,
			"POSCAR output file")
		
		,TOPIC__("-p: PARAMETERS")
		,HELP__("r",p.r,
			"restriction vector as to which dimensions should be\n"
			"considered frozen")
		,HELP__("tol",p.tol,
			"spacial tolerance\n")

		
		,TOPIC__("-p: SWITCHES")
		,HELP__("direct",p.direct,
			"output POSCAR files in direct coordinates")
		,HELP__("strip",p.strip,
			"strip indices from id strings in output POSCAR file")
	));
	}
	};
};


/** find primitive cells
 * @param inp user input
 * @param os stream to print into
 */
void p_switch(const p_input& inp, std::ostream& os=std::cout);

#endif // _P_SWITCH_

/** @}
 */
