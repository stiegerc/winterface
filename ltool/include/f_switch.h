// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup ltool
 * @{
 */

#ifndef _F_SWITCH_
#define _F_SWITCH_

#include <iostream>
#include "aux_parser.h"
#include "ll_defs.h"
#include "ll_types.h"


/** input for converting things to POSCAR files
 * As a child of aux_parser, it redefines its own versions of
 * parseKey_ and printHelp_.
 */
struct f_input: public virtual aux_parser {
public:
	/** @name filenames
	 */
	std::string wout = "";			//!< wannier90.wout file
	std::string lattice_dat = "";		//!< OMEN style lattice file
	std::string pscout = "out.psc";		//!< output POSCAR file


	/** @name switches
	 */
	bool direct = true;			//!< switch to use direct coordinates
	bool strip = false;			//!< switch to strip indices


public:
	//! parseKey_ redefinition
	struct parseKey_: public virtual aux_parser::parseKey_ {

	/** constructor defining keys
	 */
	inline explicit parseKey_(f_input& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator& i, std::sregex_token_iterator& e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			// filenames
			case "wout"_h: PARSE__(p.wout); return;
			case "lattice_dat"_h: PARSE__(p.lattice_dat); return;
			case "pscout"_h: PARSE__(p.pscout); return;

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
	inline explicit printHelp_(const f_input& p, std::ostream& os):
			aux_parser::printHelp_(p,os) {

	printHelpTuple_(os,std::make_tuple(
		TOPIC__("-f: FILENAMES")
		,HELP__("wout",p.wout,
			"wannier90 wout file")
		,HELP__("lattice_dat",p.lattice_dat,
			"OMEN lattice_dat file")
		,HELP__("pscout",p.pscout,
			"POSCAR output file")
		
		,TOPIC__("-f: SWITCHES")
		,HELP__("direct",p.direct,
			"switch to output direct coordinates or cartesian")
		,HELP__("strip",p.strip,
			"switch to strip indices off id strings")
	));
	}
	};
};


/** convert files to POSCAR
 * @param inp user input
 * @param os stream to print into
 */
void f_switch(const f_input& inp, std::ostream& os=std::cout);

#endif // _F_SWITCH_

/** @}
 */
