// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup ltool
 * @{
 */

#ifndef _C_SWITCH_
#define _C_SWITCH_

#include <iostream>
#include "aux_parser.h"
#include "ll_hbonds.h"
#include "ll_BStest.h"
#include "ll_defs.h"
#include "ll_types.h"


/** input to convert deprecated wbh format
 * As a child of aux_parser, it redefines its own versions of
 * parseKey_ and printHelp_.
 */
struct c_input: public virtual aux_parser {
public:
	/** @name filenames
	 */
	std::string bh = "bh.wad";		//!< old wbh filename
	std::string psc = "dbg/cell.psc";	//!< poscar file for the unit cell
	std::string wbh = WBH__;		//!< new wbh filename
	std::string outcar = "";		//!< VASP OUTCAR for optional Fermi energy
	

	/** @name parameters
	 */
	double Ef = nan("");			//!< Fermi energy manual specification


public:
	//! parseKey_ redefinition
	struct parseKey_: public virtual aux_parser::parseKey_ {

	/** constructor defining keys
	 */
	inline explicit parseKey_(c_input& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator& i, std::sregex_token_iterator& e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			// filenames
			case "bh"_h: PARSE__(p.bh); return;
			case "psc"_h: PARSE__(p.psc); return;
			case "wbh"_h: PARSE__(p.wbh); return;
			case "outcar"_h: PARSE__(p.outcar); return;

			// parameters
			case "Ef"_h: PARSE__(p.Ef); return;
		}
	}
	};


public:
	//! printHelp_ redefinition
	struct printHelp_: public virtual aux_parser::printHelp_ {

	/** constructor defining help messages
	 */
	inline explicit printHelp_(const c_input& p, std::ostream& os):
			aux_parser::printHelp_(p,os) {

	printHelpTuple_(os,std::make_tuple(
		TOPIC__("-f: FILENAMES")
		,HELP__("bh",p.bh,
			"filename for the input bh")
		,HELP__("psc",p.psc,
			"filename for the cell.psc file")
		,HELP__("wbh",p.wbh,
			"filename for the output wbh")
		,HELP__("outcar",p.outcar,
			"VASP OUTCAR file to determine Fermi energy")
		
		,TOPIC__("-f: PARAMETERS")
		,HELP__("Ef",p.Ef,
			"override for the Fermi energy")
	));
	}
	};
};


/** convert old bh format to wbh format.
 * @param inp user input
 * @param os stream to print into
 */
void c_switch(const c_input& inp, std::ostream& os=std::cout);

#endif // _C_SWITCH_

/** @}
 */
