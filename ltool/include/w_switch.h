// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup ltool
 * @{
 */

#ifndef _W_SWITCH_
#define _W_SWITCH_

#include <iostream>
#include "aux_parser.h"
#include "ll_hbonds.h"
#include "ll_BStest.h"
#include "ll_defs.h"
#include "ll_types.h"


/** input to generate wbh
 * As a child of aux_parser, it redefines its own versions of
 * parseKey_ and printHelp_.
 */
struct w_input: public virtual aux_parser,
       		public virtual ll_wmatching_input,
		public virtual ll_BStest_input {
public:
	/** @name filenames
	 */
	std::string wbh = WBH__;			//!< output wbh
	std::string outcar = "";			//!< OUTCAR file for Fermi energy
	std::string weig = "wannier90.eig";		//!< wannier90.eig file
	

	/** @name parameters
	 */
	double Ef = nan("");				//!< Fermi energy
	fMat C = lm__::eye<fMat>(DIM__,DIM__);		//!< expansion matrix for bandstructure test


public:
	//! parseKey_ redefinition
	struct parseKey_: public virtual aux_parser::parseKey_,
			  public virtual ll_wmatching_input::parseKey_,
			  public virtual ll_BStest_input::parseKey_ {

	/** constructor defining keys
	 */
	inline explicit parseKey_(w_input& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator& i, std::sregex_token_iterator& e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e),
			ll_wmatching_input::parseKey_(p,key,file,lcnt,i,e),
			ll_BStest_input::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			// filenames
			case "wbh"_h: PARSE__(p.wbh); return;
			case "outcar"_h: PARSE__(p.outcar); return;
			case "weig"_h: PARSE__(p.weig); return;

			// parameters
			case "Ef"_h: PARSE__(p.Ef); return;
			case "C"_h: PARSE__(p.C); return;
		}
	}
	};


public:
	//! printHelp_ redefinition
	struct printHelp_: public virtual aux_parser::printHelp_,
			   public virtual ll_wmatching_input::printHelp_,
			   public virtual ll_BStest_input::printHelp_ {

	/** constructor defining help messages
	 */
	inline explicit printHelp_(const w_input& p, std::ostream& os):
			aux_parser::printHelp_(p,os),
			ll_wmatching_input::printHelp_(p,os),
			ll_BStest_input::printHelp_(p,os) {

	printHelpTuple_(os,std::make_tuple(
		TOPIC__("-f: FILENAMES")
		,HELP__("wbh",p.wbh,
			"filename for the output wbh")
		,HELP__("outcar",p.outcar,
			"VASP OUTCAR file to determine Fermi energy")
		,HELP__("weig",p.weig,
			"Wannier90 eig file to determine band edges")
		
		,TOPIC__("-f: PARAMETERS")
		,HELP__("Ef",p.Ef,
			"override for the Fermi energy")
		,HELP__("C",p.C,
			"expansion matrix to use in the bandstructure test")
	));
	}
	};
};

#endif // _W_SWITCH_

/** @}
 */
