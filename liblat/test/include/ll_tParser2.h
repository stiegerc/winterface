// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _LL_TPARSER2_
#define _LL_TPARSER2_

#include "aux_parser.h"


struct ll_tParser2: public virtual aux_parser {
public:
	// test input
	std::string tString2;
	std::vector<std::string> tStrings2;
	double tDouble2;
	size_t tUint2;
	bool tBool2;
	std::vector<double> tDoublesF2;
	std::vector<double> tDoublesV2;
	std::vector<size_t> tUintsFNT2;
	std::vector<size_t> tUintsFT2;
	std::vector<size_t> tUintsVNT2;
	std::vector<size_t> tUintsVT2;
	std::vector<bool> tBoolsF2;
	std::vector<bool> tBoolsV2;
	lm__::fMat tMat2;
	lm__::fMat tMatF2;

public:
	// parse key struct
	struct parseKey_: public virtual aux_parser::parseKey_ {
	// constructor
	inline explicit parseKey_(ll_tParser2& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator& i, std::sregex_token_iterator& e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			case "tString2"_h: PARSE__(p.tString2); return;
			case "tStrings2"_h: PARSE__(p.tStrings2); return;
			case "tDouble2"_h: PARSE__(p.tDouble2); return;
			case "tUint2"_h: PARSE__(p.tUint2); return;
			case "tBool2"_h: PARSE__(p.tBool2); return;
			case "tDoublesF2"_h: PARSE__(p.tDoublesF2,DIM__); return;
			case "tDoublesV2"_h: PARSE__(p.tDoublesV2); return;
			case "tUintsFNT2"_h: PARSE__(p.tUintsFNT2,false,DIM__); return;
			case "tUintsFT2"_h: PARSE__(p.tUintsFT2,true,DIM__); return;
			case "tUintsVNT2"_h: PARSE__(p.tUintsVNT2,false,NPOS__); return;
			case "tUintsVT2"_h: PARSE__(p.tUintsVT2,true,NPOS__); return;
			case "tBoolsF2"_h: PARSE__(p.tBoolsF2,DIM__); return;
			case "tBoolsV2"_h: PARSE__(p.tBoolsV2); return;
			case "tMat2"_h: PARSE__(p.tMat2); return;
			case "tMatF2"_h: PARSE__(p.tMatF2); return;
		}
	}
	};

public:
	// print help struct
	struct printHelp_: public virtual aux_parser::printHelp_ {
	// constructor
	inline explicit printHelp_(const ll_tParser2& p, std::ostream& os):
			aux_parser::printHelp_(p,os) {
		printHelpTuple_(os,std::make_tuple(
		
		TOPIC__("TEST ENTRIES")
		,HELP__("tString2",p.tString2,
			"sample string")
		,HELP__("tStrings2",p.tStrings2,
			"sample array of strings")
		,HELP__("tDouble2",p.tDouble2,
			"sample double")
		,HELP__("tUint2",p.tUint2,
			"sample unsigned int")
		,HELP__("tBool2",p.tBool2,
			"should be discarded")
		,HELP__("tDoublesF2",p.tDoublesF2,
			"sample array of doubles")
		,HELP__("tDoublesV2",p.tDoublesV2,
			"sample array of doubles")
		,HELP__("tUintsFNT2",p.tUintsFNT2,DIM__,
			"sample array of unsigned ints")
		,HELP__("tUintsFT2",p.tUintsFT2,DIM__,
			"sample array of unsigned ints")
		,HELP__("tUintsVNT2",p.tUintsVNT2,NPOS__,
			"sample array of unsigned ints")
		,HELP__("tUintsVT2",p.tUintsVT2,NPOS__,
			"sample array of unsigned ints")
		,HELP__("tBoolsF2",p.tBoolsF2,DIM__,
			"sample array of bools")
		,HELP__("tBoolsV2",p.tBoolsV2,
			"sample array of bools")
		,HELP__("tMat2",p.tMat2,
			"sample matrix")
		,HELP__("tMatF2",p.tMatF2,
			"sample matrix")
		));
	}
	};
};

#endif // _LL_TPARSER2_
