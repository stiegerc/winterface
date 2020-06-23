// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _LL_TPARSER1_
#define _LL_TPARSER1_

#include "aux_parser.h"


struct ll_tParser1: public virtual aux_parser {
public:
	// test input
	std::string tString1;
	std::vector<std::string> tStrings1;
	double tDouble1;
	size_t tUint1;
	bool tBool1;
	std::vector<double> tDoublesF1;
	std::vector<double> tDoublesV1;
	std::vector<size_t> tUintsFNT1;
	std::vector<size_t> tUintsFT1;
	std::vector<size_t> tUintsVNT1;
	std::vector<size_t> tUintsVT1;
	std::vector<bool> tBoolsF1;
	std::vector<bool> tBoolsV1;
	lm__::fMat tMat1;
	lm__::fMat tMatF1;

public:
	// parse key struct
	struct parseKey_: public virtual aux_parser::parseKey_ {
	// constructor
	inline explicit parseKey_(ll_tParser1& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator& i, std::sregex_token_iterator& e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			case "tString1"_h: PARSE__(p.tString1); return;
			case "tStrings1"_h: PARSE__(p.tStrings1); return;
			case "tDouble1"_h: PARSE__(p.tDouble1); return;
			case "tUint1"_h: PARSE__(p.tUint1); return;
			case "tBool1"_h: PARSE__(p.tBool1); return;
			case "tDoublesF1"_h: PARSE__(p.tDoublesF1,DIM__); return;
			case "tDoublesV1"_h: PARSE__(p.tDoublesV1); return;
			case "tUintsFNT1"_h: PARSE__(p.tUintsFNT1,false,DIM__); return;
			case "tUintsFT1"_h: PARSE__(p.tUintsFT1,true,DIM__); return;
			case "tUintsVNT1"_h: PARSE__(p.tUintsVNT1,false,NPOS__); return;
			case "tUintsVT1"_h: PARSE__(p.tUintsVT1,true,NPOS__); return;
			case "tBoolsF1"_h: PARSE__(p.tBoolsF1,DIM__); return;
			case "tBoolsV1"_h: PARSE__(p.tBoolsV1); return;
			case "tMat1"_h: PARSE__(p.tMat1); return;
			case "tMatF1"_h: PARSE__(p.tMatF1); return;
		}
	}
	};

public:
	// print help struct
	struct printHelp_: public virtual aux_parser::printHelp_ {
	// constructor
	inline explicit printHelp_(const ll_tParser1& p, std::ostream& os):
			aux_parser::printHelp_(p,os) {
		printHelpTuple_(os,std::make_tuple(
		
		TOPIC__("TEST ENTRIES")
		,HELP__("tString1",p.tString1,
			"sample string")
		,HELP__("tStrings1",p.tStrings1,
			"sample array of strings")
		,HELP__("tDouble1",p.tDouble1,
			"sample double")
		,HELP__("tUint1",p.tUint1,
			"sample unsigned int")
		,HELP__("tBool1",p.tBool1,
			"should be discarded")
		,HELP__("tDoublesF1",p.tDoublesF1,
			"sample array of doubles")
		,HELP__("tDoublesV1",p.tDoublesV1,
			"sample array of doubles")
		,HELP__("tUintsFNT1",p.tUintsFNT1,DIM__,
			"sample array of unsigned ints")
		,HELP__("tUintsFT1",p.tUintsFT1,DIM__,
			"sample array of unsigned ints")
		,HELP__("tUintsVNT1",p.tUintsVNT1,NPOS__,
			"sample array of unsigned ints")
		,HELP__("tUintsVT1",p.tUintsVT1,NPOS__,
			"sample array of unsigned ints")
		,HELP__("tBoolsF1",p.tBoolsF1,DIM__,
			"sample array of bools")
		,HELP__("tBoolsV1",p.tBoolsV1,
			"sample array of bools")
		,HELP__("tMat1",p.tMat1,
			"sample matrix")
		,HELP__("tMatF1",p.tMatF1,
			"sample matrix")
		));
	}
	};
};

#endif // _LL_TPARSER1_
