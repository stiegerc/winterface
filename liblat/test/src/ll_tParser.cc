// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "ll_tParser.h"
#include "aux_fnvHash.h"
#include "ll_defs.h"

// constructor
ll_tParser::ll_tParser(const std::string& fileName): aux_parser(fileName) { parseFile_(); }

// parse function
void ll_tParser::parseKey_(const uint32_t key, std::sregex_token_iterator& i, std::sregex_token_iterator& e) {

	switch(key) {
		PARSE_ENTRY__(tString);
		PARSE_ENTRY__(tStrings);
		PARSE_ENTRY__(tDouble);
		PARSE_ENTRY__(tUint);
		PARSE_ENTRY__(tBool);
		PARSE_ENTRY__(tDoublesF,DIM__);
		PARSE_ENTRY__(tDoublesV);
		PARSE_ENTRY__(tUintsFNT,false,DIM__);
		PARSE_ENTRY__(tUintsFT,true,DIM__);
		PARSE_ENTRY__(tUintsVNT,false,NPOS__);
		PARSE_ENTRY__(tUintsVT,true,NPOS__);
		PARSE_ENTRY__(tBoolsF,DIM__);
		PARSE_ENTRY__(tBoolsV);
		PARSE_ENTRY__(tMat);
		PARSE_ENTRY__(tMatF);
	}

	// call base class parse key
	aux_parser::parseKey_(key,i,e);
}
