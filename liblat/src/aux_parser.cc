// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "aux_parser.h"
#include <fstream>
#include <cassert>
#include <iostream>


aux_parser::parseKey_::parseKey_(aux_parser& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator i, std::sregex_token_iterator e) {
	using namespace aux;

	switch(key) {
		case "verbosity"_h: PARSE__(p.verbosity); return;
		case "prefix"_h: PARSE__(p.prefix); return;
	};
}
aux_parser::printHelp_::printHelp_(const aux_parser& p, std::ostream& os) {
	printHelpTuple_(os,std::make_tuple(
		"legal comment symbols: " CYAN__ ICOM__ RESET__
		,TOPIC__("HELP FORMAT")
		,help_<size_t>::format(),"\n"

		TOPIC__("BASIC PARSER")
		,HELP__("verbosity",p.verbosity,
			("verbosity bitmask:\n"+
			std::bitset<4>(0).to_string()+" = "+std::to_string(0)+
				": no output\n"+
			std::bitset<4>(PRINTBIT__).to_string()+" = "+std::to_string(PRINTBIT__)+
				": text output\n"+
			std::bitset<4>(VERBOBIT__).to_string()+" = "+std::to_string(VERBOBIT__)+
				": verbose text output\n"+
			std::bitset<4>(WRITEBIT__).to_string()+" = "+std::to_string(WRITEBIT__)+
				": additional files are written\n"+
			std::bitset<4>(DEBUGBIT__).to_string()+" = "+std::to_string(DEBUGBIT__)+
				": even more files are written\n"+
			" sum any of the above for combined output\n"+
			" e.g. "+std::to_string(PRINTBIT__ | WRITEBIT__)+
			" for text and file output").c_str())
		,HELP__("prefix",p.prefix,
			"prefix to be appended in front of output files")
				
	));
}



namespace aux {
	// screen file for keys and return file as screened stringstream
	void screenFile(std::ostream& sstr, const std::string& fileName,
			const std::vector<std::string>& keys) {

		// open file
		auto file = aux::openFile<std::ifstream>(fileName);

		const std::regex rgx(RGX__);

		// read line by line, tokenize and search for key
		std::string line,line_;
		while (std::getline(file,line)) {
		
			// resize line to first occurrence of comment symbol
			line_ = line.substr(0,line.find_first_of(ICOM__));
		
			// tokenize line
			std::sregex_token_iterator i(line_.begin(),line_.end(),rgx,-1), e;
			if (i==e) { sstr << line << "\n"; continue; }
			if (!i->length()) ++i; // line starts with delimiter
			if (i==e) { sstr << line << "\n"; continue; }
			
			// search in keys for first token
			if (std::find(keys.cbegin(),keys.cend(),i->str())!=keys.cend()) {
				// key found, skip lines...
				// check if we have a matrix
				++i; if (i!=e && *i=="=") ++i;
				if (i==e) // no argument after key or '=' -> matrix
					do std::getline(file,line);
					while (!line.empty() && !file.eof());
			} else {
				// key not found, dump into stream...
				sstr << line << "\n";
				
				// check if we have a matrix
				++i; if (i!=e && *i=="=") ++i;
				if (i==e) // no argument after key or '=' -> matrix
					do {
						std::getline(file,line);
						if (file.eof()) break;
						sstr << line << "\n";
					} while (!line.empty());
			}
		}
	}
}
