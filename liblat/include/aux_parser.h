// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _AUX_PARSER_
#define _AUX_PARSER_

#include "lm_defs.h"
#include "lm_types.h"
#include "lm_tMat.h"
#include "lm_fn.h"
#include "aux_io.h"
#include <cstdint>
#include <string>
#include <regex>
#include <vector>
#include <fstream>
#include <set>
#include <iostream>
#include <bitset>


#ifndef ICOM__
#define ICOM__ "!#%"		//!< comment symbols
#endif

#ifndef RGX__
#define RGX__ "[\\s,;]+"	//!< regex for tokenizing lines
#endif

#ifndef ALLOC__
#define ALLOC__ 10		//!< default allocation for vectors of indeterminant size
#endif

#ifndef DIM__
#define DIM__ 3			//!< dimension of space
#endif

#ifndef PRINTBIT__
#define PRINTBIT__ 1		//!< bitmask bit for printing to stream
#endif
#ifndef VERBOBIT__
#define VERBOBIT__ 2		//!< bitmask for verbose printing to stream
#endif
#ifndef WRITEBIT__
#define WRITEBIT__ 4		//!< bitmask bit for writing files to disk
#endif
#ifndef DEBUGBIT__
#define DEBUGBIT__ 8		//!< bitmask bit for writing debug files to disk
#endif
#ifndef MD5BIT__
#define MD5BIT__   16		//!< bitmask bit for printing md5 sums
#endif



//! parse entry macro for 1 input
#define PARSE_1__(i1) parse(i1,file,lcnt,++i,e);
//! parse entry macro for 2 inputs
#define PARSE_2__(i1,i2) parse(i1,i2,file,lcnt,++i,e);
//! parse entry macro for 3 inputs
#define PARSE_3__(i1,i2,i3) parse(i1,i2,i3,file,lcnt,++i,e);
//! parse entry macro selector
#define PARSE_SELECT__(_1,_2,_3,NAME,...) NAME
//! parse entry variadic macro
#define PARSE__(...) PARSE_SELECT__(__VA_ARGS__, \
	PARSE_3__, PARSE_2__, PARSE_1__)(__VA_ARGS__)


//! help entry macro for 1 entry
#define HELP_1__(s,i1,text) help_<decltype(i1)>{s,text,i1,0,0}
//! help entry macro for 2 entries
#define HELP_2__(s,i1,i2,text) help_<decltype(i1)>{s,text,i1,i2,0}
//! help entry macro for 3 entries
#define HELP_3__(s,i1,i2,i3,text) help_<decltype(i1)>{s,text,i1,i2,i3}
//! help entry macro selector
#define HELP_SELECT__(_1,_2,_3,_4,_5,NAME,...) NAME
//! help entry variadic macro
#define HELP__(...) HELP_SELECT__(__VA_ARGS__, \
	HELP_3__, HELP_2__, HELP_1__)(__VA_ARGS__)


//! topic macro
#define TOPIC__(text) "\n" RED__ text RESET__



namespace aux {
	//! FNV hash from a cstring
	constexpr uint32_t fnvHash(const char* str) {
		uint32_t res = 2166136261u;
		while (*str) {
			res ^= *str++;
			res *= 16777619u;
		}
		return res;
	}
	//! literal operator for FNV hash
	constexpr uint32_t operator"" _h(const char* str, size_t) {
		return fnvHash(str);
	}
}



namespace aux {
namespace parse {

	//! id strings for types aux_parser can parse
	template<class T> inline constexpr const char* idstr() noexcept {return "UNKNOWN";}

	//! idstr specialization for std::string
	template<> inline constexpr const char*
		idstr<std::string>() noexcept {return "STRING";}
	//! idstr specialization for std::vector<std::string>
	template<> inline constexpr const char*
		idstr<std::vector<std::string>>() noexcept {return "STRINGS";}
	//! idstr specialization for double
	template<> inline constexpr const char*
		idstr<double>() noexcept {return "DOUBLE";}
	//! idstr specialization for std::vector<double>
	template<> inline constexpr const char*
		idstr<std::vector<double>>() noexcept {return "DOUBLES";}
	//! idstr specialization for size_t
	template<> inline constexpr const char*
		idstr<size_t>() noexcept {return "UINT";}
	//! idstr specialization for std::vector<size_t>
	template<> inline constexpr const char*
		idstr<std::vector<size_t>>() noexcept {return "UINTS";}
	//! idstr specialization for bool
	template<> inline constexpr const char*
		idstr<bool>() noexcept {return "BOOL";}
	//! idstr specialization for std::vector<bool>
	template<> inline constexpr const char*
		idstr<std::vector<bool>>() noexcept {return "BOOLS";}
	//! idstr specialization for fMat
	template<> inline constexpr const char*
		idstr<lm__::fMat>() noexcept {return "FMAT";}
}
}



namespace aux {

	/** Main parse file function.
	 * This function is used to parse files and write the results into input structs
	 * inheriting from aux_parser \n
	 * Template arguments are: \n
	 * - CT: the input struct
	 * @param fileName the name of the file to parse from
	 * @param lcnt linecounter, counts the number of lines parsed in the file
	 * @param kcnt keycounter, counts the number of keys parsed in the file
	 * @param os stream to print into
	 */
	template<class CT>
	inline CT parseFile(const std::string& fileName, size_t& lcnt, size_t& kcnt,
			std::ostream& os=std::cout) {
		
		// open file
		auto file = aux::openFile<std::ifstream>(fileName);

		// result container
		CT res;

		// set of parsed keys
		std::set<uint32_t> pk;

		// regex
		const std::regex rgx(RGX__);
		
		// read line by line and tokenize
		lcnt=0; // line counter
		try {
			std::string line;
			while (std::getline(file,line)) {
				++lcnt;
				
				// resize line to first occurrence of comment symbol
				line = line.substr(0,line.find_first_of(ICOM__));
				
				std::sregex_token_iterator i(line.begin(),line.end(),rgx,-1), e;
				if (i==e) continue;
				if (!i->length()) ++i; // line starts with delimiter
				if (i==e) continue;

				// hash key and check for duplicate
				const uint32_t key = aux::fnvHash(i->str().c_str());
				if (pk.find(key)!=pk.end())
					throw(std::invalid_argument(
						"duplicate key \'"+i->str()+"\'"));
				
				// parse key
				typename CT::parseKey_ prsk(res,key,file,lcnt,i,e);
				if (!prsk.rcnt())
					throw(std::invalid_argument("invalid key \'"+i->str()+"\'"));

				// insert parsed key
				pk.insert(key);
			}
		} catch(const std::exception& e) {
			throw(std::invalid_argument("file: "+fileName+", line "+
				std::to_string(lcnt)+", "+e.what()));
		}
		
		// check whole file was parsed and close file
		assert(file.eof());
		file.close();

		kcnt = pk.size(); // key counter
		if (res.verbosity) {
			os << "file \'" << fileName << "\' successfully parsed\n";
			os << kcnt << " keys parsed on " << lcnt << " lines\n\n";
		}

		return res;
	}
	/** parse file with internal counters.
	 * Template arguments are: \n
	 * - CT: the input struct
	 * @param fileName name of the file
	 * @param os stream to print into
	 */
	template<class CT>
	inline CT parseFile(const std::string& fileName, std::ostream& os=std::cout) {
		size_t d1, d2;
		return parseFile<CT>(fileName,d1,d2,os);
	}

	/** print help message function
	 * Template arguments are: \n
	 * - CT: the input struct
	 * @param inp user input struct
	 * @param os stream to print into
	 */
	template<class CT>
	inline void printHelp(const CT& inp, std::ostream& os=std::cout) {
		typename CT::printHelp_ prnt(inp,os);
	}
}



/** Basic parser struct to inherit from.
 * Inheriting from this struct enables initialization of its member variables
 * by parsing via the parseFile function.
 * An entire inheritance hierarchy is possible, allowing for splitting input
 * parameters into multiple structs. The sub structs included in aux_parser,
 * parseKey_ and printHelp_ must be redefined in each child of aux_parser and
 * the redefinitions must define their own versions. Multiple examples of this
 * process can be found in the code base. \n
 * The reasoning behind this design is that we want to group input parameters
 * for functions who depend on user input and would otherwise have a large
 * number of arguments. Additionally such structs don't have to be initialized
 * by parsing and can be used in the code like any other struct, thus the
 * parsing process does not get in the way when not needed.
 */
struct aux_parser {
public:
	/** @name basic members
	 */
	size_t verbosity = (PRINTBIT__ | WRITEBIT__);	//!< level of verbosity
	std::string prefix = "./";			//!< prefix to be added in front of output files
	
public:
	/** Parse key struct. This struct defines how keys are parsed, that is it provides the
	 * keywords in the input script and into which variable the result should be stored.
	 * Each child of aux_parser must provide their own version of this struct which must
	 * inherit from this class. The parsing process is defined in the constructor,
	 * who will then call all the other constructors in the inheritance hierarchy from top
	 * to bottom, thus parsing all the keys. Types of keys supported can be extended by
	 * writing a parse_ for a new type.
	 */
	struct parseKey_ {
	public:
		/** @name constructor
		 */
		/** Constructor by parsing from file.
		 * @param p the host struct of this struct
		 * @param key key to be parsed
		 * @param file stream to parse from
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		explicit parseKey_(aux_parser& p, const uint32_t key,
			std::ifstream& file, size_t& lcnt,
			std::sregex_token_iterator i, std::sregex_token_iterator e);

		//! key read counter
		bool rcnt() const noexcept { return rcnt_; }

		//! parse wrapper with rcnt increment for one argument
		template<class T> inline void parse(T& inp, std::ifstream& file, size_t& lcnt,
			std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			parse_(inp,file,lcnt,i,e); ++rcnt_;
		}
		//! parse wrapper with rcnt increment for two arguments
		template<class T> inline void parse(T& inp, const size_t n,
			std::ifstream& file, size_t& lcnt,
			std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			parse_(inp,n,file,lcnt,i,e); ++rcnt_;
		}
		//! parse wrapper with rcnt increment for three arguments
		template<class T> inline void parse(T& inp, const bool treat, const size_t n,
			std::ifstream& file, size_t& lcnt,
			std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			parse_(inp,treat,n,file,lcnt,i,e); ++rcnt_;
		}

	private:
		/** @name member variables
		 */
		size_t rcnt_ = 0;	//!< key read flag

	protected:
		/** @name helpers
		 */
		/** Generic exception thrower. This function advances the regex iterator past the key
		 * and the optional '=' to where the data is.
		 */
		static inline void bitch_(std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			if (i==e) // complain about empty tokens after key
				throw(std::invalid_argument("argument expected"));
			if (*i=="=") {
				++i;
				if (i==e) // complain about empty tokens after '='
				throw(std::invalid_argument("argument expected"));
			}
		}

		/** @name parse_ functions for each type.
		 */
		/** parse_ function for std::string.
		 * @param inp the std::string to store the result in
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				std::string& inp, std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			bitch_(i,e);
			assert(file||true||lcnt);

			// check for "" or ''
			inp = *i;
			if (inp.length()>2 && ((inp.front()=='\"' && inp.back()=='\"')||
					       (inp.front()=='\'' && inp.back()=='\''))) {
				inp.pop_back();
				inp.erase(inp.begin());
			}
		}
		/** parse_ function for std::vector<std::string>.
		 * @param inp the std::vector<std::string> to store the result in
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				std::vector<std::string>& inp, std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			bitch_(i,e);
			
			inp.clear();
			while (i!=e) {
				if (!i->str().empty()) {
					std::string buff;
					parse_(buff,file,lcnt,i,e);
					inp.push_back(std::move(buff));
				}
				++i;
			}
		}
		/** parse_ function for double.
		 * @param inp the double to store the result in
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				double& inp, std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			bitch_(i,e);
			assert(file||true||lcnt);
			inp = std::stod(*i);
		}
		/** parse_ function for size_t.
		 * @param inp the size_t to store the result in
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				size_t& inp, std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			bitch_(i,e);
			assert(file||true||lcnt);
			if (i->str().front()=='-')
				throw(std::invalid_argument("index must be >=0"));
			inp = std::stoul(*i);
		}
		/** parse_ function for bool.
		 * @param inp the bool to store the result in
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				bool& inp, std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			bitch_(i,e);
			assert(file||true||lcnt);
			if (i->str().front()=='t' || i->str().front()=='T') {inp = true; return;}
			if (i->str().front()=='f' || i->str().front()=='F') {inp = false; return;}
			inp =  std::stoi(*i);
		}
		/** parse_ function for std::vector<double> of fixed size.
		 * @param inp the std::vector<double> to store the result in
		 * @param n the number doubles to read
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				std::vector<double>& inp, const size_t n,
				std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			bitch_(i,e);
			assert(file||true||lcnt);

			inp.clear();
			if (n<NPOS__) {
				inp.reserve(n);
				while (i!=e && inp.size()<=n)
					inp.push_back(std::stod((i++)->str()));
			} else {
				inp.reserve(ALLOC__);
				while (i!=e) {
					if (inp.size()==inp.capacity()) inp.reserve(inp.size()+ALLOC__);
					inp.push_back(std::stod((i++)->str()));
				}
				inp.shrink_to_fit();
			}
		}
		/** parse_ function for std::vector<double>.
		 * @param inp the std::vector<double> to store the result in
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				std::vector<double>& inp,
				std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			parse_(inp,NPOS__,file,lcnt,i,e);
		}
		/** parse_ function for std::vector<size_t> of fixed size.
		 * @param inp the std::vector<size_t> to store the result in
		 * @param treat switch to sort(true) the vector or leave it(false)
		 * @param n the number doubles to read
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				std::vector<size_t>& inp, const bool treat, const size_t n,
				std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			bitch_(i,e);
			assert(file||true||lcnt);
			
			inp.clear();
			auto extRes = [&inp](const std::string& tok) -> void {
				const size_t mpos = tok.find_first_of("-:");
					
				if (!mpos) throw(std::invalid_argument("indices must be >=0"));
				if (mpos>tok.size()) 				// one number
					inp.push_back(std::stoul(tok));
				else {						// range	
					size_t l = std::stoul(tok.substr(0,mpos));
					size_t u = std::stoul(tok.substr(mpos+1));
					if (u<=l) std::swap(l,u);
					
					inp.reserve(inp.size()+u-l);
					for (size_t j=l; j<=u; ++j)
						inp.push_back(j);
				}
			};	
			if (n<NPOS__) {
				inp.reserve(n);
				while (i!=e && inp.size()<=n) extRes((i++)->str());
			} else {
				inp.reserve(ALLOC__);
				while (i!=e) {
					if (inp.size()==inp.capacity()) inp.reserve(inp.size()+ALLOC__);
					extRes((i++)->str());
				}
				inp.shrink_to_fit();
			}

			// sort and remove duplicates
			if (treat) {
				std::sort(inp.begin(),inp.end());
				inp.erase(std::unique(inp.begin(),inp.end()),inp.end());
			}
		}
		/** parse_ function for std::vector<size_t> of fixed size.
		 * @param inp the std::vector<size_t> to store the result in
		 * @param n the number doubles to read
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				std::vector<size_t>& inp, const size_t n,
				std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			parse_(inp,true,n,file,lcnt,i,e);
		}
		/** parse_ function for std::vector<size_t>.
		 * @param inp the std::vector<size_t> to store the result in
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				std::vector<size_t>& inp,
				std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			parse_(inp,true,NPOS__,file,lcnt,i,e);
		}
		/** parse_ function for std::vector<bool> of fixed size.
		 * @param inp the std::vector<bool> to store the result in
		 * @param n the number doubles to read
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				std::vector<bool>& inp, const size_t n,
				std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			bitch_(i,e);
			assert(file||true||lcnt);

			inp.clear();
			auto extRes = [&inp](const std::string& tok) -> void {
				if (tok=="t" || tok=="T") {
					inp.push_back(true);
					return;
				}
				if (tok=="f" || tok=="F") {
					inp.push_back(false);
					return;
				}
				inp.push_back(std::stoi(tok));
			};

			if (n<NPOS__) {
				inp.reserve(n);
				while (i!=e && inp.size()<=n) extRes((i++)->str());
			} else {
				inp.reserve(ALLOC__);
				while (i!=e) {
					if (inp.size()==inp.capacity()) inp.reserve(inp.size()+ALLOC__);
					extRes((i++)->str()); 
				}
				inp.shrink_to_fit();
			}
		}
		/** parse_ function for std::vector<bool>.
		 * @param inp the std::vector<bool> to store the result in
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				std::vector<bool>& inp,
				std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			parse_(inp,NPOS__,file,lcnt,i,e);
		}
		/** parse_ function for fMat.
		 * @param inp the fMat to store the result in
		 * @param file the file from which parsing is being done
		 * @param lcnt linecounter
		 * @param i regex iterator pointing to the first token of the line being parsed
		 * @param e regex iterator past the end of the tokens in the line being parsed
		 */
		inline static void parse_(
				lm__::fMat& inp,
				std::ifstream& file, size_t& lcnt,
				std::sregex_token_iterator& i, const std::sregex_token_iterator& e) {
			
			// read mat from file
			if (std::distance(i,e)) {
				std::string fileName_; parse_(fileName_,file,lcnt,i,e);
				inp = lm__::fMat(fileName_).T();
				return;
			}

			// tokenize and read doubles lambda
			std::regex rgx(RGX__);
			auto tokrd = [&lcnt,&file,&rgx]() -> std::vector<double> {
				std::string line;
				
				if (!std::getline(file,line)) return {};
				++lcnt;

				// resize line to first occurrence of comment symbol
				line = line.substr(0,line.find_first_of(ICOM__));
				
				// tokenize line
				std::sregex_token_iterator i(line.begin(),line.end(),rgx,-1), e;
				if (i==e) return {}; // return empty when encountering empty line
				if (!i->length()) ++i;
				if (i==e) return {}; // return empty when encountering empty line

				std::vector<double> res; parse_(res,NPOS__,file,lcnt,i,e);
				return res;
			};

			// parse first line
			auto ccol = tokrd();
			if (!ccol.size()) {
				inp = lm__::fMat(0,0); return;
			}

			// create result
			inp = lm__::fMat(ccol); inp.reserve(DIM__);

			// parse till empty line is encountered
			while ((ccol=tokrd()).size()) {
				if (ccol.size()!=inp.M())
					throw(std::invalid_argument("found "+std::to_string(ccol.size())+
						" entries, expected "+std::to_string(DIM__)));
				if (inp.ccap()==inp.N()) inp.reserve(inp.N()+DIM__);
				inp.push_back(lm__::fMat(ccol));
			}
			inp.shrink_to_fit();
		}
	};

public:
	/** Print help struct. This struct defines how help messages for each key should
	 * be printed. Each child of aux_parser must provide their own version of this struct
	 * which must inherit from this class. The printing process is defined in the constructor,
	 * who will then call all the other constructors in the inheritance hierarchy from top
	 * to bottom, thus printing all the messages. Additionally, the current contents of a child
	 * of aux_parser will be printed. Thus when not initialized, the default values will be
	 * printed. The recognized types can be extended by writing additional print_ functions.
	 */
	struct printHelp_ {
	public:
		/** @name constructors
		 */
		/** Constructor from struct.
		 * @param p struct holding variables to parse
		 * @param os stream to print into
		 */
		explicit printHelp_(const aux_parser& p, std::ostream& os);

	private:
		/** @name printer
		 */
		//! tuple printer
		template<class Tuple, std::size_t N>
		struct tuplePrinter_ {
			static void print(std::ostream& os, const Tuple& t) {
				tuplePrinter_<Tuple, N-1>::print(os,t);
				os << "\n" << std::get<N-1>(t);
			}
		};
		//! tuple printer root
		template<class Tuple>
		struct tuplePrinter_<Tuple, 1> {
			static void print(std::ostream& os, const Tuple& t) {
				os << "\n" << std::get<0>(t);
			}
		};
	
	protected:
		/** Basic help struct. This provides a printing interface.
		 */
		template<class T>
		struct help_ {
		public:
			/** @name member variables
			 */
			const char* key;	//!< the key in the file
			const char* text;	//!< the help message
			const T& def;		//!< the variable corresponding to the key
			const size_t M;		//!< optional size argument (number of rows)
			const size_t N;		//!< optional size argument (number of columns)
			
			//! print type dimensions
			inline std::string type() const noexcept {
				return aux::parse::idstr<T>()+
					(M ? M==NPOS__ ? " ...": " "+std::to_string(M): "")+
					(N ? N==NPOS__ ? "x...": "x"+std::to_string(N): "");
			}

			//! print help struct in defined format
			inline friend std::ostream& operator<<(std::ostream& os,
					const help_<T>& h) noexcept {
				return (os << GREEN__ << h.key << CYAN__ << " [" << h.type() << "]\n"
					   << RESET__ << h.text << "\n"
					   << CYAN__ << print_(h.def) << RESET__ << "\n");
			}

			//! print help format
			inline static constexpr const char* format() noexcept {
				return GREEN__ "key" CYAN__ " [TYPE COLSxROWS]\n"
				       RESET__ "purpose of the key and its input\n"
				       CYAN__  "value" RESET__ "\n";
			}
			
			//! defines how to print std::string
			inline static std::string print_(const std::string& inp) noexcept {
				return inp.empty() ? "\"\"": inp;
			}
			//! defines how to print double
			inline static std::string print_(const double inp) noexcept {
				std::stringstream sstr;
				sstr << std::scientific << inp;
				return sstr.str();
			}
			//! defines how to print size_t
			inline static std::string print_(const size_t inp) noexcept {
				return std::to_string(inp);
			}
			//! defines how to print bool
			inline static std::string print_(const bool inp) noexcept {
				return std::to_string(inp);
			}
			//! defines how to generic vectors
			template<class VT>
			inline static std::string print_(const std::vector<VT>& inp) noexcept {
				using namespace aux;
				if (inp.empty()) return "{}";
				std::stringstream sstr; sstr << inp;
				return sstr.str();
			}
			//! defines how to print std::vector<std::string>
			inline static std::string print_(const std::vector<std::string>& inp) noexcept {
				if (inp.empty()) return "{}";
				std::stringstream sstr;
				for (const auto& s: inp) sstr << print_(s);
				return sstr.str();
			}
			//! defines how to print std::vector<double>
			inline static std::string print_(const std::vector<double>& inp) noexcept {
				using namespace aux;
				if (inp.empty()) return "{}";
				std::stringstream sstr; sstr << std::scientific << inp;
				return sstr.str();
			}
			//! defines how to print fMat
			inline static std::string print_(const lm__::fMat& inp) noexcept {
				if (inp.empty()) return "{}";
				return lm__::T(inp).print();
			}
		};
	
	protected:
		//! format strings
		inline constexpr const char* helpFormat() noexcept {
			return help_<size_t>::format();
		}

		//! print help tuple
		template<class... Args>
		void printHelpTuple_(std::ostream& os, const std::tuple<Args...>& t) const noexcept {
			tuplePrinter_<decltype(t),sizeof...(Args)>::print(os,t);
		}
	};
};



//! streaming operator for aux_parser
inline std::ostream& operator<<(std::ostream& os, const aux_parser& inp) noexcept {
	aux::printHelp(inp,os); return os;
}
//! streaming operator for children of aux_parser
template<class T, std::enable_if_t<std::is_base_of<aux_parser,T>::value>* = nullptr>
inline std::ostream& operator<<(std::ostream& os, const T& inp) noexcept {
	aux::printHelp(inp,os); return os;
}



namespace aux {
	//! read file, screen out keys and dump result into provided stream
	void screenFile(std::ostream& sstr, const std::string& fileName,
			const std::vector<std::string>& keys={});
}

#endif // _AUX_PARSER_

/** @}
 */
