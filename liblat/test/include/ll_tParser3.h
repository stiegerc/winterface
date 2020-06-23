// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _LL_TPARSER3_
#define _LL_TPARSER3_

#include "ll_tParser1.h"
#include "ll_tParser2.h"


struct ll_tParser3: public virtual ll_tParser1, public virtual ll_tParser2 {
public:
	// nothing added

public:
	// parse key struct
	struct parseKey_: public virtual ll_tParser1::parseKey_,
			  public virtual ll_tParser2::parseKey_	{
	// constructor
	inline explicit parseKey_(ll_tParser3& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator& i, std::sregex_token_iterator& e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e),
			ll_tParser1::parseKey_(p,key,file,lcnt,i,e),
			ll_tParser2::parseKey_(p,key,file,lcnt,i,e) {}
	};

public:
	// print help struct
	struct printHelp_: public virtual ll_tParser1::printHelp_,
       			   public virtual ll_tParser2::printHelp_ {
	// constructor
	inline explicit printHelp_(const ll_tParser3& p, std::ostream& os):
			aux_parser::printHelp_(p,os),
			ll_tParser1::printHelp_(p,os),
			ll_tParser2::printHelp_(p,os) {}
	};
};

#endif // _LL_TPARSER3_
