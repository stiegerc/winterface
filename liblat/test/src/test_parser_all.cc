// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_parser_all.h"
#include "ll_tParser1.h"
#include "ll_tParser3.h"
#include "ll_testTools.h"
#include "lm_testTools.h"
#include <fstream>
#include <iomanip>
#include "ll_io.h"
#include "aux_io.h"

using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;
using namespace aux;


void test_parser_all::test_tParser() {
	const std::string path = "data/scripts/";
	
	// check parseFile exceptions
	{
		// check file not found exception
		CPPUNIT_ASSERT_THROW(parseFile<ll_tParser1>("blub"),std::invalid_argument);

		// check reading bad key script exception
		CPPUNIT_ASSERT_THROW(parseFile<ll_tParser1>(path+"bad_key"),std::invalid_argument);
	

		// strings
		// check read wout missing argument exception
		CPPUNIT_ASSERT_THROW(parseFile<ll_tParser1>(path+"bad_tString"),std::invalid_argument);

	
		// ints
		// check read uInts missing argument exception
		CPPUNIT_ASSERT_THROW(parseFile<ll_tParser1>(path+"bad_uInts_arg"),std::invalid_argument);
	
		// check read uInts bad ints exception
		CPPUNIT_ASSERT_THROW(parseFile<ll_tParser1>(path+"bad_uInts_ints"),std::invalid_argument);
	
		// check read uInts bad <0 ints exception
		CPPUNIT_ASSERT_THROW(parseFile<ll_tParser1>(path+"bad_uInts_nints"),std::invalid_argument);


		// bools
		// check read tBools missing argument exception
		CPPUNIT_ASSERT_THROW(parseFile<ll_tParser1>(path+"bad_tBools_arg"),std::invalid_argument);

		// check read tBools bad bools exception
		CPPUNIT_ASSERT_THROW(parseFile<ll_tParser1>(path+"bad_tBools_bools"),std::invalid_argument);
	

		// mats
		// check read tMat from bad file
		CPPUNIT_ASSERT_THROW(parseFile<ll_tParser1>(path+"bad_tMat_file"),std::invalid_argument);

		// check read tMat variable entries from script
		CPPUNIT_ASSERT_THROW(parseFile<ll_tParser1>(path+"bad_tMat_var"),std::invalid_argument);
	
		// check read tMat bad num from script
		CPPUNIT_ASSERT_THROW(parseFile<ll_tParser1>(path+"bad_tMat_num"),std::invalid_argument);
	}
	{
		CPPUNIT_ASSERT_NO_THROW(parseFile<ll_tParser1>(path+"comments_only"));
		CPPUNIT_ASSERT_NO_THROW(parseFile<ll_tParser1>(path+"read_mat"));
	}
	
	
	// helper lambdas
	const auto rndString = [](const size_t l, const size_t u) -> std::string {
		std::string res;
		
		const size_t N = genRndST(l,u);
		const auto rnd = rand<fMat>(1,N,97.0,122.99);
		res.reserve(N);
		for (const auto i: rnd)
			res.push_back(char(i));
		return res;
	};
	const auto rndStrings = [&rndString](const size_t N,
			const size_t l, const size_t u) -> std::vector<std::string> {
		std::vector<std::string> res; res.reserve(genRndST(1,N));
		while (res.size()<res.capacity())
			res.push_back(rndString(l,u));
		return res;
	};
	const auto rndDouble = [](const double l, const double u) -> double { return genRndDouble(l,u); };
	const auto rndUint = [](const size_t l, const size_t u) -> size_t { return genRndST(l,u); };
	const auto rndBool = []() -> bool { return genRndST(0,1); };
	const auto rndDoubles = [](const size_t sl, const size_t su, const double l, const double u)
			-> std::vector<double> {
		std::vector<double> res(genRndST(sl,su));
		for (auto& i: res) i = genRndDouble(l,u);
		return res;
	};
	const auto rndUints = [](const size_t sl, const size_t su, const size_t l, const size_t u)
			-> std::vector<size_t> {
		std::vector<size_t> res(genRndST(sl,su));
		for (auto& i: res) i = genRndST(l,u);
		return res;
	};
	const auto rndBools = [](const size_t sl, const size_t su) -> std::vector<bool> {
		std::vector<bool> res(genRndST(sl,su));
		for (size_t i=0; i!=res.size(); ++i)
			res[i] = genRndST(0,1);
		return res;
	};
	const auto rndMat = [](const size_t lm, const size_t um,
			 const size_t ln, const size_t un,
			 const double l, const double u) -> fMat {
		return rand<fMat>(genRndST(lm,um),genRndST(ln,un),l,u);
	};
	const auto treatUints = [](std::vector<size_t> inp) -> std::vector<size_t> {
		std::sort(inp.begin(),inp.end());
		inp.erase(std::unique(inp.begin(),inp.end()),inp.end());
		return inp;
	};


	// generate random input
	const std::string tString1 = rndString(10,20);
	const std::vector<std::string> tStrings1 = rndStrings(genRndST(1,5),3,8);
	const double tDouble1 = rndDouble(-10.0,10.0);
	const size_t tUint1 = rndUint(0,999);
	const bool tBool1 = rndBool();
	const std::vector<double> tDoublesF1 = rndDoubles(DIM__,DIM__,-1.0,1.0);
	const std::vector<double> tDoublesV1 = rndDoubles(1,10,-1.0,1.0);
	const std::vector<size_t> tUintsFNT1 = rndUints(DIM__,DIM__,0,DIM__);
	const std::vector<size_t> tUintsFT1 = rndUints(DIM__,DIM__,0,DIM__);
	const std::vector<size_t> tUintsVNT1 = rndUints(1,10,0,10);
	const std::vector<size_t> tUintsVT1 = rndUints(1,10,0,10);
	const std::vector<bool> tBoolsF1 = rndBools(DIM__,DIM__);
	const std::vector<bool> tBoolsV1 = rndBools(1,10);
	const auto tMat1 = rndMat(1,8,1,6,-1.0,1.0);
	const auto tMatF1 = rndMat(1,8,1,6,-1.0,1.0);

	const size_t lcnt = 20+tMat1.N();
	const size_t kcnt = 16;

	// write random script file
	{
		std::ofstream file;
		file.open(path+"rnd1");

		file << "#randomly generated script\n\n";
		file << "tString1 = " << tString1 << "\n";
		file << "tStrings1 = " << tStrings1 << "\n";
		file << "tDouble1 = " << tDouble1 << "\n";
		file << "tUint1 = " << tUint1 << "\n";
		file << "tBool1 = " << tBool1 << "\n";
		file << "tDoublesF1 = " << tDoublesF1 << "\n";
		file << "tDoublesV1 = " << tDoublesV1 << "\n";
		file << "tUintsFNT1 = " << tUintsFNT1 << "\n";
		file << "tUintsFT1 = " << tUintsFT1 << "\n";
		file << "tUintsVNT1 = " << tUintsVNT1 << "\n";
		file << "tUintsVT1 = " << tUintsVT1 << "\n";
		file << "tBoolsF1 = " << tBoolsF1 << "\n";
		file << "tBoolsV1 = " << tBoolsV1 << "\n";
		file << "tMat1\n" << lm__::T(tMat1) << "\n\n";
		file << "tMatF1 = " << path << "tMatF1" << "\n";
		file << "\nverbosity = 0";

		file.close();
		lm__::T(tMatF1).printToFile(path+"tMatF1");
	}


	// read in random file and check
	{
		// set mtol_
		set_mtol(1e-5);

		CPPUNIT_ASSERT_NO_THROW(parseFile<ll_tParser1>(path+"rnd1"));
		
		size_t lcnt_, kcnt_;
		const auto P = parseFile<ll_tParser1>(path+"rnd1",lcnt_,kcnt_);

		CPPUNIT_ASSERT_EQUAL(lcnt,lcnt_);
		CPPUNIT_ASSERT_EQUAL(kcnt,kcnt_);
		
		CPPUNIT_ASSERT_EQUAL(tString1,P.tString1);
		CPPUNIT_ASSERT(tStrings1==P.tStrings1);
		CPPUNIT_ASSERT_DELTA(tDouble1,P.tDouble1,mtol());
		CPPUNIT_ASSERT_EQUAL(tUint1,P.tUint1);
		CPPUNIT_ASSERT_EQUAL(tBool1,P.tBool1);

		CPPUNIT_ASSERT_EQUAL(tDoublesF1.size(),P.tDoublesF1.size());
		for (size_t i=0; i!=tDoublesF1.size(); ++i)
			CPPUNIT_ASSERT_DELTA(tDoublesF1[i],P.tDoublesF1[i],mtol());
		
		CPPUNIT_ASSERT_EQUAL(tDoublesV1.size(),P.tDoublesV1.size());
		for (size_t i=0; i!=tDoublesV1.size(); ++i)
			CPPUNIT_ASSERT_DELTA(tDoublesV1[i],P.tDoublesV1[i],mtol());
		
		CPPUNIT_ASSERT(tUintsFNT1==P.tUintsFNT1);
		CPPUNIT_ASSERT(treatUints(tUintsFT1)==P.tUintsFT1);
		CPPUNIT_ASSERT(tUintsVNT1==P.tUintsVNT1);
		CPPUNIT_ASSERT(treatUints(tUintsVT1)==P.tUintsVT1);
	
		CPPUNIT_ASSERT(tBoolsF1==P.tBoolsF1);
		CPPUNIT_ASSERT(tBoolsV1==P.tBoolsV1);
		
		CPPUNIT_ASSERT(tMat1==P.tMat1);
		CPPUNIT_ASSERT(tMatF1==P.tMatF1);

		// reset mtol_
		reset_mtol();
	}



	// generate 2nd set of random input
	const std::string tString2 = rndString(10,20);
	const std::vector<std::string> tStrings2 = rndStrings(genRndST(1,5),3,8);
	const double tDouble2 = rndDouble(-10.0,10.0);
	const size_t tUint2 = rndUint(0,999);
	const bool tBool2 = rndBool();
	const std::vector<double> tDoublesF2 = rndDoubles(DIM__,DIM__,-1.0,1.0);
	const std::vector<double> tDoublesV2 = rndDoubles(1,10,-1.0,1.0);
	const std::vector<size_t> tUintsFNT2 = rndUints(DIM__,DIM__,0,DIM__);
	const std::vector<size_t> tUintsFT2 = rndUints(DIM__,DIM__,0,DIM__);
	const std::vector<size_t> tUintsVNT2 = rndUints(1,10,0,10);
	const std::vector<size_t> tUintsVT2 = rndUints(1,10,0,10);
	const std::vector<bool> tBoolsF2 = rndBools(DIM__,DIM__);
	const std::vector<bool> tBoolsV2 = rndBools(1,10);
	const auto tMat2 = rndMat(1,8,1,6,-1.0,1.0);
	const auto tMatF2 = rndMat(1,8,1,6,-1.0,1.0);

	// write random script file using both sets
	{
		std::ofstream file;
		file.open(path+"rnd3");

		file << "#randomly generated script\n\n";
		file << "tString1 = " << tString1 << "\n";
		file << "tStrings1 = " << tStrings1 << "\n";
		file << "tDouble1 = " << tDouble1 << "\n";
		file << "tUint1 = " << tUint1 << "\n";
		file << "tBool1 = " << tBool1 << "\n";
		file << "tDoublesF1 = " << tDoublesF1 << "\n";
		file << "tDoublesV1 = " << tDoublesV1 << "\n";
		file << "tUintsFNT1 = " << tUintsFNT1 << "\n";
		file << "tUintsFT1 = " << tUintsFT1 << "\n";
		file << "tUintsVNT1 = " << tUintsVNT1 << "\n";
		file << "tUintsVT1 = " << tUintsVT1 << "\n";
		file << "tBoolsF1 = " << tBoolsF1 << "\n";
		file << "tBoolsV1 = " << tBoolsV1 << "\n";
		file << "tMat1\n" << lm__::T(tMat1) << "\n\n";
		file << "tMatF1 = " << path << "tMatF1" << "\n";
		file << "\nverbosity = 0\n";
		file << "tString2 = " << tString2 << "\n";
		file << "tStrings2 = " << tStrings2 << "\n";
		file << "tDouble2 = " << tDouble2 << "\n";
		file << "tUint2 = " << tUint2 << "\n";
		file << "tBool2 = " << tBool2 << "\n";
		file << "tDoublesF2 = " << tDoublesF2 << "\n";
		file << "tDoublesV2 = " << tDoublesV2 << "\n";
		file << "tUintsFNT2 = " << tUintsFNT2 << "\n";
		file << "tUintsFT2 = " << tUintsFT2 << "\n";
		file << "tUintsVNT2 = " << tUintsVNT2 << "\n";
		file << "tUintsVT2 = " << tUintsVT2 << "\n";
		file << "tBoolsF2 = " << tBoolsF2 << "\n";
		file << "tBoolsV2 = " << tBoolsV2 << "\n";
		file << "tMat2\n" << lm__::T(tMat2) << "\n\n";
		file << "tMatF2 = " << path << "tMatF2" << "\n";

		file.close();
		lm__::T(tMatF2).printToFile(path+"tMatF2");
	}


	// read in random file and check
	{
		// set mtol_
		set_mtol(1e-5);

		CPPUNIT_ASSERT_NO_THROW(parseFile<ll_tParser3>(path+"rnd3"));
		
		const auto P = parseFile<ll_tParser3>(path+"rnd3");
		
		CPPUNIT_ASSERT_EQUAL(tString1,P.tString1);
		CPPUNIT_ASSERT(tStrings1==P.tStrings1);
		CPPUNIT_ASSERT_DELTA(tDouble1,P.tDouble1,mtol());
		CPPUNIT_ASSERT_EQUAL(tUint1,P.tUint1);
		CPPUNIT_ASSERT_EQUAL(tBool1,P.tBool1);

		CPPUNIT_ASSERT_EQUAL(tDoublesF1.size(),P.tDoublesF1.size());
		for (size_t i=0; i!=tDoublesF1.size(); ++i)
			CPPUNIT_ASSERT_DELTA(tDoublesF1[i],P.tDoublesF1[i],mtol());
		
		CPPUNIT_ASSERT_EQUAL(tDoublesV1.size(),P.tDoublesV1.size());
		for (size_t i=0; i!=tDoublesV1.size(); ++i)
			CPPUNIT_ASSERT_DELTA(tDoublesV1[i],P.tDoublesV1[i],mtol());
		
		CPPUNIT_ASSERT(tUintsFNT1==P.tUintsFNT1);
		CPPUNIT_ASSERT(treatUints(tUintsFT1)==P.tUintsFT1);
		CPPUNIT_ASSERT(tUintsVNT1==P.tUintsVNT1);
		CPPUNIT_ASSERT(treatUints(tUintsVT1)==P.tUintsVT1);
	
		CPPUNIT_ASSERT(tBoolsF1==P.tBoolsF1);
		CPPUNIT_ASSERT(tBoolsV1==P.tBoolsV1);
		
		CPPUNIT_ASSERT(tMat1==P.tMat1);
		CPPUNIT_ASSERT(tMatF1==P.tMatF1);

		CPPUNIT_ASSERT_EQUAL(tString2,P.tString2);
		CPPUNIT_ASSERT(tStrings2==P.tStrings2);
		CPPUNIT_ASSERT_DELTA(tDouble2,P.tDouble2,mtol());
		CPPUNIT_ASSERT_EQUAL(tUint2,P.tUint2);
		CPPUNIT_ASSERT_EQUAL(tBool2,P.tBool2);

		CPPUNIT_ASSERT_EQUAL(tDoublesF2.size(),P.tDoublesF2.size());
		for (size_t i=0; i!=tDoublesF2.size(); ++i)
			CPPUNIT_ASSERT_DELTA(tDoublesF2[i],P.tDoublesF2[i],mtol());
		
		CPPUNIT_ASSERT_EQUAL(tDoublesV2.size(),P.tDoublesV2.size());
		for (size_t i=0; i!=tDoublesV2.size(); ++i)
			CPPUNIT_ASSERT_DELTA(tDoublesV2[i],P.tDoublesV2[i],mtol());
		
		CPPUNIT_ASSERT(tUintsFNT2==P.tUintsFNT2);
		CPPUNIT_ASSERT(treatUints(tUintsFT2)==P.tUintsFT2);
		CPPUNIT_ASSERT(tUintsVNT2==P.tUintsVNT2);
		CPPUNIT_ASSERT(treatUints(tUintsVT2)==P.tUintsVT2);
	
		CPPUNIT_ASSERT(tBoolsF2==P.tBoolsF2);
		CPPUNIT_ASSERT(tBoolsV2==P.tBoolsV2);
		
		CPPUNIT_ASSERT(tMat2==P.tMat2);
		CPPUNIT_ASSERT(tMatF2==P.tMatF2);


		// reset mtol_
		reset_mtol();
	}
}

void test_parser_all::test_screenFile() {

	const std::string path = "data/scripts/";

	// dump file into stringstream lambda
	const auto dfsstr = [](const std::string& fileName) -> std::stringstream {

		auto file = aux::openFile<std::ifstream>(fileName);

		std::string line;
		std::stringstream res;
		while (std::getline(file,line))
			res << line << "\n";
		file.close();
		return res;
	};

	auto ref = dfsstr(path+"rnd1");

	// try without any keys
	{
		std::stringstream sstr;
		aux::screenFile(sstr,path+"rnd1",{});
		CPPUNIT_ASSERT(ref.str()==sstr.str());
	}

	// try with a random subset of inserted keys
	{
		// generate random subset of keys in file rnd
		std::vector<std::string> keys = {
			"tString1","tDouble1","tUint1","tBool1",
			"tDoublesF1","tDoublesV1","tUintsFNT1","tUintsFT1","tUintsVNT1",
			"tUintsVT1","tBoolsF1","tBoolsV1","tMat1","tMatF1","verbosity"	
		};
		
		std::random_device rd;
    		std::mt19937 g(rd());
		std::shuffle(keys.begin(),keys.end(),g);
		
		const size_t n = genRndST(0,keys.size());
		const std::vector<std::string> keys_rm(keys.begin(),keys.begin()+n);
		const std::vector<std::string> keys_no_rm(keys.begin()+n,keys.end());

		const std::regex rgx(RGX__);
		
		// check no keys in keys_rm and all keys in keys_no_rm occur in the screened file
		{
			std::stringstream sstr;
			aux::screenFile(sstr,path+"rnd1",keys_rm);
			const auto filestr = sstr.str();

			// tokenize filestr
			std::sregex_token_iterator i(filestr.begin(),filestr.end(),rgx,-1), e;
			
			CPPUNIT_ASSERT(std::none_of(keys_rm.cbegin(),keys_rm.cend(),
				[i,e](const std::string& key) -> bool {
					return std::find(i,e,key)!=e;
				})
			);
			CPPUNIT_ASSERT(std::all_of(keys_no_rm.cbegin(),keys_no_rm.cend(),
				[i,e](const std::string& key) -> bool {
					return std::find(i,e,key)!=e;
				})
			);
		}

		// write screened file to 'rnd_screened', add removed keys to bottom
		{
			const std::string fileName = path+"rnd_screened";

			auto file = aux::openFile<std::ofstream>(fileName);
			
			aux::screenFile(file,path+"rnd1",keys_rm);
			
			file << "\n\nKEYS REMOVED:\n" << keys_rm;
			file.close();
		}
	}
}


const char* test_parser_all::test_id() noexcept {
	return "test_parser_all";
}

CppUnit::Test* test_parser_all::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_parser_all>(
		"test_tParser", &test_parser_all::test_tParser));
	suite->addTest(new CppUnit::TestCaller<test_parser_all>(
		"test_screenFile", &test_parser_all::test_screenFile));

	return suite;
}
