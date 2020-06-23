// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_BStest_all.h"
#include "ll_BStest.h"
#include "ll_hbonds.h"
#include "ll_io.h"
#include "ll_fn.h"
#include "ll_testTools.h"
#include "lm_testTools.h"
#include <sstream>


using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;


void test_BStest_all::test_exceptions() {
	const std::string path = "data/w90/mos2/";
	std::stringstream sstr;
	
	ll_wmatching_input winp;
	winp.wout = path+WOUT__;
	winp.hrdat = path+HR__;
	winp.mode = {"atomic"};
	winp.tol = 0.0;
	winp.rm_unmatched = 1;
	winp.verbosity = 0;
	winp.prefix = "outp/";

	const ll_hbonds W(winp,std::cout);
	
	fMat EXP = eye<fMat>(DIM__,DIM__);
	rv r(DIM__,false);
	vb_cb E = {0.0,0.0};

	ll_BStest_input inp;
	inp.verbosity = 3;

	// mesh
	{
		// negative rho_k
		inp.rho_k = -1.0;
		CPPUNIT_ASSERT_THROW(meshBStest(W,EXP,r,E,inp,sstr),
			std::invalid_argument);
		
		// very high rho_k
		inp.rho_k = 1e8;
		CPPUNIT_ASSERT_THROW(meshBStest(W,EXP,r,E,inp,sstr),
			std::invalid_argument);
		
		inp.rho_k = 1000.0;

		// bad expansion, empty
		EXP = fMat();
		CPPUNIT_ASSERT_THROW(meshBStest(W,EXP,r,E,inp,sstr),
			std::invalid_argument);
		
		// bad expansion, det 0
		EXP = ones<fMat>(DIM__,DIM__); EXP[0] = 0.0;
		CPPUNIT_ASSERT_THROW(meshBStest(W,EXP,r,E,inp,sstr),
			std::invalid_argument);
		
		// bad expansion, not integers
		EXP = ones<fMat>(DIM__,DIM__); EXP[0] = 1.4;
		CPPUNIT_ASSERT_THROW(meshBStest(W,EXP,r,E,inp,sstr),
			std::invalid_argument);
	}
	
	// trace
	{
		// very high Nk
		inp.Nk = 1e8;
		CPPUNIT_ASSERT_THROW(traceBStest(W,EXP,r,inp,sstr),
			std::invalid_argument);
		
		inp.Nk = 1000.0;

		// bad expansion, empty
		EXP = fMat();
		CPPUNIT_ASSERT_THROW(traceBStest(W,EXP,r,inp,sstr),
			std::invalid_argument);
		
		// bad expansion, det 0
		EXP = ones<fMat>(DIM__,DIM__); EXP[0] = 0.0;
		CPPUNIT_ASSERT_THROW(traceBStest(W,EXP,r,inp,sstr),
			std::invalid_argument);
		
		// bad expansion, not integers
		EXP = ones<fMat>(DIM__,DIM__); EXP[0] = 1.4;
		CPPUNIT_ASSERT_THROW(traceBStest(W,EXP,r,inp,sstr),
			std::invalid_argument);

		EXP = eye<fMat>(DIM__,DIM__);

		// bad kpts, empty
		inp.kpts = fMat();
		CPPUNIT_ASSERT_THROW(traceBStest(W,EXP,r,inp,sstr),
			std::invalid_argument);
		
		// bad kpts, empty
		inp.kpts = rand<fMat>(DIM__-1,genRndST());
		CPPUNIT_ASSERT_THROW(traceBStest(W,EXP,r,inp,sstr),
			std::invalid_argument);
	}

}
void test_BStest_all::test_random_expansions() {
	const std::string path = "data/w90/mos2/";

	ll_wmatching_input inp;
	inp.wout = path+WOUT__;
	inp.hrdat = path+HR__;
	inp.mode = {"atomic"};
	inp.tol = 0.0;
	inp.rm_unmatched = 1;
	inp.verbosity = 0;
	inp.prefix = "outp/";

	const ll_hbonds wbh(inp,std::cout);
	
	const rv r = {false,true,false};
	const auto E = findBandEdges(readEf(path+OUTCAR__),readEig(path+WEIG__));

	// mesh
	{
		std::stringstream sstr;

		ll_BStest_input inp;
		inp.bzbounds = {.0,.5};
		inp.rho_k = genRndDouble(1000.0,3000.0);
		inp.hrdat = path+HR__;
		inp.re = true;
		inp.Lvb = 2.0;
		inp.Lcb = 2.0;
		inp.toldev = 30.0;
		inp.Nthreads = 4;
		inp.verbosity = 3;
		inp.prefix = "outp/";

		fMat EXP;
		do EXP = randi<fMat>(2,2,-3,3);
		while (det(EXP)==0.0);

		EXP = fMat({EXP(0,0),0.0,EXP(1,0),
			        0.0, 1.0,    0.0,
			    EXP(0,1),0.0,EXP(1,1)},3,3);

		meshBStest(wbh,EXP,r,E,inp,sstr);

		const std::string outp = sstr.str();
		const size_t pos1 = outp.find("max deviations in window",0);
		const size_t pos2 = outp.find("]eV:",pos1)+12;
		const size_t pos3 = outp.find("meV",pos2);
		const double err_vb = std::stod(outp.substr(pos2,pos3-pos2));

		const size_t pos4 = outp.find("]eV:",pos3)+12;
		const size_t pos5 = outp.find("meV",pos2);
		const double err_cb = std::stod(outp.substr(pos4,pos5-pos4));

		CPPUNIT_ASSERT(std::abs(err_vb)<3.0);
		CPPUNIT_ASSERT(std::abs(err_cb)<3.0);
	}
}


const char* test_BStest_all::test_id() noexcept {
	return "test_BStest";
}

CppUnit::Test* test_BStest_all::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_BStest_all>(
		"test_exceptions", &test_BStest_all::test_exceptions));
	suite->addTest(new CppUnit::TestCaller<test_BStest_all>(
		"test_random_expansions", &test_BStest_all::test_random_expansions));
	
	return suite;
}
