// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_omen_prepper_all.h"
#include "ll_omen.h"
#include <sstream>

using namespace lm__;
using namespace ll__;
using namespace aux;


void test_omen_prepper_all::test_ctor_exceptions() {
	const std::string wpath = "data/w90/mos2/";
	const std::string woutp = "outp/";

	ll_wf_input inp; inp.verbosity = 3;
	std::stringstream sstr;

	// bad wout
	CPPUNIT_ASSERT_THROW(omen::prepper(inp,sstr),std::invalid_argument);
	inp.wout = wpath+"wannier90.wout";

	// bad hrdat
	CPPUNIT_ASSERT_THROW(omen::prepper(inp,sstr),std::invalid_argument);
	inp.hrdat = wpath+"wannier90_hr.dat";

	// bad mode
	inp.mode = {"nothing"};
	CPPUNIT_ASSERT_THROW(omen::prepper(inp,sstr),std::invalid_argument);
	inp.mode = {"atomic"};


	// bad C, bad size
	inp.C = fMat({1.0, 0.0,
		      0.0,-1.0,
		      0.0, 0.0},3);
	CPPUNIT_ASSERT_THROW(omen::prepper(inp,sstr),std::invalid_argument);
	
	// bad C, not integers
	inp.C = fMat({1.1, 0.2, 1.2,
		      0.3,-1.1, 2.3,
		      0.9, 0.4, 1.3},3);
	CPPUNIT_ASSERT_THROW(omen::prepper(inp,sstr),std::invalid_argument);
	inp.C = fMat();
	

	// bad S, bad size
	inp.S = fMat({1.0, 0.0,
		      0.0,-1.0,
		      0.0, 0.0},3);
	CPPUNIT_ASSERT_THROW(omen::prepper(inp,sstr),std::invalid_argument);
	
	// bad S, bad det zero
	inp.S = fMat({1.0, 0.0, 0.0,
		      0.0,-1.0, 0.0,
		      0.0, 0.0, 0.0},3);
	CPPUNIT_ASSERT_THROW(omen::prepper(inp,sstr),std::invalid_argument);

	// bad S, bad det negative
	inp.S = fMat({1.0, 0.0, 0.0,
		      0.0,-1.0, 0.0,
		      0.0, 0.0, 1.0},3);
	CPPUNIT_ASSERT_THROW(omen::prepper(inp,sstr),std::invalid_argument);
	inp.S = fMat();

	// bad C, resulting basis not orthorhombic
	inp.C = fMat({1.0, 0.0, 0.0,
		      0.0, 1.0, 0.0,
		      0.0, 0.0, 1.0},3,3);
	CPPUNIT_ASSERT_THROW(omen::prepper(inp,sstr),std::invalid_argument);
	inp.C = fMat();

	// bad xyz, bad size
	inp.xyz = fMat({1.0, 0.0,
		        0.0,-1.0,
		        0.0, 0.0},3);
	CPPUNIT_ASSERT_THROW(omen::prepper(inp,sstr),std::invalid_argument);
	
	// bad xyz, not orthogonal
	inp.xyz = fMat({1.0, 0.0, 1.0,
		        0.0,-1.0, 0.0,
		        0.0, 0.0, 1.0},3);
	CPPUNIT_ASSERT_THROW(omen::prepper(inp,sstr),std::invalid_argument);
	inp.xyz = eye<fMat>(3,3);
}

void test_omen_prepper_all::test_ctor_properties() {

	const std::string wpath = "data/w90/mos2/";
	const std::string woutp = "outp/";

	ll_wf_input inp; inp.verbosity = 3;
	std::stringstream sstr;

	inp.wout = wpath+"wannier90.wout";
	inp.hrdat = wpath+"wannier90_hr.dat";
	inp.weig = wpath+"wannier90.eig";
	inp.outcar = wpath+"OUTCAR";
	inp.prefix = woutp;

	omen::prepper P(inp,sstr);

	// cell is orthorhombic, positively aligned with cartesian axis
	CPPUNIT_ASSERT(!P.ORcell.empty());
	CPPUNIT_ASSERT_EQUAL(size_t(DIM__),P.ORcell.dim());
	CPPUNIT_ASSERT(P.ORcell.B().ob());
	CPPUNIT_ASSERT_EQUAL(double(1.0),P.ORcell.sign());
	CPPUNIT_ASSERT(diag(diag(P.ORcell.B())) == P.ORcell.B());

	// check types in wbh are a subset of those in cell
	{
		const auto t1 = P.ORcell.id();
		const auto t2 = P.W.cell().id();
		CPPUNIT_ASSERT(std::includes(t1.cbegin(),t1.cend(),
					     t2.cbegin(),t2.cend()));
	}

	// check r, I
	CPPUNIT_ASSERT_EQUAL(P.ORcell.dim(),P.r.size());

	// check C and NNE
	CPPUNIT_ASSERT(P.C == round(P.C));
	CPPUNIT_ASSERT(P.W.cell().B().prod(P.C).ob());
	CPPUNIT_ASSERT(P.NNE == round(P.NNE));
	CPPUNIT_ASSERT(P.NNE == diag(diag(P.NNE)));
}


const char* test_omen_prepper_all::test_id() noexcept {
	return "test_omen_prepper_all";
}

CppUnit::Test* test_omen_prepper_all::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_omen_prepper_all>(
		"test_ctor_exceptions", &test_omen_prepper_all::test_ctor_exceptions));
	suite->addTest(new CppUnit::TestCaller<test_omen_prepper_all>(
		"test_ctor_properties", &test_omen_prepper_all::test_ctor_properties));
	
	return suite;
}
