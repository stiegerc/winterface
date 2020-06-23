// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_fn_wannier_matching.h"
#include "ll_fn.h"
#include "ll_io.h"
#include "ll_cell.h"
#include "aux_sort.h"
#include "ll_testTools.h"
#include "lm_testTools.h"

using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;


void test_fn_wannier_matching::test_matchCenters() {
	const std::string ipath = "data/w90/mos2/";
	const std::string prefix = "outp/mos2_";

	const auto B = readB(ipath+"wannier90.wout");
	const auto Ap_ = readAp(ipath+"wannier90.wout");
	const auto Wp_ = readWp(ipath+"wannier90.wout");

	auto Ap = Ap_.Ap();
	const auto Wp = B.leftDivide(Wp_.Wp());
	const auto I = matchCenters(B,Ap,Wp,NN);
	CPPUNIT_ASSERT_EQUAL(Ap_.N(),Ap.N());

	printPOSCAR(prefix+"Ap.psc",B,Ap_.Ap(),Ap_.id(),1.0,true,16,"atomic positions");
	printPOSCAR(prefix+"Wp.psc",B,Wp,aTv(Wp.N(),0),{},1.0,true,16,"wannier positions");	
	printPOSCAR(prefix+"Wp_matched.psc",B,Wp,wiToT(I),
			{"H","P","N","S","O","C"},1.0,true,16,"wannier positions matched");

	const auto D = genDmat(B,Ap,Wp%1.0,NN);
	const auto M = mmin(D);
	for(size_t i=0; i!=M.pos.size(); ++i)
		CPPUNIT_ASSERT(std::find(I[M.pos[i]].cbegin(),I[M.pos[i]].cend(),i)!=I[M.pos[i]].cend());

	CPPUNIT_ASSERT_EQUAL(size_t(5),I[0].size());	// Mo
	CPPUNIT_ASSERT_EQUAL(size_t(5),I[1].size());	// Mo
	CPPUNIT_ASSERT_EQUAL(size_t(3),I[2].size());	// S
	CPPUNIT_ASSERT_EQUAL(size_t(3),I[3].size());	// S
	CPPUNIT_ASSERT_EQUAL(size_t(3),I[4].size());	// S
	CPPUNIT_ASSERT_EQUAL(size_t(3),I[5].size());	// S
}

void test_fn_wannier_matching::test_clusterize() {
	const std::string ipath = "data/w90/mos2/";
	const std::string prefix = "outp/mos2_";

	const double eps = .5;
	const size_t minpts = 2;

	const auto B = readB(ipath+"wannier90.wout");
	const auto Ap_ = readAp(ipath+"wannier90.wout");
	const auto Wp_ = readWp(ipath+"wannier90.wout");

	auto Ap = Ap_.Ap();
	const auto Wp = B.leftDivide(Wp_.Wp());
	const auto I = matchCenters(B,Ap,Wp,NN);
	const auto J = clusterize(B,Wp,eps,minpts,NN);

	printPOSCAR(prefix+"Wp_cluster.psc",B,Wp,wiToT(J),
			{"H","P","N","S","O","C"},1.0,true,16,"wannier positions clustered");

	for (auto i: I)
		CPPUNIT_ASSERT(std::find(J.begin(),J.end(),i)!=I.end());
}

void test_fn_wannier_matching::test_genCenters() {
	const std::string ipath = "data/w90/mos2/";
	const std::string prefix = "outp/mos2_";

	const double eps = .5;
	const size_t minpts = 2;

	const auto B = readB(ipath+"wannier90.wout");
	const auto Ap_ = readAp(ipath+"wannier90.wout");
	const auto Wp_ = readWp(ipath+"wannier90.wout");

	auto Ap = Ap_.Ap();
	const auto Wp = B.leftDivide(Wp_.Wp());
	const auto J = clusterize(B,Wp,eps,minpts,NN);
	
	const auto Cp =  genCenters(B,Wp,Wp_.s(),J);
	printPOSCAR(prefix+"Cp.psc",B,Cp,aTv(J.size(),1),
			{"As","As","As","As","As","As"},1.0,true,16,"cluster centers");

	const auto D = genDmat(B,Ap,Cp,NN);
	const auto m_D = mmin(D);
	CPPUNIT_ASSERT(all(m_D.mat.lt(.1)));
}


const char* test_fn_wannier_matching::test_id() noexcept {
	return "test_fn_wannier_matching";
}

CppUnit::Test* test_fn_wannier_matching::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_fn_wannier_matching>(
		"test_matchCenters", &test_fn_wannier_matching::test_matchCenters));
	suite->addTest(new CppUnit::TestCaller<test_fn_wannier_matching>(
		"test_clusterize", &test_fn_wannier_matching::test_clusterize));
	suite->addTest(new CppUnit::TestCaller<test_fn_wannier_matching>(
		"test_genCenters", &test_fn_wannier_matching::test_genCenters));

	return suite;
}
