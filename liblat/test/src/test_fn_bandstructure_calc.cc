// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_fn_bandstructure_calc.h"
#include "ll_fn.h"
#include "ll_io.h"
#include "ll_cell.h"
#include "ll_testTools.h"
#include "lm_testTools.h"


using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;


void test_fn_bandstructure_calc::test_findBandEdges() {
	const std::string path = "data/w90/";

	// mos2_dl
	{
		const double Ef = readEf(path+"mos2/OUTCAR");
		const auto E = readEig(path+"mos2/wannier90.eig");
		const auto bs = findBandEdges(Ef,E);

		const double vb = max(E.getl(E.lt(Ef)));
		const double cb = min(E.getl(E.gt(Ef)));
	
		CPPUNIT_ASSERT_EQUAL(vb,bs.vb);
		CPPUNIT_ASSERT_EQUAL(cb,bs.cb);
	}
	
	// inas_mono
	{
		const double Ef = readEf(path+"inas_mono/OUTCAR");
		const auto E = readEig(path+"inas_mono/wannier90.eig");
		const auto bs = findBandEdges(Ef,E);

		const double vb = max(E.getl(E.lt(Ef)));
		const double cb = min(E.getl(E.gt(Ef)));
	
		CPPUNIT_ASSERT_EQUAL(vb,bs.vb);
		CPPUNIT_ASSERT_EQUAL(cb,bs.cb);
	}
}
void test_fn_bandstructure_calc::test_calcBS() {
	const size_t Np = 301;

	// mos2_dl
	{
		const std::string path = "data/w90/chk/mos2_dl/";
		const std::string prefix = "outp/mos2_";
		
		const auto H = readHr(path+"wannier90_hr.dat");
		
		const fMat P({-.5,.0,.0,.5,.0,.0},3);
		const auto k = genPath(P,Np);

		const auto E = calcBS(H,k.path());

		E.printToFile(prefix+"E.dat");
		k.path().printToFile(prefix+"kpath.dat");
		k.pos().printToFile(prefix+"kpos.dat");
	}
	
	// inas_mono
	{
		const std::string path = "data/w90/chk/inas_mono/";
		const std::string prefix = "outp/inas_mono_";
		
		const auto H = readHr(path+"wannier90_hr.dat");
		
		const fMat P({-.5,.0,.0,.5,.0,.0},3);
		const auto k = genPath(P,Np);

		const auto E = calcBS(H,k.path());

		E.printToFile(prefix+"E.dat");
		k.path().printToFile(prefix+"kpath.dat");
		k.pos().printToFile(prefix+"kpos.dat");
	}
}
void test_fn_bandstructure_calc::test_calcFoldedBS() {
	const size_t Np = 301;
	const std::string path = "data/w90/mos2/";
	const std::string prefix = "outp/mos2_folded_";

	const auto H = readHr(path+"wannier90_hr.dat");
	const auto B = readB(path+"wannier90.wout");
	const fMat P({-.5,.0,.0,.5,.0,.0},3);
	const auto k = genPath(P,Np);
	const size_t Gi = std::find(k.ccBegin(),
			k.ccEnd(),fMat({.0,.0,.0})).i();

	// primitive cell using calcBS
	const auto Eprim = calcBS(H,k.path());
	Eprim.printToFile(prefix+"E_prim.dat");
	k.path().printToFile(prefix+"kpath_prim.dat");
	k.pos().printToFile(prefix+"kpos_prim.dat");
	
	// primitive cell using calcFoldedBS
	const auto Eprimf = calcFoldedBS(H,k.path(),B,B);
	Eprimf.printToFile(prefix+"E_primf.dat");
	k.path().printToFile(prefix+"kpath_primf.dat");
	k.pos().printToFile(prefix+"kpos_primf.dat");

	// compare prim, primf
	CPPUNIT_ASSERT(Eprim==Eprimf);

	// orthorhombic cell
	{
		const auto RB = findBasis(B,eye<fMat>(3));
		const auto E = calcFoldedBS(H,k.path(),RB,B);

		E.printToFile(prefix+"E_orth.dat");
		k.path().printToFile(prefix+"kpath_orth.dat");
		k.pos().printToFile(prefix+"kpos_orth.dat");
	
		// compare energies at gamma point
		for (const auto e: Eprim.cAt(Gi))
			CPPUNIT_ASSERT(any(E.cAt(Gi).eq(e)));

		// compare number of bands
		const double f = std::abs(det(B)/det(RB));
		CPPUNIT_ASSERT_EQUAL(Eprim.M(),size_t(std::round(f*E.M())));
	}

	// random cell
	{
		fMat R;
		do R = randi<fMat>(3,3,-2,2);
		while (std::abs(det(R))<mtol());
		
		const auto RB = B.prod(R);
		const auto E = calcFoldedBS(H,k.path(),RB,B);

		E.printToFile(prefix+"E_rand.dat");
		k.path().printToFile(prefix+"kpath_rand.dat");
		k.pos().printToFile(prefix+"kpos_rand.dat");
		
		// compare energies at gamma point
		for (const auto e: Eprim.cAt(Gi))
			CPPUNIT_ASSERT(any(E.cAt(Gi).eq(e)));

		// compare number of bands
		const double f = std::abs(det(B)/det(RB));
		CPPUNIT_ASSERT_EQUAL(Eprim.M(),size_t(std::round(f*E.M())));
	}
}


const char* test_fn_bandstructure_calc::test_id() noexcept {
	return "test_fn_bandstructure_calc";
}

CppUnit::Test* test_fn_bandstructure_calc::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_fn_bandstructure_calc>(
		"test_findBandEdges", &test_fn_bandstructure_calc::test_findBandEdges));
	suite->addTest(new CppUnit::TestCaller<test_fn_bandstructure_calc>(
		"test_calcBS", &test_fn_bandstructure_calc::test_calcBS));
	suite->addTest(new CppUnit::TestCaller<test_fn_bandstructure_calc>(
		"test_calcFoldedBS", &test_fn_bandstructure_calc::test_calcFoldedBS));
	
	return suite;
}
