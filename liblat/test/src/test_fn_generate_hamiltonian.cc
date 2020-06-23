// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_fn_generate_hamiltonian.h"
#include "ll_fn.h"
#include "ll_io.h"
#include "ll_cell.h"
#include "ll_testTools.h"
#include "lm_testTools.h"


using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;


void test_fn_generate_hamiltonian::test_genHam_from_trans() {
	const std::string path = "data/w90/";

	// InAs_mono
	{
		const auto WT = readWannierTransf(path+"inas_mono/wannier90_u.mat",
						  path+"inas_mono/wannier90_u_dis.mat");
		const auto E = readEig(path+"inas_mono/wannier90.eig");
		const auto Href = readHr(path+"inas_mono/wannier90_hr.dat");

		const auto H = genHam(WT,E,Href.R());
		
		CPPUNIT_ASSERT(Href.R()==H.R());
		CPPUNIT_ASSERT(Href.size()==H.size());
		CPPUNIT_ASSERT(std::all_of(H.cbegin(),H.cend(),[&Href](const cMat& i){
			return size(i)==size(Href.front());}));
		
		const double htol = 1e-6;
		for (size_t i=0; i<H.N(); ++i)
			CPPUNIT_ASSERT(all(abs(abs(H[i])-abs(Href[i])).leq(htol)));
	}

	// MoS2_dl
	{
		const auto WT = readWannierTransf(path+"mos2/wannier90_u.mat");
		const auto E = readEig(path+"mos2/wannier90.eig");
		const auto Href = readHr(path+"mos2/wannier90_hr.dat");
		
		const auto H = genHam(WT,E,Href.R());
		
		CPPUNIT_ASSERT(Href.R()==H.R());
		CPPUNIT_ASSERT(Href.size()==H.size());
		CPPUNIT_ASSERT(std::all_of(H.cbegin(),H.cend(),[&Href](const cMat& i){
			return size(i)==size(Href.front());}));
		
		const double htol = 1e-6;
		for (size_t i=0; i<H.N(); ++i)
			CPPUNIT_ASSERT(all(abs(abs(H[i])-abs(Href[i])).leq(htol)));
	}
}

void test_fn_generate_hamiltonian::test_genHam_from_wbh() {
	const std::string path = "data/w90/mos2/";

	const auto hr = readHr(path+"wannier90_hr.dat");

	const ll_hbonds W(path+WBH__);

	const auto ghr = genHam<cMat>(W.cell(),W.r(),W,true);
	
	const auto I = W.inds();
	const auto Norb = W.Norb(I);

	CPPUNIT_ASSERT(!ghr.empty());
	CPPUNIT_ASSERT_EQUAL(W.dim(),ghr.dim());
	CPPUNIT_ASSERT(ghr.R() == round(ghr.R()));

	CPPUNIT_ASSERT(hr.R() == ghr.R());
	
	const auto Ehr = eigh(std::accumulate(hr.cbegin(),hr.cend(),
		zeros<cMat>(hr.Nw(),hr.Nw())));
	const auto Eghr = eigh(std::accumulate(ghr.cbegin(),ghr.cend(),
		zeros<cMat>(ghr.Nw(),ghr.Nw())));
	CPPUNIT_ASSERT(Ehr==Eghr);

	CPPUNIT_ASSERT_EQUAL(hr.Nw(),ghr.Nw());
	CPPUNIT_ASSERT(std::all_of(ghr.cbegin(),ghr.cend(),
		[&hr](const auto& i)->bool{ return size(hr.front())==size(i);}));
	CPPUNIT_ASSERT(std::all_of(ghr.cbegin(),ghr.cbegin(),
		[](const auto& i)->bool{return i!=0.0;}));
}

const char* test_fn_generate_hamiltonian::test_id() noexcept {
	return "test_fn_generate_hamiltonian";
}

CppUnit::Test* test_fn_generate_hamiltonian::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_fn_generate_hamiltonian>(
		"test_genHam_from_trans", &test_fn_generate_hamiltonian::test_genHam_from_trans));
	suite->addTest(new CppUnit::TestCaller<test_fn_generate_hamiltonian>(
		"test_genHam_from_wbh", &test_fn_generate_hamiltonian::test_genHam_from_wbh));
	
	return suite;
}
