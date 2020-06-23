// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_fn_basis_finding.h"
#include "lm_lambda.h"
#include "ll_fn.h"
#include "ll_cell.h"
#include "ll_testTools.h"
#include "lm_testTools.h"


using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;


void test_fn_basis_finding::test_findBasis() {
	// check findBasis using bad input
	{
		const auto tCell = genZincblende(3,1.0);
		
		const auto T = fMat({1.32321321,4.32132312,.213123212,
				     .3232145321,2.321312,32132.321321,
				     321312.32132,1212.2331,321312.121},3,3);
		CPPUNIT_ASSERT_THROW(findBasis(tCell.B(),T),std::runtime_error);
	}
	
	// check findBasis using Zincblende
	{
		const auto tCell = genZincblende(3,1.0);
		
		const size_t imax = 100;
		const auto r = rand<fMat>(3,imax,0.5,1.0);
		for (size_t i=0; i<imax; ++i) {
			const auto R = tCell.B().prod(randi<fMat>(3,3,-5,5)).prod(diag(r.cAt(i)));
			if (ll__::vol(R)<0.25) continue;
		
			CPPUNIT_ASSERT_NO_THROW(findBasis(tCell.B(),R));
			auto CK = findBasis(tCell.B(),R);
			
			// check if all columns of CK are parallel to those in R
			for (size_t j=0; j!=3; ++j)
				CPPUNIT_ASSERT(parallel(CK.cAt(j),R.cAt(j))>0.0);

			tCell.B().leftDivideEq(CK);
			CPPUNIT_ASSERT((CK-round(CK))==0.0);
		}
	}
	
	// check findBasis using diChalcogenide
	{
		const auto tCell = genDiChalcogenide(1.0);
		
		const size_t imax = 100;
		const auto r = rand<fMat>(3,imax,0.5,1.0);
		for (size_t i=0; i<imax; ++i) {
			const auto R = tCell.B().prod(randi<fMat>(3,3,-5,5).prod(diag(r.cAt(i))));
			if (ll__::vol(R)<3.0*std::sqrt(3.0)) continue;
			
			auto CK = findBasis(tCell.B(),R);
			
			// check if all columns of CK are parallel to those in R
			for (size_t j=0; j<3; ++j)
				CPPUNIT_ASSERT(parallel(CK.cAt(j),R.cAt(j))>0.0);
			
			// check if all integers in basis B
			tCell.B().leftDivideEq(CK);
			CPPUNIT_ASSERT((CK-round(CK))==0.0);
		}
	}
	
	// check findBasis using higher dimensional "Zincblende" (for shits and giggles)
	{
		const size_t jmax = 3;
		const size_t imax = 10;
		
		const auto S = randi<fMat>(1,jmax,4,7);
		for (size_t j=0; j<jmax; ++j) {
			const auto s = size_t(S[j]);
			
			const auto tCell = genZincblende(s,1.0);
			const auto r = rand<fMat>(s,imax,0.5,1.0);
			
			for (size_t i=0; i<imax; ++i) {
				const auto R = tCell.B().prod(randi<fMat>(s,s,-5,5)).prod(diag(r.cAt(i)));
				if (ll__::volnorm(R)<tCell.volnorm()) continue;
			
				auto CK = findBasis(tCell.B(),R);
								
				// check if all columns of CK are parallel to those in R
				for (size_t j=0; j<3; ++j)
					CPPUNIT_ASSERT(parallel(CK.cAt(j),R.cAt(j))>0.0);
				
				// check if all integers in basis B
				tCell.B().leftDivideEq(CK);
				
				CPPUNIT_ASSERT((CK-round(CK))==0.0);
			}
		}
	}

	// check findBasis using random cell
	{
		const auto tCell = genRandom();
		
		const size_t imax = 100;
		const auto r = rand<fMat>(3,imax,0.5,1.0);
		for (size_t i=0; i<imax; ++i) {
			const auto R = tCell.B().prod(randi<fMat>(3,3,-5,5).prod(diag(r.cAt(i))));
			if (ll__::vol(R)<3.0*std::sqrt(3.0)) continue;
			
			auto CK = findBasis(tCell.B(),R);
			
			// check if all columns of CK are parallel to those in R
			for (size_t j=0; j<3; ++j)
				CPPUNIT_ASSERT(parallel(CK.cAt(j),R.cAt(j))>0.0);
			
			// check if all integers in basis B
			tCell.B().leftDivideEq(CK);
			CPPUNIT_ASSERT((CK-round(CK))==0.0);
		}
	}
}
void test_fn_basis_finding::test_orthogonalize() {

	const auto tCell = genRandom();
	const auto r = genRndrv(tCell.dim());
	
	const auto OB = orthogonalize(tCell.B(),r);

	const auto i = inds(r), ni = ninds(r);
	if (ni.empty()) CPPUNIT_ASSERT(OB.ob());
	else {
		if (i.empty()) CPPUNIT_ASSERT(tCell.B()==OB);
		else CPPUNIT_ASSERT(zeros<fMat>(ni.size(),i.size())==
				    T(OB.get({},ni)).prod(OB.get({},i)));
		CPPUNIT_ASSERT(tCell.B().get({},ni)==OB.get({},ni));
	}
	CPPUNIT_ASSERT(mnorm(tCell.B())==mnorm(OB));
	if (tCell.B().ob()) CPPUNIT_ASSERT(tCell.B()==OB);
}

const char* test_fn_basis_finding::test_id() noexcept {
	return "test_fn_basis_finding";
}

CppUnit::Test* test_fn_basis_finding::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_fn_basis_finding>(
		"test_findBasis", &test_fn_basis_finding::test_findBasis));
	suite->addTest(new CppUnit::TestCaller<test_fn_basis_finding>(
		"test_orthogonalize", &test_fn_basis_finding::test_orthogonalize));
	
	return suite;
}
