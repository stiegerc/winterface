// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_cell_conversion.h"
#include "lm_testTools.h"
#include "ll_testTools.h"
#include "aux_io.h"
#include "ll_fn.h"
#include "ll_io.h"
#include "ll_lambda.h"
#include <iostream>


using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;
using namespace aux;


void test_cell_conversion::test_getRB() {
	const auto D =  genRndST(1,5);

	const auto tCell1 = genRandom(D);
	const auto RB = tCell1.getRB();

	CPPUNIT_ASSERT(size(RB)==size(tCell1.B()));
	CPPUNIT_ASSERT(T(RB).prod(tCell1.B())==2.0*M_PI*eye<fMat>(D));
}

void test_cell_conversion::test_getAp() {
	const auto D =  genRndST(1,5);
	const auto tCell1 = genRandom(D);

	// invalid type gives empty
	CPPUNIT_ASSERT(tCell1.getAp(NPOS__).empty());

	// valid types
	for (const auto t: tCell1.types()) {
		const auto Ap = tCell1.getAp(t);
		CPPUNIT_ASSERT_EQUAL(tCell1.dim(),Ap.M());
		CPPUNIT_ASSERT_EQUAL(tCell1.Ntype(t),Ap.N());
		
		for (auto i=Ap.ccBegin(),e=Ap.ccEnd(); i!=e; ++i) {
			const auto j = std::find(tCell1.Ap().ccBegin(),tCell1.Ap().ccEnd(),*i);
			CPPUNIT_ASSERT_EQUAL(t,tCell1.type(j.i()));
		}
	}

	// random aTv
	{
		const auto ts = genRndaTv(tCell1);
		const auto N = tCell1.Ntype(ts);
		const auto NN = std::accumulate(N.cbegin(),N.cend(),aC(0));
		const auto Ap = tCell1.getAp(ts);
		
		CPPUNIT_ASSERT_EQUAL(NN,Ap.N());
		
		auto i = Ap.ccBegin();
		for (const auto t: ts) {
			const auto Ap_ = tCell1.getAp(t);
			for (auto j=Ap_.ccBegin(),e=Ap_.ccEnd(); j!=e; ++j,++i)
				CPPUNIT_ASSERT(*i==*j);
		}
			
	}
}

void test_cell_conversion::test_getcAp() {
	const auto D =  genRndST(1,5);
	const auto tCell1 = genRandom(D);
	
	
	// valid types
	for (const auto t: tCell1.types()) {
		const auto Ap = tCell1.getAp(t);
		const auto cAp = tCell1.getcAp(t);
		CPPUNIT_ASSERT(tCell1.B().prod(Ap)==cAp);
	}

	// random aTv
	{
		const auto ts = genRndaTv(tCell1);
		const auto Ap = tCell1.getAp(ts);
		const auto cAp = tCell1.getcAp(ts);
		CPPUNIT_ASSERT(tCell1.B().prod(Ap)==cAp);
	}	
}

void test_cell_conversion::test_getSubCell() {
	const auto D =  genRndST(1,5);
	const auto tCell1 = genRandom(D);
	
	// try one type
	for (const auto t: tCell1.types()) {
		const auto sCell = tCell1.getSubCell({t});

		CPPUNIT_ASSERT(tCell1.B()==sCell.B());
		CPPUNIT_ASSERT(tCell1.getAp(t)==sCell.Ap());
		CPPUNIT_ASSERT_EQUAL(size_t(1),sCell.types().size());
		CPPUNIT_ASSERT_EQUAL(aT(0),sCell.types().front());
		CPPUNIT_ASSERT_EQUAL(tCell1.Ntype(t),sCell.Ntype(aT(0)));
		CPPUNIT_ASSERT_EQUAL(size_t(1),sCell.id().size());
		
		CPPUNIT_ASSERT_EQUAL(tCell1.stripId(tCell1.id(t)),
				     sCell.stripId(sCell.id(0)));
	}
	
	// random types
	{
		const auto ts = genRndaTv(tCell1);
		const auto sCell = tCell1.getSubCell(ts);
		
		CPPUNIT_ASSERT(tCell1.B()==sCell.B());
		CPPUNIT_ASSERT(tCell1.getAp(ts)==sCell.Ap());
		CPPUNIT_ASSERT(rg(0,1,ts.size())==sCell.types());
		CPPUNIT_ASSERT(tCell1.Ntype(ts)==sCell.Ntype());

		for (size_t t=0; t!=ts.size(); ++t)
			CPPUNIT_ASSERT_EQUAL(tCell1.stripId(tCell1.id(ts[t])),
					     sCell.stripId(sCell.id(t)));
	}
}

void test_cell_conversion::test_getBonds() {
	
	const size_t D = DIM__;
	const auto NN = genNNmat(rv(D,false));
		
	// empty cell has no bonds
	{
		const double f = -1.1;
		ll_cell tCell;

		const auto b = tCell.getBonds(f,NN);
		CPPUNIT_ASSERT(b.empty());
		CPPUNIT_ASSERT(tCell==b.cell());
		CPPUNIT_ASSERT_EQUAL(double(0.0),b.radius());
	}

	// check primitive Zincblende
	{
		auto tCell = genZincblende(D,1.0);
		const double f = -1.1;
			
		const double d = sqrt(3.0)/4.0;
		const auto b = tCell.getBonds(f,NN);

		// check cell is the same
		CPPUNIT_ASSERT(tCell==b.cell());

		// check size is 2 (#index pairs)
		CPPUNIT_ASSERT_EQUAL(size_t(2),b.size());

		// check index pairs are (0:1) and (1:0)
		CPPUNIT_ASSERT_EQUAL(size_t(0),b[0].i1());
		CPPUNIT_ASSERT_EQUAL(size_t(1),b[0].i2());
		CPPUNIT_ASSERT_EQUAL(size_t(1),b[1].i1());
		CPPUNIT_ASSERT_EQUAL(size_t(0),b[1].i2());

		// check inversion symmetry
		CPPUNIT_ASSERT(b.bond(0,1)+b.bond(1,0)==0.0);
		CPPUNIT_ASSERT(std::all_of(b[0].ccBegin(),b[0].ccEnd(),
			[&b](const auto& i)->bool{
				const auto itr = std::lower_bound(b[1].ccBegin(),b[1].ccEnd(),-i,vcmp);
				return itr!=b[1].ccEnd() && *itr==-i;
			})
		);
		
		// check bonds all have same length
		{
			const auto bb = b.simple();
			CPPUNIT_ASSERT(mnorm(bb.bonds())==d);
		}

		// check each index has 4 neighbours
		for (const auto i: b.inds())
			CPPUNIT_ASSERT_EQUAL(size_t(4),b.Nindex(i)/2);
	}

	// compare primitive to cubic Zincblende
	{
		auto tCell1 = genZincblende(D,1.0);
		auto tCell2 = tCell1.copy().changeBasis(eye<fMat>(3));
		size_t F = std::round(tCell2.vol()/tCell1.vol());
			
		const double f = -1.1;
		const auto b1 = tCell1.getBonds(f,NN);
		const auto b2 = tCell2.getBonds(f,NN);

		// check #indices scales with volume ratio F
		const auto I1 = b1.inds();
		const auto I2 = b2.inds();
		CPPUNIT_ASSERT_EQUAL(F*I1.size(),I2.size());

		// check each index has 4 neighbours
		for (const auto i: I2)
			CPPUNIT_ASSERT_EQUAL(size_t(4),b2.Nindex(i)/2);

		// compare bonds, each bond vector in bb1 should appear F times in bb2
		const auto bb1 = b1.simple();
		const auto bb2 = b2.simple();
		CPPUNIT_ASSERT_EQUAL(F*bb1.size(),bb2.size());
		for (auto i=bb1.ccBegin(),ie=bb1.ccEnd(); i!=ie; ++i) {
			for (auto j=std::lower_bound(bb2.ccBegin(),bb2.ccEnd(),*i,vcmp),je=j+F;
					j!=je; ++j) {
				CPPUNIT_ASSERT(*i==*j);
				CPPUNIT_ASSERT_EQUAL(bb1[size_t(i)],
					i_i(tCell2.type(bb2[size_t(j)].i1()),
					    tCell2.type(bb2[size_t(j)].i2())));
			}
		}
	}

	// check slightly randomized Zincblende
	{
		const auto B = .5*(ones<fMat>(3)-eye<fMat>(3));
			
		fMat Ap(D,0); Ap.reserve(2);
		Ap.push_back(zeros<fMat>(3,1));
		Ap.push_back(rand<fMat>(3,1,.24,.26));
			
		ll_cell tCell(std::move(B),std::move(Ap),{1,1});
		const double f = -1.1;
		
		const double d = norm(tCell.B().prod(fMat({.26,.26,.26})));
		const auto b = tCell.getBonds(f,NN);

		// check cell is the same
		CPPUNIT_ASSERT(tCell==b.cell());

		// check size is 2 (#index pairs)
		CPPUNIT_ASSERT_EQUAL(size_t(2),b.size());

		// check index pairs are (0:1) and (1:0)
		CPPUNIT_ASSERT_EQUAL(size_t(0),b[0].i1());
		CPPUNIT_ASSERT_EQUAL(size_t(1),b[0].i2());
		CPPUNIT_ASSERT_EQUAL(size_t(1),b[1].i1());
		CPPUNIT_ASSERT_EQUAL(size_t(0),b[1].i2());

		// check inversion symmetry
		CPPUNIT_ASSERT(b.bond(0,1)+b.bond(1,0)==0.0);
		CPPUNIT_ASSERT(std::all_of(b[0].ccBegin(),b[0].ccEnd(),
			[&b](const auto& i)->bool{
				const auto itr = std::lower_bound(b[1].ccBegin(),b[1].ccEnd(),-i,vcmp);
				return itr!=b[1].ccEnd() && *itr==-i;
			})
		);

		// check bonds all have length <= d
		{
			const auto bb = b.simple();
			CPPUNIT_ASSERT(all(mnorm(bb.bonds()).leq(d)));
		}

		// check each index has 4 neighbours
		for (const auto i: b.inds())
			CPPUNIT_ASSERT_EQUAL(size_t(4),b.Nindex(i)/2);
	}
}


const char* test_cell_conversion::test_id() noexcept {
	return "test_cell_conversion";
}

CppUnit::Test* test_cell_conversion::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_cell_conversion>(
		"test_getRB", &test_cell_conversion::test_getRB));
	suite->addTest(new CppUnit::TestCaller<test_cell_conversion>(
		"test_getAp", &test_cell_conversion::test_getAp));
	suite->addTest(new CppUnit::TestCaller<test_cell_conversion>(
		"test_getcAp", &test_cell_conversion::test_getcAp));
	suite->addTest(new CppUnit::TestCaller<test_cell_conversion>(
		"test_getSubCell", &test_cell_conversion::test_getSubCell));
	suite->addTest(new CppUnit::TestCaller<test_cell_conversion>(
		"test_getBonds", &test_cell_conversion::test_getBonds));
	
	return suite;
}
