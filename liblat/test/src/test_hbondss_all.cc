// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_hbondss_all.h"
#include "ll_testTools.h"
#include "ll_hbondss.h"
#include "ll_io.h"
#include "aux_io.h"
#include <iostream>


using namespace lm__;
using namespace ll__;
using namespace lm__::test;
using namespace ll__::test;
using namespace aux;

void test_hbondss_all::test_ctor_exceptions() {

	const std::string path = "data/wbhs/";
	ll_hbondss_input inp;
	inp.verbosity = 0;
	
	// bad filename
	inp.wbh = {"blabla"};
	CPPUNIT_ASSERT_THROW(ll_hbondss(inp,std::cout),
		std::invalid_argument);

	// good filenames
	inp.wbh = {path+"ws2_mos2_interf.wad",
		   path+"ws2_mos2_dbll.wad",
		   path+"mos2_ws2_dbll.wad"};
	CPPUNIT_ASSERT_NO_THROW(ll_hbondss(inp,std::cout));
}

void test_hbondss_all::test_all() {

	set_mtol(WTOL__);

	const std::string path = "data/wbhs/";
	
	ll_hbondss_input inp;
	inp.verbosity = 0;
	inp.wbh = {path+"ws2_mos2_interf.wad",
		   path+"ws2_mos2_dbll.wad",
		   path+"mos2_ws2_dbll.wad"};
	inp.sptol = .1;
	inp.perimeter_radius = -5;

	const ll_hbondss W(inp,std::cout);


	// index composing
	{
		const size_t rnd = genRndST(0,NPOS__);
		const size_t j = ll_hbondss::j_(rnd);
		const size_t i = ll_hbondss::i_(rnd);
		CPPUNIT_ASSERT_EQUAL(rnd,ll_hbondss::compInd_(j,i));
	}


	// information
	CPPUNIT_ASSERT_EQUAL(size_t(DIM__),W.dim());
	{
		std::vector<double> R; R.reserve(3);
		for (const auto& w: W)
			R.push_back(w.radius());
		CPPUNIT_ASSERT_DELTA(W.radius(),
			*std::max_element(R.cbegin(),R.cend()),1e-3);
	}
	{
		const fMat r0 = W[0].rmat();
		const fMat r1 = W[1].rmat();
		const fMat r2 = W[2].rmat();
		CPPUNIT_ASSERT_EQUAL(r0&r1&r2,W.rmat());
	}
	{
		const rv r0 = W[0].r();
		const rv r1 = W[1].r();
		const rv r2 = W[2].r();

		rv ref(DIM__,true);
		for (size_t d=0; d!=DIM__; ++d)
			ref[d] = r0[d] && r1[d] && r2[d];

		CPPUNIT_ASSERT(ref==W.r());
	}
	{
		const auto I = W.inds();
		CPPUNIT_ASSERT_EQUAL(I.size(),W[0].cell().N()+
					      W[1].cell().N()+
					      W[2].cell().N());
		for (size_t j=0; j!=W.size(); ++j) {
			auto cI = W[j].inds();
			for (auto& i: cI)
				i = ll_hbondss::compInd_(j,i);

			CPPUNIT_ASSERT(std::all_of(cI.cbegin(),cI.cend(),
			[&I](const size_t i) -> bool {
				return std::find(I.cbegin(),I.cend(),i)!=I.cend();
			}));
		}
	}
	{
		for (const auto& i: W.inds())
			CPPUNIT_ASSERT(W.Norb(i)!=0);

		const auto Norb = W.Norb(W.inds());
		CPPUNIT_ASSERT(std::all_of(Norb.cbegin(),Norb.cend(),
			[](const size_t n)->bool{ return n; }));
	}
	{
		for (size_t i=0; i!=W.size(); ++i) {
			std::vector<size_t> J(W[i].cell().N());
			std::transform(W[i].cell().id().cbegin(),
				       W[i].cell().id().cend(),J.begin(),
				[&W,i](const auto& s) -> size_t {
					return W.ind(std::to_string(i)+':'+s);
			});

			CPPUNIT_ASSERT(std::all_of(J.cbegin(),J.cend(),
				[&W,i](const size_t j) -> bool {
					const auto dd = W.decompInd_(j);
					return dd.j==i && dd.i<W[i].cell().N();
				}));
		}
	}


	// perimeter connections
	{
		set_mtol(WTOL__);

		const auto LM = readOlf(path+"layer_matrix");
		const auto T = W.ind(LM.id);

		// bad type exceptions
		{
			auto BT = T;
			size_t ind = genRndST(0,T.size()-1);
			
			// bad j
			{
				const size_t j = genRndST(W.size(),genRndST(1,10)*W.size());
				BT[ind] = ll_hbondss::compInd_(j,0);
				CPPUNIT_ASSERT_THROW(W.getPerimeterConnections(LM.Ap,BT,inp,std::cout),
					std::invalid_argument);
			}

			// bad i
			{
				const size_t j = genRndST(0,W.size()-1);
				const size_t i = genRndST(W[j].cell().N(),genRndST(1,10)*W[j].cell().N());
				
				BT[ind] = ll_hbondss::compInd_(j,i);
				CPPUNIT_ASSERT_THROW(W.getPerimeterConnections(LM.Ap,BT,inp,std::cout),
					std::invalid_argument);
			}
		}

		// structure incomplete exception
		{
			auto BT = T;
			const size_t j = genRndST(1,W.size()-1);
			for (auto& t: BT)
				if (ll_hbondss::j_(t)==j) t = 0;
			CPPUNIT_ASSERT_THROW(W.getPerimeterConnections(LM.Ap,BT,inp,std::cout),
				std::invalid_argument);
		}

		// check no exception thrown
		CPPUNIT_ASSERT_NO_THROW(W.getPerimeterConnections(LM.Ap,T,inp,std::cout));


		// check properties of extracted bonds
		const auto B = W.getPerimeterConnections(LM.Ap,T,inp,std::cout);

		CPPUNIT_ASSERT_EQUAL(B.RT.size(),B.bnds12.N());

		// bond length
		for (auto i=B.bnds12.ccBegin(),e=B.bnds12.ccEnd(); i!=e; ++i)
			CPPUNIT_ASSERT(norm(*i)<-inp.perimeter_radius);

		// entries are across sections only
		for (auto i=B.RT.cbegin(),e=B.RT.cend(); i!=e; ++i)
			CPPUNIT_ASSERT(i->j1!=i->j2);

		// entries are unique
		{
			auto b = B.bnds12.ccBegin();
			for (auto i=B.RT.cbegin(),e=B.RT.cend(); i!=e; ++i,++b) {
				auto b_ = b+1;
				for (auto i_=i+1; i_!=e; ++i_,++b_)
					CPPUNIT_ASSERT(i->dat!=i_->dat || *b!=*b_);
			}
		}

		// entries are legal
		{
			auto b = B.bnds12.ccBegin();
			for (auto i=B.RT.cbegin(),e=B.RT.cend(); i!=e; ++i,++b) {
				CPPUNIT_ASSERT(i->j1<W.size());
				CPPUNIT_ASSERT(i->j2<W.size());
				CPPUNIT_ASSERT(W[i->j1].getInteraction({i->i1,i->a2}, (*b))!=W.eH());
				CPPUNIT_ASSERT(W[i->j2].getInteraction({i->i2,i->a1},-(*b))!=W.eH());
			}
		}
		

		reset_mtol();
	}

	
	// searching
	{
		// set tolerance
		{
			const double tol = .1;
			auto cW = W; cW.setQueryTolDirect(tol);
			for (const auto& w: cW)
				CPPUNIT_ASSERT_DELTA(tol,w.queryTol(),1e-10);

			cW.setQueryTolCartesian(tol);
			for (const auto& w: cW)
				CPPUNIT_ASSERT_DELTA(w.cell().directTol(tol),w.queryTol(),1e-10);
		}

		// try searching for all bonds in the sub wbh
		{
			auto cW = W; cW.setQueryTolCartesian(.1);

			for (size_t j=0; j!=cW.size(); ++j) {
				const auto S = cW[j].simple();
				
				auto ib = S.ccBegin();
				for (auto ii = S.cbegin(), ie = S.cend(); ii!=ie; ++ii,++ib) {
					const size_t i1 = ll_hbondss::compInd_(j,ii->i1());
					const size_t i2 = ll_hbondss::compInd_(j,ii->i2());
					CPPUNIT_ASSERT(W.getInteraction({i1,i2},*ib)!=W.eH());
				}
			}
		}
	}


	reset_mtol();
}


const char* test_hbondss_all::test_id() noexcept {
	return "test_hbondss_all";
}

CppUnit::Test* test_hbondss_all::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());

	suite->addTest(new CppUnit::TestCaller<test_hbondss_all>(
		"test_ctor_exceptions", &test_hbondss_all::test_ctor_exceptions));
	suite->addTest(new CppUnit::TestCaller<test_hbondss_all>(
		"test_all", &test_hbondss_all::test_all));
	
	return suite;
}
