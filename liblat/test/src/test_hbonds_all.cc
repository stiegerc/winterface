// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_hbonds_all.h"
#include "ll_testTools.h"
#include "ll_hbonds.h"
#include "ll_types.h"
#include "ll_io.h"
#include "ll_fn.h"
#include "aux_io.h"
#include <iostream>


using namespace lm__;
using namespace ll__;
using namespace lm__::test;
using namespace ll__::test;
using namespace aux;

void test_hbonds_all::test_ctor_exceptions() {
	const std::string path = "data/w90/mos2/";

	// bad input files
	{
		ll_wmatching_input inp;
		inp.wout = "blub";
		inp.hrdat = path+"wannier90_hr.dat";
		inp.mode = {"atomic"};
		inp.tol = 0.0;
		inp.rm_unmatched = 1;
		inp.verbosity = 0;
		
		CPPUNIT_ASSERT_THROW(ll_hbonds(inp,std::cout),
			std::invalid_argument);
		
		inp.wout = path+"wannier90.wout";
		inp.hrdat = "blub";
		
		CPPUNIT_ASSERT_THROW(ll_hbonds(inp,std::cout),
			std::invalid_argument);
	}
}

void test_hbonds_all::test_all() {
	const std::string path = "data/w90/mos2/";
	
	// default
	{
		ll_hbonds W;
		CPPUNIT_ASSERT(W.empty());
		CPPUNIT_ASSERT(W.eH()==cMat());
	}

	// from input file
	ll_wmatching_input inp;
	inp.wout = path+"wannier90.wout";
	inp.hrdat = path+"wannier90_hr.dat";
	inp.mode = {"atomic"};
	inp.tol = 0.0;
	inp.rm_unmatched = 1;
	inp.verbosity = 0;

	const auto hr = readHr(inp.hrdat);
	
	const auto tol = genRndDouble(1e-4,.1);
	const auto keep = [tol](const cMat& inp)->bool{
		for (const auto i: inp)
			if (cmph(i,tol)) return true;
		return false;
	};
	const auto hr_tol = readHr(path+"wannier90_hr.dat",{},keep);

	const auto rref = hrDim(inp.hrdat);

	// test with everything kept
	const ll_hbonds W(inp,std::cout);
	{
		CPPUNIT_ASSERT_EQUAL(W.cell().N()*W.cell().N()*hr.N(),W.cardinality());
		for (const auto& i: W.inds())
			CPPUNIT_ASSERT_EQUAL(2*W.cell().N()*hr.N(),W.Nindex(i));

		CPPUNIT_ASSERT(std::all_of(W.cbegin(),W.cend(),
			[&hr](const auto& inp)->bool{return inp.R()==hr.R();}));

		// check bonds are symmetric
		const auto b = W.simple();
		for (auto i=b.ccBegin(),j=b.ccEnd()-1,ie=b.ccEnd(); i!=ie; ++i,--j)
			CPPUNIT_ASSERT(*i == -(*j));

		// check all interactions are found
		{
			auto j = b.cbegin();
			for (auto i=b.ccBegin(),ie=b.ccEnd(); i!=ie; ++i,++j)
				CPPUNIT_ASSERT(W.getInteraction(*j,*i) != W.eH());
		}

		// check interactions match spacial symmetry
		{
			auto j = b.cbegin();
			for (auto i=b.ccBegin(),ie=b.ccEnd(); i!=ie; ++i,++j)
				CPPUNIT_ASSERT(W.getInteraction(*j,*i) ==
					     T(W.getInteraction(conj(*j),-(*i))));
		}

		// check total number of hamiltonian elements in wbh is same as in hr
		CPPUNIT_ASSERT_EQUAL(hr.Nw()*hr.Nw()*hr.N(),
			std::accumulate(W.cbegin(),W.cend(),size_t(0),
				[](const size_t acc, const auto& inp)->size_t{
					return acc +
						std::accumulate(inp.cbegin(),inp.cend(),size_t(0),
						[](const size_t acc, const auto& inp)->size_t{
							return acc+inp.size();
						});
				})
		);
			
		// check all data of wannier90 hamiltonian is found somewhere in bonds hamiltonian
		for (const auto& H: hr)
		for (const auto& h: H)
			CPPUNIT_ASSERT(
			std::any_of(W.cbegin(),W.cend(),[&h](const auto& i)->bool {
				return std::any_of(i.cbegin(),i.cend(),[&h](const auto& ch)->bool {
					return std::any_of(ch.cbegin(),ch.cend(),
							[&h](const auto& hh)->bool {
						return (h==hh);
					});
				});
			}));

		// test r
		CPPUNIT_ASSERT(rref==W.r());
	}

	// test with cutoff tolerance
	inp.tol = tol;
	const ll_hbonds W_(inp,std::cout);
	{
		CPPUNIT_ASSERT(W.cell().N()*W.cell().N()*hr.N()>=W_.cardinality());
		for (const auto& i: W_.inds())
			CPPUNIT_ASSERT(2*W.cell().N()*hr.N()>=W_.Nindex(i));

		// check all Rvecs in wbh_ are also found in hr
		CPPUNIT_ASSERT(std::all_of(W_.cbegin(),W_.cend(),
			[&hr](const auto& inp)->bool{
				return inp.N()<=hr.N() &&
					std::all_of(inp.ccBegin(),inp.ccEnd(),[&hr](const auto& cr)->bool{
					const auto itr = std::lower_bound(hr.ccBegin(),hr.ccEnd(),cr,vcmp);
					return itr!=hr.ccEnd() && *itr==cr;
					});
			;}));
		
		// check bonds are symmetric
		const auto b = W_.simple();
		for (auto i=b.ccBegin(),j=b.ccEnd()-1,ie=b.ccEnd(); i!=ie; ++i,--j)
			CPPUNIT_ASSERT(*i == -(*j));

		// check all interactions are found
		{
			auto j = b.cbegin();
			for (auto i=b.ccBegin(),ie=b.ccEnd(); i!=ie; ++i,++j)
				CPPUNIT_ASSERT(W_.getInteraction(*j,*i) != W_.eH());
		}

		// check interactions match spacial symmetry
		{
			auto j = b.cbegin();
			for (auto i=b.ccBegin(),ie=b.ccEnd(); i!=ie; ++i,++j)
				CPPUNIT_ASSERT(W_.getInteraction(*j,*i) ==
					     T(W_.getInteraction(conj(*j),-(*i))));
		}
		
		// test r
		CPPUNIT_ASSERT(rref==W_.r());
	}

	// compare full vs cutoff
	{
		set_mtol(W_.maxtol());

		const auto b = W_.simple();

		// check all bonds/interactions from wbh_ are also found in wbh
		{
			auto j = b.cbegin();
			for (auto i=b.ccBegin(),ie=b.ccEnd(); i!=ie; ++i,++j)
				CPPUNIT_ASSERT(W_.getInteraction(*j,*i) ==
					       W.getInteraction(*j,*i));
		}
		reset_mtol();	
	}

	// test accepting none
	{
		inp.tol = 1e12;
		const ll_hbonds BHempty(inp,std::cout);
		CPPUNIT_ASSERT(BHempty.empty());
	}

	// try writing to disk, then reading back in and compare and check exceptions
	{
		const std::string prefix = "outp/";
		W.writeToFile(prefix+WBH__);

		const ll_hbonds W_read(prefix+WBH__);
		CPPUNIT_ASSERT_EQUAL(W.size(),W_read.size());
		CPPUNIT_ASSERT_EQUAL(W.cardinality(),W_read.cardinality());

		CPPUNIT_ASSERT(W.Norb()==W_read.Norb());

		CPPUNIT_ASSERT(W.cell()==W_read.cell());
		for (size_t i=0; i!=W.size(); ++i)
			CPPUNIT_ASSERT_EQUAL(W[i],W_read[i]);

		// bad file exceptions
		CPPUNIT_ASSERT_THROW(ll_hbonds(prefix+"hyugsa.wad"),std::invalid_argument);
		CPPUNIT_ASSERT_THROW(W.writeToFile(prefix+"/not_exist/BH.wad"),std::invalid_argument);

		// bad header exception
		CPPUNIT_ASSERT_THROW(ll_hbonds(path+"bad_header.wad"),std::runtime_error);
	}
}


const char* test_hbonds_all::test_id() noexcept {
	return "test_hbonds_all";
}

CppUnit::Test* test_hbonds_all::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());

	suite->addTest(new CppUnit::TestCaller<test_hbonds_all>(
		"test_ctor_exceptions", &test_hbonds_all::test_ctor_exceptions));
	suite->addTest(new CppUnit::TestCaller<test_hbonds_all>(
		"test_all", &test_hbonds_all::test_all));
	
	return suite;
}
