// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_bonds_all.h"
#include "ll_testTools.h"
#include "ll_bonds.h"
#include "ll_types.h"
#include "ll_fn.h"
#include "ll_io.h"
#include "aux_io.h"
#include <iostream>


using namespace lm__;
using namespace ll__;
using namespace lm__::test;
using namespace ll__::test;
using namespace aux;

void test_bonds_all::test_i_i_R() {

	const auto i = genRndMEST();

	// i_i
	{
		i_i t1;
		CPPUNIT_ASSERT_EQUAL(NPOS__,t1.i1());
		CPPUNIT_ASSERT_EQUAL(NPOS__,t1.i2());

		i_i t2(i.M);
		CPPUNIT_ASSERT_EQUAL(i.M,t2.i1());
		CPPUNIT_ASSERT_EQUAL(NPOS__,t2.i2());

		i_i t3(i.M,i.N);
		CPPUNIT_ASSERT_EQUAL(i.M,t3.i1());
		CPPUNIT_ASSERT_EQUAL(i.N,t3.i2());

		CPPUNIT_ASSERT(i_i(i.M,i.M)==i_i(i.M,i.M));
		CPPUNIT_ASSERT(i_i(i.M,i.M)<=i_i(i.M,i.M));
		CPPUNIT_ASSERT(i_i(i.M,i.M)>=i_i(i.M,i.M));
		CPPUNIT_ASSERT(i_i(i.M,i.M)!=i_i(i.M,i.N));
		CPPUNIT_ASSERT(i_i(i.M,i.M)!=i_i(i.N,i.M));
		CPPUNIT_ASSERT(i_i(i.M,i.M)!=i_i(i.N,i.N));
		if (i.M<i.N) {
			CPPUNIT_ASSERT(i_i(i.M,i.M)<i_i(i.N,i.M));
			CPPUNIT_ASSERT(i_i(i.M,i.M)<i_i(i.N,i.N));
			CPPUNIT_ASSERT(i_i(i.M,i.N)<i_i(i.N,i.M));
			CPPUNIT_ASSERT(i_i(i.M,i.N)<i_i(i.N,i.N));
		} else {
			CPPUNIT_ASSERT(i_i(i.M,i.M)>i_i(i.N,i.M));
			CPPUNIT_ASSERT(i_i(i.M,i.M)>i_i(i.N,i.N));
			CPPUNIT_ASSERT(i_i(i.M,i.N)>i_i(i.N,i.M));
			CPPUNIT_ASSERT(i_i(i.M,i.N)>i_i(i.N,i.N));
		}
	}

	// i_i_R
	{
		const size_t d = genRndST();
		const auto mat = rand<fMat>(genRndST(),genRndST());

		i_i_R t1;
		CPPUNIT_ASSERT_EQUAL(NPOS__,t1.i1());
		CPPUNIT_ASSERT_EQUAL(NPOS__,t1.i2());
		CPPUNIT_ASSERT_EQUAL(fMat(DIM__,0),t1.R());

		i_i_R t2(i.M,i.N,d);
		CPPUNIT_ASSERT_EQUAL(i.M,t2.i1());
		CPPUNIT_ASSERT_EQUAL(i.N,t2.i2());
		CPPUNIT_ASSERT_EQUAL(fMat(d,0),t2.R());

		i_i_R t3(i.M,i.N,mat);
		CPPUNIT_ASSERT_EQUAL(i.M,t3.i1());
		CPPUNIT_ASSERT_EQUAL(i.N,t3.i2());
		CPPUNIT_ASSERT_EQUAL(mat,t3.R());
	}
}

void test_bonds_all::test_ctor_info_print() {
	
	const auto tCell = genZincblende();

	// empty
	{
		ll_bonds<> b;
		CPPUNIT_ASSERT(b.empty());
	}

	// only cell
	{
		ll_bonds<> b(tCell);
		CPPUNIT_ASSERT(tCell==b.cell());
	}

	// from cell, keep all
	{
		ll_bonds<> tBonds(tCell,NN,[](const fMat& b, const i_i& j)->bool
					{return true;});

		CPPUNIT_ASSERT(tCell==tBonds.cell());
		CPPUNIT_ASSERT_EQUAL(size_t(4),tBonds.size());
		CPPUNIT_ASSERT_EQUAL(4*NN.N(),tBonds.cardinality());
		
		{
			const fMat ref({ .0 ,  .0,  .0,
					 .25, .25, .25,
					-.25,-.25,-.25,
					 .0 ,  .0,  .0},3);
			CPPUNIT_ASSERT(ref.cAt(0)==tBonds.bond(0,0));
			CPPUNIT_ASSERT(ref.cAt(1)==tBonds.bond(0,1));
			CPPUNIT_ASSERT(ref.cAt(2)==tBonds.bond(1,0));
			CPPUNIT_ASSERT(ref.cAt(3)==tBonds.bond(1,1));
		}
		{
			std::vector<i_i> ref = {{0,0},{0,1},{1,0},{1,1}};
			CPPUNIT_ASSERT_EQUAL(ref.size(),tBonds.size());
			
			auto i = ref.cbegin();
			for (auto j=tBonds.cbegin(),e=tBonds.cend(); j!=e; ++j,++i) {
				CPPUNIT_ASSERT(*i == *j);
				CPPUNIT_ASSERT(NN == j->R());
			}
		}

		std::vector<size_t> Iref = {0,1};
		std::vector<size_t> Nref = {tBonds.cardinality(),tBonds.cardinality()};
		CPPUNIT_ASSERT(Iref==tBonds.inds());
		CPPUNIT_ASSERT_EQUAL(Nref[0],tBonds.Nindex(0));
		CPPUNIT_ASSERT_EQUAL(Nref[1],tBonds.Nindex(1));
		CPPUNIT_ASSERT(Nref==tBonds.Nindex(Iref));

		fMat b(tBonds.dim(),0); b.reserve(tBonds.cardinality());
		for (auto i=tCell.ccBegin(),e=tCell.ccEnd(); i!=e; ++i)
			for (auto j=tCell.ccBegin(); j!=e; ++j)
				for (auto n=NN.ccBegin(),ne=NN.ccEnd(); n!=ne; ++n)
					b.push_back(tCell.B().prod(*j+*n-*i));

		CPPUNIT_ASSERT_DELTA(max(mnorm(b)),tBonds.radius(),10.0*mtol());
		for (size_t i=0; i!=b.M(); ++i)
			CPPUNIT_ASSERT_DELTA(nmax(b).mat[i],tBonds.boundary()[0],10.0*mtol());
	}

	// from cell, with cutoff
	{
		const auto r = norm(tCell.B().prod(fMat({.25,.25,.25},tCell.dim())));
		ll_bonds<> tBonds(tCell,NN,r);
		
		CPPUNIT_ASSERT(tCell==tBonds.cell());
	
		CPPUNIT_ASSERT_EQUAL(size_t(2),tBonds.size());
		CPPUNIT_ASSERT_EQUAL(4*tCell.N(),tBonds.cardinality());

		{
			const fMat ref({.25, .25, .25,
				       -.25,-.25,-.25},3);
			CPPUNIT_ASSERT(ref.cAt(0)==tBonds.bond(0,1));
			CPPUNIT_ASSERT(ref.cAt(1)==tBonds.bond(1,0));
		}
		{
			std::vector<i_i> ref = {{0,1},{1,0}};
			CPPUNIT_ASSERT_EQUAL(ref.size(),tBonds.size());
			
			auto i = ref.cbegin();
			for (auto j=tBonds.cbegin(),e=tBonds.cend(); j!=e; ++j,++i) {
				CPPUNIT_ASSERT(*i == *j);
				CPPUNIT_ASSERT(size_t(4) == j->N());
			}
		}

		std::vector<size_t> Iref = {0,1};
		std::vector<size_t> Nref = {8,8};
		CPPUNIT_ASSERT(Iref==tBonds.inds());
		CPPUNIT_ASSERT_EQUAL(Nref[0],tBonds.Nindex(0));
		CPPUNIT_ASSERT_EQUAL(Nref[1],tBonds.Nindex(1));
		CPPUNIT_ASSERT(Nref==tBonds.Nindex(Iref));

		CPPUNIT_ASSERT_DELTA(r,tBonds.radius(),10.0*mtol());
		
		const auto b = tCell.B().prod(fMat({.25,.25,.25},3));
		for (size_t i=0; i!=b.M(); ++i)
			CPPUNIT_ASSERT_DELTA(b[i],tBonds.boundary()[i],10.0*mtol());
	}
}

void test_bonds_all::test_mod_conv() {
	const size_t D = DIM__;
	const auto tCell = genRandom(D);
	const auto tBonds1 = tCell.getBonds(-1.1,NN);

	// rotate
	{
		const auto R = genRndR(D);
		auto tBonds2 = tBonds1;

		tBonds2.rotate(R);
		CPPUNIT_ASSERT(tCell.copy().rotate(R) == tBonds2.cell());
	}

	// scale
	{
		const double f = genRndDouble(.1,2.0);
		auto tBonds2 = tBonds1;
		tBonds2.scale(f);
		CPPUNIT_ASSERT(tCell.copy().scale(f)==tBonds2.cell());
	}
	{
		const fMat f = rand<fMat>(3,1,.1,2.0);
		auto tBonds2 = tBonds1;
		tBonds2.scale(f);
		CPPUNIT_ASSERT(tCell.copy().scale(f)==tBonds2.cell());
	}

	// stress
	{
		const fMat S = eye<fMat>(3,3)-rand<fMat>(3,3,-.05,.05);
		auto tBonds2 = tBonds1;
		tBonds2.stress(S);
		CPPUNIT_ASSERT(tCell.copy().stress(S)==tBonds2.cell());
	}

	// bonds has same size, is sorted along vectors
	const auto b = tBonds1.simple();
	CPPUNIT_ASSERT_EQUAL(b.size(),tBonds1.cardinality());
	CPPUNIT_ASSERT(std::is_sorted(b.ccBegin(),b.ccEnd(),vcmp));

	// bonds has inversion symmetry
	for (auto fi=b.ccBegin(), fe=b.ccEnd(), ri=b.ccEnd()-1; fi!=fe; ++fi,--ri)
		CPPUNIT_ASSERT(*fi==-*ri);
	
	// unique indices are the same
	std::vector<size_t> Iref = tBonds1.inds();
	std::vector<size_t> Ick; Ick.reserve(2*b.size());
	for (const auto& i: b)
		Ick.push_back(i.i1()), Ick.push_back(i.i2());
	std::sort(Ick.begin(),Ick.end());
	Ick.resize(std::distance(Ick.begin(),std::unique(Ick.begin(),Ick.end())));
	CPPUNIT_ASSERT(Iref==Ick);

	// check all vectors in bonds are legal bonds for tCell
	auto j=b.ccBegin();
	for (auto i=b.cbegin(),ie=b.cend(); i!=ie; ++i,++j) {
		const auto p = tCell.cAt(i->i1())+tCell.B().leftDivide(*j);
		const auto itr = std::find(tCell.ccBegin(),tCell.ccEnd(),p%1.0);
		CPPUNIT_ASSERT_EQUAL(size_t(itr),i->i2());
	}
}

void test_bonds_all::test_getBondCenters() {

	const size_t D = DIM__;
	const auto NN = genNNmat(rv(D,false));

	// check ordered id lambda
	const std::regex r1("[_]+"), r2("[:()]+");
	const auto ordId = [&r1,&r2](const idv& id)->bool {
		for (const auto& s: id) {
			std::vector<std::string> buff;
			for (std::sregex_token_iterator i(s.cbegin(),s.cend(),r1,-1), e; i!=e; ++i)
				if (!i->str().empty()) buff.push_back(i->str());
			if (!std::is_sorted(buff.cbegin(),buff.cend()))
				return false;

			for (const auto& g: buff) {
				std::vector<std::string> buff;
				for (std::sregex_token_iterator i(g.cbegin(),g.cend(),r2,-1), e; i!=e; ++i)
					if (!i->str().empty()) buff.push_back(i->str());
				if (buff.size()!=2) return false;
				if (!std::is_sorted(buff.cbegin(),buff.cend())) return false;
			}
		}
		return true;
	};
		
	// empty cell has no bond centers
	{
		auto b = ll_bonds<>(ll_cell()).getBondCenters();
		CPPUNIT_ASSERT(b.empty());
	}

	// check primitive Zincblende
	{
		const auto tCell = genZincblende(D,1.0);
		
		const double f = -1.1;
		
		const auto bnd = tCell.getBonds(f,NN);
		const auto BC = bnd.getBondCenters();
		CPPUNIT_ASSERT(cunique(BC));
		for (const auto t: BC.types())
			CPPUNIT_ASSERT(std::is_sorted(BC.ccBegin(t),BC.ccEnd(t),vcmp));
		
		// check all bond centers included
		const auto sbnd = bnd.simple();
		for (size_t i=0; i!=sbnd.N(); ++i) {
			
			const auto cbc = (tCell.cAt(sbnd[i].i1()) +
				.5*tCell.B().leftDivide(sbnd.cAt(i)))%1.0;

			const auto itr = std::find(BC.ccBegin(),BC.ccEnd(),cbc);
			
			CPPUNIT_ASSERT(BC.ccEnd() != itr); 
			
			size_t i1=sbnd[i].i1(), i2=sbnd[i].i2();
			if (i2<i1) std::swap(i1,i2);
			const std::string cid = "_"+tCell.id(tCell.type(i1))
				+"_"+tCell.id(tCell.type(i2))+"_";
			
			std::string aid = BC.id(BC.type(size_t(itr)));
			CPPUNIT_ASSERT(std::string::npos != aid.find_first_of(cid));
		}
		CPPUNIT_ASSERT(ordId(BC.id()));
	}

	// check expanded Zincblende
	{
		auto C = randi<fMat>(D,D,-2,2);
		do C = randi<fMat>(D,D,-2,2);
		while (std::abs(det(C))<mtol());
		const auto tCell = genZincblende(D,1.0).expand(C);
	
		const double f = -1.1;
		
		const auto bnd = tCell.getBonds(f,NN);
		const auto BC = bnd.getBondCenters();
		CPPUNIT_ASSERT(cunique(BC));
		for (const auto t: BC.types())
			CPPUNIT_ASSERT(std::is_sorted(BC.ccBegin(t),BC.ccEnd(t),vcmp));

		// check all bond centers included
		const auto sbnd = bnd.simple();
		for (size_t i=0; i!=sbnd.N(); ++i) {
			
			const auto cbc = (tCell.cAt(sbnd[i].i1()) +
				.5*tCell.B().leftDivide(sbnd.cAt(i)))%1.0;
			const auto itr = std::find(BC.ccBegin(),BC.ccEnd(),cbc);
			
			CPPUNIT_ASSERT(BC.ccEnd() != itr); 
			
			size_t i1=sbnd[i].i1(), i2=sbnd[i].i2();
			if (i2<i1) std::swap(i1,i2);
			const std::string cid = "_"+tCell.id(tCell.type(i1))
				+"_"+tCell.id(tCell.type(i2))+"_";
			
			std::string aid = BC.id(BC.type(size_t(itr)));
			CPPUNIT_ASSERT(std::string::npos != aid.find_first_of(cid));
		}
		CPPUNIT_ASSERT(ordId(BC.id()));
	}
	
	// check primitive Random
	{
		const auto tCell = genRandom(D);
		
		const double f = -1.1;
		
		const auto bnd = tCell.getBonds(f,NN);
		const auto BC = bnd.getBondCenters();
		CPPUNIT_ASSERT(cunique(BC));
		for (const auto t: BC.types())
			CPPUNIT_ASSERT(std::is_sorted(BC.ccBegin(t),BC.ccEnd(t),vcmp));
		
		// check all bond centers included
		const auto sbnd = bnd.simple();
		for (size_t i=0; i!=sbnd.N(); ++i) {
			
			const auto cbc = (tCell.cAt(sbnd[i].i1()) +
				.5*tCell.B().leftDivide(sbnd.cAt(i)))%1.0;
			const auto itr = std::find(BC.ccBegin(),BC.ccEnd(),cbc);
			
			CPPUNIT_ASSERT(BC.ccEnd() != itr); 
			
			size_t i1=sbnd[i].i1(), i2=sbnd[i].i2();
			if (i2<i1) std::swap(i1,i2);
			const std::string cid = "_"+tCell.id(tCell.type(i1))
				+"_"+tCell.id(tCell.type(i2))+"_";
			
			std::string aid = BC.id(BC.type(size_t(itr)));
			CPPUNIT_ASSERT(std::string::npos != aid.find_first_of(cid));
		}
		CPPUNIT_ASSERT(ordId(BC.id()));
	}

	// check expanded Zincblende
	{
		auto C = randi<fMat>(D,D,-2,2);
		do C = randi<fMat>(D,D,-2,2);
		while (std::abs(det(C))<mtol());
		const auto tCell = genRandom(D).expand(C);
	
		const double f = -1.1;
		
		const auto bnd = tCell.getBonds(f,NN);
		const auto BC = bnd.getBondCenters();
		CPPUNIT_ASSERT(cunique(BC));
		for (const auto t: BC.types())
			CPPUNIT_ASSERT(std::is_sorted(BC.ccBegin(t),BC.ccEnd(t),vcmp));

		// check all bond centers included
		const auto sbnd = bnd.simple();
		for (size_t i=0; i!=sbnd.N(); ++i) {
			
			const auto cbc = (tCell.cAt(sbnd[i].i1()) +
				.5*tCell.B().leftDivide(sbnd.cAt(i)))%1.0;
			const auto itr = std::find(BC.ccBegin(),BC.ccEnd(),cbc);
			
			CPPUNIT_ASSERT(BC.ccEnd() != itr); 
			
			size_t i1=sbnd[i].i1(), i2=sbnd[i].i2();
			if (i2<i1) std::swap(i1,i2);
			const std::string cid = "_"+tCell.id(tCell.type(i1))
				+"_"+tCell.id(tCell.type(i2))+"_";
			
			std::string aid = BC.id(BC.type(size_t(itr)));
			CPPUNIT_ASSERT(std::string::npos != aid.find_first_of(cid));
		}
		CPPUNIT_ASSERT(ordId(BC.id()));
	}

	// check cubic case with duplicate bond centers
	{
		fMat Ap({.0,.0,.0,
			 .5,.0,.0,
			 .0,.5,.0,
			 .0,.0,.5,
			 .5,.5,.0,
			 .5,.0,.5,
			 .0,.5,.5,
			 .5,.5,.5},3);
		const auto tCell = ll_cell(eye<fMat>(3,3),Ap,aCv(Ap.N(),1),idv(Ap.N(),"Cu"));
		
		const double f = -1.1;
		
		const auto bnd = tCell.getBonds(f,NN);
		const auto BC = bnd.getBondCenters();
		CPPUNIT_ASSERT(cunique(BC));
		for (const auto t: BC.types())
			CPPUNIT_ASSERT(std::is_sorted(BC.ccBegin(t),BC.ccEnd(t),vcmp));
		

		// check all bond centers included
		const auto sbnd = bnd.simple();
		for (size_t i=0; i!=sbnd.N(); ++i) {
			
			const auto cbc = (tCell.cAt(sbnd[i].i1()) +
				.5*tCell.B().leftDivide(sbnd.cAt(i)))%1.0;
			const auto itr = std::find(BC.ccBegin(),BC.ccEnd(),cbc);
		
			CPPUNIT_ASSERT(BC.ccEnd() != itr); 
			
			size_t i1=sbnd[i].i1(), i2=sbnd[i].i2();
			if (i2<i1) std::swap(i1,i2);
			const std::string cid = "_"+tCell.id(tCell.type(i1))
				+"_"+tCell.id(tCell.type(i2))+"_";
			
			std::string aid = BC.id(BC.type(size_t(itr)));
			CPPUNIT_ASSERT(std::string::npos != aid.find_first_of(cid));
		}
		CPPUNIT_ASSERT(ordId(BC.id()));
	}
}

void test_bonds_all::test_search() {
	const size_t D = DIM__;
	const auto tCell = genRandom(D);
	const auto tBonds = tCell.getBonds(-1.1,NN);
	const auto b = tBonds.simple();
/*
	// set tolerance
	tBonds.setQueryTolDirect(WTOL__);
	CPPUNIT_ASSERT_EQUAL(double(WTOL__),tBonds.queryTol());

	// try searching for index pairs
	{
		for (size_t i=0; i!=tBonds.size(); ++i) {
			const auto j = tBonds.search(tBonds[i]);
			CPPUNIT_ASSERT(tBonds.cbegin()+i == j);
		}

		const auto I = tBonds.inds();
		const auto j = tBonds.search(i_i(genRndST(I.back(),I.back()+100),
						 genRndST(I.back(),I.back()+100)));
		CPPUNIT_ASSERT(tBonds.cend() == j);
	}

	// try searching for bonds
	{
		auto j = b.ccBegin();
		for (const auto& i: b) {
			const auto itr = tBonds.search(i);
			CPPUNIT_ASSERT(itr != tBonds.cend());

			const auto jtr = tBonds.search(itr,*j);
			CPPUNIT_ASSERT(jtr != tBonds.e());

			const auto cb = tBonds.cell().B().prod(tBonds.bond(*itr) + *jtr);
			CPPUNIT_ASSERT(*j == cb);
			++j;

			const auto jtr_ = tBonds.search(itr,rand<fMat>(D,1,10.0,20.0));
			CPPUNIT_ASSERT_EQUAL(tBonds.e(),jtr_);
		}
	}
*/
	// try searching for slightly randomized bonds
	{
		tBonds.setQueryTolDirect(tBonds.maxtol());

		auto j = b.ccBegin();
		for (const auto& i: b) {
			const auto itr = tBonds.search(i);
			CPPUNIT_ASSERT(itr != tBonds.cend());

			const auto bc = *j + tBonds.cell().B().prod(rand<fMat>(D,1,-.2,.2));
			const auto jtr = tBonds.search(itr,bc);
			CPPUNIT_ASSERT(jtr != tBonds.e());

			const auto cb = tBonds.cell().B().prod(tBonds.bond(*itr) + *jtr);
			CPPUNIT_ASSERT(*j == cb);
			++j;
		}
	}
}

const char* test_bonds_all::test_id() noexcept {
	return "test_bonds_all";
}

CppUnit::Test* test_bonds_all::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
//	suite->addTest(new CppUnit::TestCaller<test_bonds_all>(
//		"test_i_i_R", &test_bonds_all::test_i_i_R));
//	suite->addTest(new CppUnit::TestCaller<test_bonds_all>(
//		"test_ctor_info_print", &test_bonds_all::test_ctor_info_print));
//	suite->addTest(new CppUnit::TestCaller<test_bonds_all>(
//		"test_mod_conv", &test_bonds_all::test_mod_conv));
//	suite->addTest(new CppUnit::TestCaller<test_bonds_all>(
//		"test_getBondCenters", &test_bonds_all::test_getBondCenters));
	suite->addTest(new CppUnit::TestCaller<test_bonds_all>(
		"test_search", &test_bonds_all::test_search));
	
	return suite;
}
