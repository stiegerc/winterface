// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_cell_ctor.h"
#include "ll_cell.h"
#include "ll_types.h"
#include "ll_io.h"
#include "ll_testTools.h"
#include "lm_testTools.h"
#include "aux_io.h"
#include "libmat.h"
#include <iostream>


using namespace lm__;
using namespace ll__;
using namespace lm__::test;
using namespace ll__::test;
using namespace aux;

void test_cell_ctor::test_default() {
	const ll_cell tCell;
	CPPUNIT_ASSERT(tCell.empty());
	CPPUNIT_ASSERT(tCell.B().empty());
	CPPUNIT_ASSERT(tCell.Ap().empty());
	CPPUNIT_ASSERT(!tCell.dim());
	CPPUNIT_ASSERT(!tCell.N());
	CPPUNIT_ASSERT(tCell.id().empty());
}

void test_cell_ctor::test_B_Ap_() {
	const fMat B({0.0,.5,.5,.5,0.0,.5,.5,.5,0.0},3,3);
	const fMat Ap({0.0,0.0,0.0,.25,.25,.25},3,2);

	// with specifying N, no id
	{
		const aCv N = {1,1};
		const ll_cell tCell(B,Ap,N);

		CPPUNIT_ASSERT(!tCell.empty());
		CPPUNIT_ASSERT(!tCell.B().empty());
		CPPUNIT_ASSERT(!tCell.Ap().empty());
		CPPUNIT_ASSERT(B==tCell.B());
		CPPUNIT_ASSERT(Ap==tCell.Ap());
		CPPUNIT_ASSERT(tCell.id().empty());
		CPPUNIT_ASSERT(aTv({0,1})==tCell.types());
		CPPUNIT_ASSERT(aCv({1,1})==tCell.Ntype());
	}
	
	// with specifying N, with id
	{
		const aCv N = {1,1};
		const idv id = {"Pd","S"};
		const ll_cell tCell(B,Ap,N,id);

		CPPUNIT_ASSERT(!tCell.empty());
		CPPUNIT_ASSERT(!tCell.B().empty());
		CPPUNIT_ASSERT(!tCell.Ap().empty());
		CPPUNIT_ASSERT(B==tCell.B());
		CPPUNIT_ASSERT(Ap==tCell.Ap());
		CPPUNIT_ASSERT(aTv({0,1})==tCell.types());
		CPPUNIT_ASSERT(aCv({1,1})==tCell.Ntype());
		CPPUNIT_ASSERT(id==tCell.id());
	}

	// try moving stuff in, no id
	{
		auto tB=B;
		auto tAp=Ap;
		const aCv N = {1,1};

		const double* tB_ptr = tB.data();
		const double* tAp_ptr = tAp.data();
			
		const ll_cell tCell(std::move(tB),std::move(tAp),N);
		CPPUNIT_ASSERT(!tCell.empty());
		CPPUNIT_ASSERT(!tCell.B().empty());
		CPPUNIT_ASSERT(!tCell.Ap().empty());
		CPPUNIT_ASSERT(tCell.id().empty());
		CPPUNIT_ASSERT(B==tCell.B());
		CPPUNIT_ASSERT(Ap==tCell.Ap());
		CPPUNIT_ASSERT_EQUAL(tB_ptr,tCell.B().data());
		CPPUNIT_ASSERT_EQUAL(tAp_ptr,tCell.Ap().data());
	}
	
	// try moving stuff in, with id
	{
		auto tB=B;
		auto tAp=Ap;
		const aCv N = {1,1};
		idv id = {"Br","I"};

		const double* tB_ptr = tB.data();
		const double* tAp_ptr = tAp.data();
		const std::string* id_ptr = id.data();

		const ll_cell tCell(std::move(tB),std::move(tAp),N,std::move(id));
		CPPUNIT_ASSERT(!tCell.empty());
		CPPUNIT_ASSERT(!tCell.B().empty());
		CPPUNIT_ASSERT(!tCell.Ap().empty());
		CPPUNIT_ASSERT(B==tCell.B());
		CPPUNIT_ASSERT(Ap==tCell.Ap());
		CPPUNIT_ASSERT_EQUAL(tB_ptr,tCell.B().data());
		CPPUNIT_ASSERT_EQUAL(tAp_ptr,tCell.Ap().data());
		CPPUNIT_ASSERT_EQUAL(id_ptr,tCell.id().data());
	}

	// try inserting unordered data, with id
	{
		const fMat Ap({0.0,0.0,0.0,
				0.25,0.25,0.25,
				0.0,0.5,0.5,
				0.25,0.75,0.75,
				0.5,0.0,0.5,
				0.75,0.25,0.75,
				0.5,0.5,0.0,
				0.75,0.75,0.25},3,8);
		const idv id = {"Co","Ni","Co","Ni","Co","Ni","Co","Ni"};

		const ll_cell tCell(eye<fMat>(3),Ap,id);
		CPPUNIT_ASSERT(tCell.B()==eye<fMat>(3));
		
		const fMat Ap_s({0.0,0.0,0.0,
				 0.0,0.5,0.5,
				 0.5,0.0,0.5,
				 0.5,0.5,0.0,
				 0.25,0.25,0.25,
				 0.25,0.75,0.75,
				 0.75,0.25,0.75,
				 0.75,0.75,0.25},3,8);
		CPPUNIT_ASSERT(Ap_s==tCell.Ap());
		
		CPPUNIT_ASSERT(aTv({0,1})==tCell.types());
		CPPUNIT_ASSERT(aCv({4,4})==tCell.Ntype());
		CPPUNIT_ASSERT(idv({"Co","Ni"})==tCell.id());
	}

	// from B,Ap,id:  all unique ids
	{
		const auto rnd = genRndCell(3,20,40);

		idv id; id.reserve(rnd.N());
		while (id.size()<id.capacity()) {
			const auto cid = genRndAtomicString(genRndST(1,3));
			if (std::find(id.cbegin(),id.cend(),cid)==id.cend())
				id.push_back(cid);
		}

		ll_cell tCell(rnd.B(),rnd.Ap(),id);
		CPPUNIT_ASSERT(tCell.B()==rnd.B());
		
		const auto NT = tCell.Ntype();
		CPPUNIT_ASSERT_EQUAL(tCell.N(),NT.size());
		CPPUNIT_ASSERT(std::all_of(NT.cbegin(),NT.cend(),
			[](const aC i)->bool{return i==1;}));
		
		auto rAp = rnd.Ap();
		const auto I = aux::sorted_order(id.cbegin(),id.cend());
		aux::reorder(id.begin(),I); aux::reorder(rAp.cBegin(),I);
		
		CPPUNIT_ASSERT(rAp==tCell.Ap());
		
		const auto id_ = tCell.id();
		CPPUNIT_ASSERT(std::equal(id.cbegin(),id.cend(),id_.cbegin()));
	}

	// from B,Ap,id:  all the same id
	{
		const auto rnd = genRndCell(3,20,40);

		idv id(rnd.N(),genRndAtomicString(genRndST(1,3)));

		ll_cell tCell(rnd.B(),rnd.Ap(),id);
		CPPUNIT_ASSERT(tCell.B()==rnd.B());
		for (auto i=tCell.ccBegin(),e=tCell.ccEnd(); i!=e; ++i)
			CPPUNIT_ASSERT(std::find(rnd.ccBegin(),rnd.ccEnd(),*i)!=rnd.ccEnd());

		CPPUNIT_ASSERT(idv({id.front()})==tCell.id());
		CPPUNIT_ASSERT(aTv(tCell.N(),0)==tCell.type());
		CPPUNIT_ASSERT(aCv({tCell.N()})==tCell.Ntype());
	}

	// from B,Ap,id:  some duplicate ids
	{
		const auto rnd = genRndCell(3,20,40);

		idv id; id.reserve(rnd.N());

		const size_t N = genRndST(1,id.capacity());
		while (id.size()<N) {
			const auto cid = genRndAtomicString(genRndST(1,3));
			if (std::find(id.cbegin(),id.cend(),cid)==id.cend())
				id.push_back(cid);
		}
		while (id.size()<id.capacity())
			id.push_back(id[genRndST(0,id.size()-1)]);
		std::random_shuffle(id.begin(),id.end());

		ll_cell tCell(rnd.B(),rnd.Ap(),id);
		CPPUNIT_ASSERT(tCell.B()==rnd.B());
		for (auto i=tCell.ccBegin(),e=tCell.ccEnd(); i!=e; ++i)
			CPPUNIT_ASSERT(std::find(rnd.ccBegin(),rnd.ccEnd(),*i)!=rnd.ccEnd());
		
		auto uid = id;
		std::sort(uid.begin(),uid.end());
		uid.resize(std::distance(uid.begin(),std::unique(uid.begin(),uid.end())));
		CPPUNIT_ASSERT(std::equal(uid.cbegin(),uid.cend(),tCell.id().cbegin()));
	}
}

void test_cell_ctor::test_file() {
	// check exceptions
	{
		CPPUNIT_ASSERT_THROW(ll_cell("xyz"),std::invalid_argument);
	}

	// check POSCAR1
	{
		const ll_cell tCell("data/dft/POSCAR1");
		
		const double a0=3.57;
		const fMat B({0.0,0.5,0.5,
				0.5,0.0,0.5,
				0.5,0.5,0.0},3,3);
			
		const fMat Ap({0.0,0.0,0.0,
				0.25,0.25,0.25},3,2);
			
		CPPUNIT_ASSERT(tCell.B()==(a0*B));
		CPPUNIT_ASSERT(tCell.Ap()==Ap);
		CPPUNIT_ASSERT(aTv({0,1})==tCell.types());
		CPPUNIT_ASSERT(aCv({1,1})==tCell.Ntype());
	}

	// check POSCAR2
	{
		const ll_cell tCell("data/dft/POSCAR2");
			
		const double a0=3.57;
		const fMat B({0.0,0.5,0.5,
				0.5,0.0,0.5,
				0.5,0.5,0.0},3,3);

		const fMat Ap({0.0,0.0,0.0,
				0.25,0.25,0.25},3,2);
			
		CPPUNIT_ASSERT(tCell.B()==(a0*B));
		CPPUNIT_ASSERT(tCell.Ap()==Ap);
		CPPUNIT_ASSERT(aTv({0,1})==tCell.types());
		CPPUNIT_ASSERT(aCv({1,1})==tCell.Ntype());
	}

	// check POSCAR3
	{
		const ll_cell tCell1("data/dft/POSCAR3");
		
		const fMat B({1.7,0.0,2.94449,
				0.0,29.893856,0.0,
				-1.7,0.0,2.94449},3,3);
				
		const fMat Ap({0.0000000000000004, 0.0561988389855091, 0.0000000000000000,
			       0.3333333333333338, 0.2747673635679518, 0.3333333333333335,
			       0.0000000000000004, 0.2185685245824425, 0.0000000000000000,
			       0.0000000000000004, 0.3309662025534611, 0.0000000000000000,
			       0.3333333333333338,-0.0000000000000001, 0.3333333333333335,
			       0.3333333333333338, 0.1123976779710184, 0.3333333333333335},3,6);
			
		CPPUNIT_ASSERT(tCell1.Ap()==Ap);
		CPPUNIT_ASSERT(aTv({0,1})==tCell1.types());
		CPPUNIT_ASSERT(aCv({2,4})==tCell1.Ntype());
		CPPUNIT_ASSERT(idv({"Mo","S"})==tCell1.id());
	}

	// check wannier90 wout
	{
		const ll_cell tCell1("data/w90/mos2/wannier90.wout");

		const fMat B({3.400000,0.000000,0.000000,
			      0.000000,29.893856,0.000000,
			     -1.700000,0.000000,2.944486},3,3);
		const fMat Ap({0.00000,1.68000,0.00000,
			      -0.00000,8.21385,1.96299,
			      -0.00000,0.00000,1.96299,
			       0.00000,9.89385,0.00000,
			       0.00000,6.53387,0.00000,
			      -0.00000,3.36001,1.96299},3,6);
		const idv id = {"Mo","S"};
		const aCv N = {2,4};

		const ll_cell ref(std::move(B),B.leftDivide(Ap)%1.0,std::move(N),std::move(id));
		CPPUNIT_ASSERT(ref==tCell1);
	}

	// check lattice_dat
	{
		const std::string fileName = "outp/gopferteli";
		
		ll_cell tCell1("data/w90/mos2/wannier90.wout");
		fMat C;
		do {
			C = randi<fMat>(tCell1.dim(),tCell1.dim(),-2.0,2.0);
		} while (std::abs(lm__::det(C))<= .9);
		tCell1.expand(C);
		
		printOlf(fileName,tCell1.B(),tCell1.getcAp(),tCell1.id(tCell1.type()),0,.0);

		const ll_cell tCell2(fileName,tCell1.dim());
		CPPUNIT_ASSERT_EQUAL(tCell1,tCell2);
	}
}

void test_cell_ctor::test_copy() {
	const fMat B({0.0,.5,.5,.5,0.0,.5,.5,.5,0.0},3,3);
	const fMat Ap({0.0,0.0,0.0,.25,.25,.25},3,2);
	const aCv T = {1,1};
		
	const ll_cell tCell1(B,Ap,T);
	const ll_cell tCell2=tCell1;
		
	CPPUNIT_ASSERT(!tCell2.empty());
	CPPUNIT_ASSERT(!tCell2.B().empty());
	CPPUNIT_ASSERT(!tCell2.Ap().empty());
	CPPUNIT_ASSERT(tCell2.B()==tCell1.B());
	CPPUNIT_ASSERT(tCell2.Ap()==tCell1.Ap());
	CPPUNIT_ASSERT(tCell1.types()==tCell2.types());
	CPPUNIT_ASSERT(tCell1.Ntype()==tCell2.Ntype());
	CPPUNIT_ASSERT(tCell2.B().begin()!=tCell1.B().begin());
	CPPUNIT_ASSERT(tCell2.Ap().begin()!=tCell1.Ap().begin());
	CPPUNIT_ASSERT(tCell1.id()==tCell2.id());
}

void test_cell_ctor::test_move() {
	const fMat B({0.0,.5,.5,.5,0.0,.5,.5,.5,0.0},3,3);
	const fMat Ap({0.0,0.0,0.0,.25,.25,.25},3,2);
	const aCv N = {1,1};
	const idv id = {"Ta","La"};

	ll_cell tCell1(B,Ap,N,id);
	const ll_cell tCell2=tCell1;
		
	const auto B_ptr = tCell1.B().data();
	const auto Ap_ptr = tCell1.Ap().data();
	const auto id_ptr = tCell1.id().data();

	const ll_cell tCell3 = std::move(tCell1);

	CPPUNIT_ASSERT(!tCell3.empty());
	CPPUNIT_ASSERT(!tCell3.B().empty());
	CPPUNIT_ASSERT(!tCell3.Ap().empty());
	CPPUNIT_ASSERT(tCell3.B()==tCell2.B());
	CPPUNIT_ASSERT(tCell3.Ap()==tCell2.Ap());
	CPPUNIT_ASSERT(tCell2.types()==tCell3.types());
	CPPUNIT_ASSERT(tCell2.Ntype()==tCell3.Ntype());
	CPPUNIT_ASSERT_EQUAL(B_ptr,tCell3.B().data());
	CPPUNIT_ASSERT_EQUAL(Ap_ptr,tCell3.Ap().data());
	CPPUNIT_ASSERT_EQUAL(id_ptr,tCell3.id().data());
}


const char* test_cell_ctor::test_id() noexcept {
	return "test_cell_ctor";
}

CppUnit::Test* test_cell_ctor::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_cell_ctor>(
		"test_default", &test_cell_ctor::test_default));
	suite->addTest(new CppUnit::TestCaller<test_cell_ctor>(
		"test_B_Ap_", &test_cell_ctor::test_B_Ap_));
	suite->addTest(new CppUnit::TestCaller<test_cell_ctor>(
		"test_file", &test_cell_ctor::test_file));
	suite->addTest(new CppUnit::TestCaller<test_cell_ctor>(
		"test_copy", &test_cell_ctor::test_copy));
	suite->addTest(new CppUnit::TestCaller<test_cell_ctor>(
		"test_move", &test_cell_ctor::test_move));
	
	return suite;
}
