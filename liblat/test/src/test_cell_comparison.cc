// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_cell_comparison.h"
#include "lm_testTools.h"
#include "ll_testTools.h"
#include "aux_io.h"
#include <iostream>


using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;
using namespace aux;

void test_cell_comparison::test_getPmat() {
	const auto D = genRndST(1,5);
	
	// some empty
	{
		const cell tCell1;
		const cell tCell2;

		CPPUNIT_ASSERT(tCell1.getPmat(tCell2).empty());

		const auto tCell3 = genRandom(D);
		CPPUNIT_ASSERT(tCell1.getPmat(tCell3).empty());
		CPPUNIT_ASSERT(tCell3.getPmat(tCell1).empty());
	}
	
	// dim mismatch
	{
		const auto tCell1 = genRandom(D);
		const auto tCell2 = genRandom(D+1);

		CPPUNIT_ASSERT(tCell1.getPmat(tCell2).empty());
	}
	
	{
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;

		const auto n1 = genRndST(0,D-1);
		const auto n2 = genRndST(0,D-1);

		tCell1.swapDim(n1,n2);
		tCell1.permute(tCell1.getPmat(tCell2));

		CPPUNIT_ASSERT(tCell1==tCell2);
	}
	{
		const auto C = mut_excl_cells(D);
		CPPUNIT_ASSERT(C.cell1.getPmat(C.cell2).empty());
	}
}

void test_cell_comparison::test_getAvec() {
	const auto D = genRndST(1,5);

	// some empty
	{
		const cell tCell1;
		const cell tCell2;
		
		CPPUNIT_ASSERT(tCell1.getAvec(tCell2).empty());
		const auto tCell3 = genRandom(D);
		CPPUNIT_ASSERT(tCell1.getAvec(tCell3).empty());
		CPPUNIT_ASSERT(tCell3.getAvec(tCell1).empty());
	}

	// dim mismatches
	{
		const auto tCell1 = genRandom(D);
		const auto tCell2 = genRandom(D+1);


		CPPUNIT_ASSERT(tCell1.getAvec(tCell1,eye<fMat>(D+1)).empty());
		CPPUNIT_ASSERT(tCell1.getAvec(tCell2,eye<fMat>(D)).empty());
		CPPUNIT_ASSERT(tCell1.getAvec(tCell2).empty());
	}

	// no shift
	{
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;
		
		const auto A = tCell1.getAvec(tCell2);
		CPPUNIT_ASSERT(A==0.0);
	}

	// with shift
	{	
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;
		
		const auto sh = rand<fMat>(D,1,.2,.8);
		tCell1.shift(sh);
		
		const auto A = tCell1.getAvec(tCell2);
		CPPUNIT_ASSERT(!A.empty());
		
		tCell1.shift(A);
		CPPUNIT_ASSERT(tCell2==tCell1);
	}
	
	// with permutation
	{
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;
		
		const auto n1 = genRndST(0,D-1);
		const auto n2 = genRndST(0,D-1);

		tCell1.swapDim(n1,n2);
		
		const auto A = tCell1.getAvec(tCell2,tCell1.getPmat(tCell2));
		CPPUNIT_ASSERT(A==0.0);
	}

	// with permutation and shift
	{
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;
		
		const auto sh = rand<fMat>(D,1,.2,.8);
		tCell1.shift(sh);
		
		const auto n1 = genRndST(0,D-1);
		const auto n2 = genRndST(0,D-1);

		tCell1.swapDim(n1,n2);
		
		const auto P = tCell1.getPmat(tCell2);
		const auto A = tCell1.getAvec(tCell2,P);

		tCell1.shift(A);
		tCell1.permute(P);
	
		CPPUNIT_ASSERT(tCell2==tCell1);
	}

	// no Avec for foreigners
	{
		const auto C = mut_excl_cells(D);
		
		CPPUNIT_ASSERT(C.cell1.getAvec(C.cell2).empty());
	}
}

void test_cell_comparison::test_sameLattice() {
	const auto D = genRndST(1,5);

	// some empty
	{
		const cell tCell1;
		const cell tCell2;

		CPPUNIT_ASSERT(tCell1.sameLattice(tCell2));

		const auto tCell3 = genRandom(D);
		CPPUNIT_ASSERT(!tCell1.sameLattice(tCell3));
		CPPUNIT_ASSERT(!tCell3.sameLattice(tCell1));
	}

	// dim mismatch
	{
		const auto tCell1 = genRandom(D);
		const auto tCell2 = genRandom(D+1);

		CPPUNIT_ASSERT(!tCell1.sameLattice(tCell2));
	}

	// same cell
	{
		const auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;

		CPPUNIT_ASSERT(tCell1.sameLattice(tCell2));
	}
	
	// foreigners
	{
		const auto C = mut_excl_cells(D);
		CPPUNIT_ASSERT(!C.cell1.sameLattice(C.cell2));;
	}
	
	// with shift
	{	
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;
		
		const auto sh = rand<fMat>(D,1,.2,.8);
		tCell1.shift(sh);
		
		CPPUNIT_ASSERT(tCell1.sameLattice(tCell2));
	}

	// with permutation
	{
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;
		
		const auto n1 = genRndST(0,D-1);
		const auto n2 = genRndST(0,D-1);

		tCell1.swapDim(n1,n2);
		
		CPPUNIT_ASSERT(tCell1.sameLattice(tCell2));
	}

	// with shift and permutation
	{	
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;
		
		const auto sh = rand<fMat>(D,1,.2,.8);
		tCell1.shift(sh);
		
		const auto n1 = genRndST(0,D-1);
		const auto n2 = genRndST(0,D-1);

		tCell1.swapDim(n1,n2);
		
		CPPUNIT_ASSERT(tCell1.sameLattice(tCell2));
	}
}

void test_cell_comparison::test_operator_equal_unequal() {
	const auto D = genRndST(1,5);

	// same cell
	{
		const auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;

		CPPUNIT_ASSERT(tCell1==tCell2);
		CPPUNIT_ASSERT(!(tCell1!=tCell2));
	}

	// foreigners
	{
		const auto C = mut_excl_cells(D);
		CPPUNIT_ASSERT(!(C.cell1==C.cell2));
		CPPUNIT_ASSERT(C.cell1!=C.cell2);
	}

	// with shift
	{	
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;
		
		const auto sh = rand<fMat>(D,1,.2,.8);
		tCell1.shift(sh);
		
		CPPUNIT_ASSERT(!(tCell1==tCell2));
		CPPUNIT_ASSERT(tCell1!=tCell2);
	}

	// with permutation
	{
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;
	

		const auto n1 = genRndST(0,D-1);
		const auto n2 = genRndST(0,D-1);

		tCell1.swapDim(n1,n2);
		
		if (n1!=n2) {
			CPPUNIT_ASSERT(!(tCell1==tCell2));
			CPPUNIT_ASSERT(tCell1!=tCell2);
		} else {
			CPPUNIT_ASSERT(tCell1==tCell2);
			CPPUNIT_ASSERT(!(tCell1!=tCell2));
		}
	}

	// with shift and permutation
	{	
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;
		
		const auto sh = rand<fMat>(D,1,.2,.8);
		tCell1.shift(sh);
		
		const auto n1 = genRndST(0,D-1);
		const auto n2 = genRndST(0,D-1);

		tCell1.swapDim(n1,n2);
		
		CPPUNIT_ASSERT(!(tCell1==tCell2));
		CPPUNIT_ASSERT(tCell1!=tCell2);
	}
}


const char* test_cell_comparison::test_id() noexcept {
	return "test_cell_comparison";
}

CppUnit::Test* test_cell_comparison::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_cell_comparison>(
		"test_getPmat", &test_cell_comparison::test_getPmat));
	suite->addTest(new CppUnit::TestCaller<test_cell_comparison>(
		"test_getAvec", &test_cell_comparison::test_getAvec));
	suite->addTest(new CppUnit::TestCaller<test_cell_comparison>(
		"test_sameLattice", &test_cell_comparison::test_sameLattice));
	suite->addTest(new CppUnit::TestCaller<test_cell_comparison>(
		"test_operator_equal_unequal", &test_cell_comparison::test_operator_equal_unequal));
	
	return suite;
}
