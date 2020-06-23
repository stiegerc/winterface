// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_cell_assign.h"
#include "ll_testTools.h"
#include "aux_io.h"
#include <iostream>

using namespace ll__;
using namespace ll__::test;
using namespace aux;

void test_cell_assign::test_swap() {
	auto tCell1 = genRandom();
	auto tCell2 = genRandom();

	const auto Bptr1 = tCell1.B_.data();
	const auto Apptr1 = tCell1.mat_.data();
	const auto Tptr1 = tCell1.T_.data();
	const auto idptr1 = tCell1.id_.data();

	const auto Bptr2 = tCell2.B_.data();
	const auto Apptr2 = tCell2.mat_.data();
	const auto Tptr2 = tCell2.T_.data();
	const auto idptr2 = tCell2.id_.data();

	const auto B1 = tCell1.B_;
	const auto Ap1 = tCell1.mat_;
	const auto T1 = tCell1.T_;
	const auto id1 = tCell1.id_;
	
	const auto B2 = tCell2.B_;
	const auto Ap2 = tCell2.mat_;
	const auto T2 = tCell2.T_;
	const auto id2 = tCell2.id_;

	swap(tCell1,tCell2);

	CPPUNIT_ASSERT_EQUAL(Bptr1,tCell2.B_.data());
	CPPUNIT_ASSERT_EQUAL(Apptr1,tCell2.mat_.data());
	CPPUNIT_ASSERT(Tptr1==tCell2.T_.data());
	CPPUNIT_ASSERT_EQUAL(Bptr2,tCell1.B_.data());
	CPPUNIT_ASSERT_EQUAL(Apptr2,tCell1.mat_.data());
	CPPUNIT_ASSERT(Tptr2==tCell1.T_.data());
	CPPUNIT_ASSERT_EQUAL(idptr1,tCell2.id_.data());
	CPPUNIT_ASSERT_EQUAL(idptr2,tCell1.id_.data());

	CPPUNIT_ASSERT(B1==tCell2.B_);
	CPPUNIT_ASSERT(Ap1==tCell2.mat_);
	CPPUNIT_ASSERT(T1==tCell2.T_);
	CPPUNIT_ASSERT(B2==tCell1.B_);
	CPPUNIT_ASSERT(Ap2==tCell1.mat_);
	CPPUNIT_ASSERT(T2==tCell1.T_);
	CPPUNIT_ASSERT(id1==tCell2.id_);
	CPPUNIT_ASSERT(id2==tCell1.id_);
}

void test_cell_assign::test_operator_equal() {
	
	// test copy assign
	{
		auto tCell1 = genRandom();
		const auto tCell2 = genRandom();

		tCell1 = tCell2;
		CPPUNIT_ASSERT(tCell1==tCell2);
	}

	// test move assign
	{
		auto tCell1 = genRandom();
		auto tCell2 = genRandom();
		const auto tCell3 = tCell2;

		double* B_ptr = tCell2.B_.data();
		double* Ap_ptr = tCell2.mat_.data();
		size_t* T_ptr = tCell2.T_.data();
		std::string* id_ptr = tCell2.id_.data();

		tCell1 = std::move(tCell2);
		CPPUNIT_ASSERT(tCell1==tCell3);
		CPPUNIT_ASSERT_EQUAL(B_ptr,tCell1.B_.data());
		CPPUNIT_ASSERT_EQUAL(Ap_ptr,tCell1.mat_.data());
		CPPUNIT_ASSERT_EQUAL(T_ptr,tCell1.T_.data());
		CPPUNIT_ASSERT_EQUAL(id_ptr,tCell1.id_.data());
	}
}


const char* test_cell_assign::test_id() noexcept {
	return "test_cell_assign";
}

CppUnit::Test* test_cell_assign::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_cell_assign>(
		"test_swap", &test_cell_assign::test_swap));
	suite->addTest(new CppUnit::TestCaller<test_cell_assign>(
		"test_operator_equal", &test_cell_assign::test_operator_equal));
	
	return suite;
}
