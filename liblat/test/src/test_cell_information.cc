// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_cell_information.h"
#include "ll_testTools.h"
#include "aux_io.h"
#include <iostream>

using namespace ll__;
using namespace ll__::test;
using namespace lm__;
using namespace lm__::test;
using namespace aux;

void test_cell_information::test_empty() {
	{
		const auto tCell = cell();
		CPPUNIT_ASSERT(tCell.empty());
	}
	{
		auto tCell = genRandom();
		const auto B = tCell.moveB();
		CPPUNIT_ASSERT(!tCell.empty());
	}
	{
		auto tCell = genRandom();
		const auto Ap = tCell.moveAp();
		CPPUNIT_ASSERT(tCell.empty());
	}
	{
		auto tCell = genRandom();
		const auto B = tCell.moveB();
		const auto Ap = tCell.moveAp();
		CPPUNIT_ASSERT(tCell.empty());
	}
}

void test_cell_information::test_dim() {
	const auto D = genRndST(1,5);
	const auto tCell = genRandom(D);
	CPPUNIT_ASSERT_EQUAL(D,tCell.dim());
}

void test_cell_information::test_N_Nspecies() {
	const auto D = genRndST(1,5);
	const auto tCell1 = genRandom(D);

	CPPUNIT_ASSERT_EQUAL(tCell1.Ap().N(),tCell1.N());
	CPPUNIT_ASSERT_EQUAL(size_t(tCell1.T_.size()-1),tCell1.Nspecies());

	auto tCell2 = tCell1;
}

void test_cell_information::test_vol_sign() {
	const auto tCell = genRandom();
	const auto det_ = det(tCell.B());
	const double sign = std::signbit(det_) ? -1.0: 1.0;

	CPPUNIT_ASSERT_DELTA(std::abs(det_),tCell.vol(),mtol());
	CPPUNIT_ASSERT_EQUAL(sign,tCell.sign());
}

void test_cell_information::test_find() {
	const auto D = genRndST(1,5);
	auto tCell = genRandom(D);

	for (const auto t: tCell.types())
		for (auto i=tCell.ccBegin(t),e=tCell.ccEnd(t); i!=e; ++i) {
			CPPUNIT_ASSERT(tCell.find(
					*i
					+randi<fMat>(tCell.dim(),1,-100,100)
					-genRndDouble(0,mtol())
					,t)!=tCell.ccEnd(t));
		}
}

void test_cell_information::test_lVec() {
	const auto D = genRndST(1,5);
	auto tCell = genRandom(D);

	const auto C = rnd_i(D,-3,3);
	
	CPPUNIT_ASSERT(tCell.lVec(C));
	for (auto i=C.ccBegin(), e=C.ccEnd(); i!=e; ++i)
		CPPUNIT_ASSERT(tCell.lVec(*i));

	const auto NB = tCell.B().prod(C);
	const auto LV = NB.leftDivide(tCell.B());
	tCell.changeBasis(NB);

	CPPUNIT_ASSERT(tCell.lVec(LV));
	for (auto i=LV.ccBegin(),e=LV.ccEnd(); i!=e; ++i)
		CPPUNIT_ASSERT(tCell.lVec(*i));
}

void test_cell_information::test_primitive() {
	const auto D = genRndST(1,5);
	auto tCell = genRandom(D);
	
	CPPUNIT_ASSERT(tCell.primitive());

	const auto C = rnd_i(D,-3,3);
	const auto NB = tCell.B().prod(C);
	tCell.changeBasis(NB);

	CPPUNIT_ASSERT(!tCell.primitive());
}

void test_cell_information::test_validBasis() {
	const auto D = genRndST(1,5);
	auto tCell = genRandom(D);

	CPPUNIT_ASSERT(tCell.validBasis(tCell.B()));
	
	const auto OB = tCell.B();
	const auto C = rnd_i(D,-3,3);
	const auto NB = tCell.B().prod(C);

	tCell.changeBasis(NB);
	CPPUNIT_ASSERT(tCell.validBasis(OB));
}


const char* test_cell_information::test_id() noexcept {
	return "test_cell_information";
}

CppUnit::Test* test_cell_information::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_cell_information>(
		"test_empty", &test_cell_information::test_empty));
	suite->addTest(new CppUnit::TestCaller<test_cell_information>(
		"test_dim", &test_cell_information::test_dim));
	suite->addTest(new CppUnit::TestCaller<test_cell_information>(
		"test_N_Nspecies", &test_cell_information::test_N_Nspecies));
	suite->addTest(new CppUnit::TestCaller<test_cell_information>(
		"test_vol_sign", &test_cell_information::test_vol_sign));
	suite->addTest(new CppUnit::TestCaller<test_cell_information>(
		"test_find", &test_cell_information::test_find));
	suite->addTest(new CppUnit::TestCaller<test_cell_information>(
		"test_lVec", &test_cell_information::test_lVec));
	suite->addTest(new CppUnit::TestCaller<test_cell_information>(
		"test_primitive", &test_cell_information::test_primitive));
	suite->addTest(new CppUnit::TestCaller<test_cell_information>(
		"test_validBasis", &test_cell_information::test_validBasis));
	
	return suite;
}
