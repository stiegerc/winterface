// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_cell_type_information.h"
#include "ll_testTools.h"
#include "aux_io.h"
#include <algorithm>
#include <iostream>


using namespace ll__;
using namespace ll__::test;
using namespace lm__;
using namespace lm__::test;
using namespace aux;

void test_cell_type_information::test_types_inds_() {
	const auto D = genRndST();
	const auto tCell = genRndCell(D);

	// types
	CPPUNIT_ASSERT(rg(0,1,tCell.Nspecies())==tCell.types());

	// inds
	CPPUNIT_ASSERT(rg(0,1,tCell.N())==tCell.inds());
	
	// validType
	for (aT t=0; t!=tCell.Nspecies(); ++t)
		CPPUNIT_ASSERT(tCell.validType(t));
	CPPUNIT_ASSERT(!tCell.validType(tCell.Nspecies()));

	// validInd
	for (size_t i=0; i!=tCell.N(); ++i)
		CPPUNIT_ASSERT(tCell.validInd(i));
	CPPUNIT_ASSERT(!tCell.validType(tCell.N()));
}

void test_cell_type_information::test_type() {
	const auto D = genRndST();
	const auto tCell = genRndCell(D);

	for (size_t i=0; i!=tCell.N(); ++i) {
		const auto t = tCell.type(i);
		CPPUNIT_ASSERT(tCell.validType(t));
		CPPUNIT_ASSERT(i>=tCell.T_[t] && i<tCell.T_[t+1]);
	}

	const auto T = tCell.type(tCell.inds());
	CPPUNIT_ASSERT(std::all_of(T.cbegin(),T.cend(),[&tCell](const aT t){return tCell.validType(t);}));
	CPPUNIT_ASSERT(tCell.type(tCell.inds())==tCell.type());
}

void test_cell_type_information::test_id_() {
	const auto D = genRndST();
	const auto tCell = genRndCell(D);

	// id
	{
		idv id_; id_.reserve(tCell.Nspecies());
		for (const auto t: tCell.types())
			id_.push_back(tCell.id(t));
		CPPUNIT_ASSERT(id_==tCell.id());

		CPPUNIT_ASSERT(tCell.type()==tCell.type(tCell.id(tCell.type())));
	}

	// stripId
	{
		const std::string s1 = "AAA-123231";
		const std::string ref1 = "AAA";
		CPPUNIT_ASSERT_EQUAL(ref1,ll_cell::stripId(s1));
		
		const std::string s2 = "(A-1:B-2)_(B:AAA)_(H-4:J)";
		const std::string ref2 = "(A:B)_(B:AAA)_(H:J)";
		CPPUNIT_ASSERT_EQUAL(ref2,ll_cell::stripId(s2));
	}
}

void test_cell_type_information::test_Ntype() {
	const auto D = genRndST();
	const auto tCell = genRndCell(D);

	aCv N; N.reserve(tCell.Nspecies());
	for (aT t=0; t!=tCell.Nspecies(); ++t)
		N.push_back(tCell.Ntype(t));
	CPPUNIT_ASSERT_EQUAL(std::accumulate(N.cbegin(),N.cend(),aC(0)),tCell.N());
	CPPUNIT_ASSERT(N==tCell.Ntype(tCell.types()));
	CPPUNIT_ASSERT(N==tCell.Ntype());
}

void test_cell_type_information::test_leastFreqType_mostFreqType() {
	const auto D = genRndST();
	const auto tCell = genRndCell(D);

	const auto N = tCell.Ntype();
	CPPUNIT_ASSERT_EQUAL(*std::min_element(N.begin(),N.end()),tCell.Ntype(tCell.leastFreqType()));
	CPPUNIT_ASSERT_EQUAL(*std::max_element(N.begin(),N.end()),tCell.Ntype(tCell.mostFreqType()));
}

void test_cell_type_information::test_ind() {
	const auto D = genRndST();
	const auto tCell = genRndCell(D);
	
	for (auto t: tCell.types()) {
		const auto I = tCell.ind(t);
		CPPUNIT_ASSERT(std::all_of(I.cbegin(),I.cend(),[&tCell](const size_t i){return i<tCell.N();}));
		CPPUNIT_ASSERT(rg(tCell.T_[t],1,tCell.T_[t+1]-tCell.T_[t])==I);
	}
}


const char* test_cell_type_information::test_id() noexcept {
	return "test_cell_type_information";
}

CppUnit::Test* test_cell_type_information::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_cell_type_information>(
		"test_types_inds_", &test_cell_type_information::test_types_inds_));
	suite->addTest(new CppUnit::TestCaller<test_cell_type_information>(
		"test_type", &test_cell_type_information::test_type));
	suite->addTest(new CppUnit::TestCaller<test_cell_type_information>(
		"test_id_", &test_cell_type_information::test_id_));
	suite->addTest(new CppUnit::TestCaller<test_cell_type_information>(
		"test_Ntype", &test_cell_type_information::test_Ntype));
	suite->addTest(new CppUnit::TestCaller<test_cell_type_information>(
		"test_leastFreqType_mostFreqType", &test_cell_type_information::test_leastFreqType_mostFreqType));
	suite->addTest(new CppUnit::TestCaller<test_cell_type_information>(
		"test_ind", &test_cell_type_information::test_ind));
	
	return suite;
}
