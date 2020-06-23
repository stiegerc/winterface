// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_el_arithmetic.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


template<class TT, class FT, class CT>
void test_tMat_el_arithmetic<TT,FT,CT>::test_operator_minus() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	const auto tMat1 = rnd<tMat>(M,N);
	const auto tMat2 = -tMat1;
	for (size_t i=0; i!=tMat1.L(); ++i)
		CPPUNIT_ASSERT_EQUAL(tMat1[i],-tMat2[i]);
}

template<class TT, class FT, class CT>
void test_tMat_el_arithmetic<TT,FT,CT>::test_operator_plus_pluseq() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// +
	{	
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ftMat = rnd<fMat>(M,N);
		const auto ctMat = rnd<cMat>(M,N);
	
		const auto farg = ftMat[0];
		const auto carg = ctMat[0];
	
		const auto tMat2 = tMat1+farg;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat2.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]+farg,tMat2[i]);
			
		const auto tMat3 = tMat1+carg;
		CPPUNIT_ASSERT(tMat3.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]+carg,tMat3[i]);
			
		const auto tMat4 = tMat1+ftMat;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat4.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]+ftMat[i],tMat4[i]);
			
		const auto tMat5 = tMat1+ctMat;
		CPPUNIT_ASSERT(tMat5.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]+ctMat[i],tMat5[i]);
	}

	// +=
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ftMat = rnd<fMat>(M,N);
		const auto ctMat = rnd<cMat>(M,N);
	
		const auto farg = ftMat[0];
		const auto carg = ctMat[0];
	
		auto tMat2 = tMat1; tMat2+=farg;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat2.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]+farg,tMat2[i]);
	
		auto tMat3 = tMat1; tMat3+=carg;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat3.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(cc<TT>(tMat1[i]+carg),tMat3[i]);
			
		auto tMat4 = tMat1; tMat4+=ftMat;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat4.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]+ftMat[i],tMat4[i]);
			
		auto tMat5 = tMat1; tMat5+=ctMat;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat5.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(cc<TT>(tMat1[i]+ctMat[i]),tMat5[i]);
	}
}

template<class TT, class FT, class CT>
void test_tMat_el_arithmetic<TT,FT,CT>::test_operator_minus_minuseq() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	// -
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ftMat = rnd<fMat>(M,N);
		const auto ctMat = rnd<cMat>(M,N);
	
		const auto farg = ftMat[0];
		const auto carg = ctMat[0];
	
		const auto tMat2 = tMat1-farg;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat2.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]-farg,tMat2[i]);
			
		const auto tMat3 = tMat1-carg;
		CPPUNIT_ASSERT(tMat3.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]-carg,tMat3[i]);
			
		const auto tMat4 = tMat1-ftMat;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat4.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]-ftMat[i],tMat4[i]);
			
		const auto tMat5 = tMat1-ctMat;
		CPPUNIT_ASSERT(tMat5.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]-ctMat[i],tMat5[i]);
	}

	// -=
	{	
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ftMat = rnd<fMat>(M,N);
		const auto ctMat = rnd<cMat>(M,N);
	
		const auto farg = ftMat[0];
		const auto carg = ctMat[0];
	
		auto tMat2 = tMat1; tMat2-=farg;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat2.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]-farg,tMat2[i]);
			
		auto tMat3 = tMat1; tMat3-=carg;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat3.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(cc<TT>(tMat1[i]-carg),tMat3[i]);
			
		auto tMat4 = tMat1; tMat4-=ftMat;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat4.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]-ftMat[i],tMat4[i]);
			
		auto tMat5 = tMat1; tMat5-=ctMat;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat5.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(cc<TT>(tMat1[i]-ctMat[i]),tMat5[i]);
	}
}

template<class TT, class FT, class CT>
void test_tMat_el_arithmetic<TT,FT,CT>::test_operator_prod_prodeq() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	// *
	{	
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ftMat = rnd<fMat>(M,N);
		const auto ctMat = rnd<cMat>(M,N);
	
		const auto farg = ftMat[0];
		const auto carg = ctMat[0];
	
		const auto tMat2 = tMat1*farg;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat2.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]*farg,tMat2[i]);
			
		const auto tMat3 = tMat1*carg;
		CPPUNIT_ASSERT(tMat3.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]*carg,tMat3[i]);
			
		const auto tMat4 = tMat1*ftMat;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat4.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]*ftMat[i],tMat4[i]);
			
		const auto tMat5 = tMat1*ctMat;
		CPPUNIT_ASSERT(tMat5.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]*ctMat[i],tMat5[i]);
	}

	// *=
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ftMat = rnd<fMat>(M,N);
		const auto ctMat = rnd<cMat>(M,N);
	
		const auto farg = ftMat[0];
		const auto carg = ctMat[0];
	
		auto tMat2 = tMat1; tMat2*=farg;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat2.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]*farg,tMat2[i]);
			
		auto tMat3 = tMat1; tMat3*=carg;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat3.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(cc<TT>(tMat1[i]*carg),tMat3[i]);
			
		auto tMat4 = tMat1; tMat4*=ftMat;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat4.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]*ftMat[i],tMat4[i]);
			
		auto tMat5 = tMat1; tMat5*=ctMat;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat5.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(cc<TT>(tMat1[i]*ctMat[i]),tMat5[i]);
	}
}

template<class TT, class FT, class CT>
void test_tMat_el_arithmetic<TT,FT,CT>::test_operator_div_diveq() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	// /
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ftMat = rnd<fMat>(M,N);
		const auto ctMat = rnd<cMat>(M,N);
	
		const auto farg = ftMat[0];
		const auto carg = ctMat[0];
	
		const auto tMat2 = tMat1/farg;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat2.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]/farg,tMat2[i]);
			
		const auto tMat3 = tMat1/carg;
		CPPUNIT_ASSERT(tMat3.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]/carg,tMat3[i]);
			
		const auto tMat4 = tMat1/ftMat;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat4.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]/ftMat[i],tMat4[i]);
			
		const auto tMat5 = tMat1/ctMat;
		CPPUNIT_ASSERT(tMat5.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]/ctMat[i],tMat5[i]);
	}

	// /=
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto ftMat = rnd<fMat>(M,N);
		const auto ctMat = rnd<cMat>(M,N);
	
		const auto farg = ftMat[0];
		const auto carg = ctMat[0];
	
		auto tMat2 = tMat1; tMat2/=farg;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat2.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]/farg,tMat2[i]);
			
		auto tMat3 = tMat1; tMat3/=carg;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat3.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_DELTA(cc<TT>(tMat1[i]/carg),tMat3[i],delta);
			
		auto tMat4 = tMat1; tMat4/=ftMat;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat4.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(tMat1[i]/ftMat[i],tMat4[i]);
			
		auto tMat5 = tMat1; tMat5/=ctMat;
		CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat5.cpx());
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_DELTA(cc<TT>(tMat1[i]/ctMat[i]),tMat5[i],delta);
	}
}

template<class TT, class FT, class CT>
void test_tMat_el_arithmetic<TT,FT,CT>::test_operator_plus_tlhs() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	const auto tMat1 = rnd<tMat>(M,N);
	const auto ftMat = rnd<fMat>(M,N);
	const auto ctMat = rnd<cMat>(M,N);
	
	const auto farg = ftMat[0];
	const auto carg = ctMat[0];
	
	const auto tMat2 = farg+tMat1;
	CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat2.cpx());
	for (size_t i=0; i!=tMat1.L(); ++i)
		CPPUNIT_ASSERT_EQUAL(tMat1[i]+farg,tMat2[i]);
			
	const auto tMat3 = carg+tMat1;
	CPPUNIT_ASSERT(tMat3.cpx());
	for (size_t i=0; i!=tMat1.L(); ++i)
		CPPUNIT_ASSERT_EQUAL(tMat1[i]+carg,tMat3[i]);
}

template<class TT, class FT, class CT>
void test_tMat_el_arithmetic<TT,FT,CT>::test_operator_minus_tlhs() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	const auto tMat1 = rnd<tMat>(M,N);
	const auto ftMat = rnd<fMat>(M,N);
	const auto ctMat = rnd<cMat>(M,N);
	
	const auto farg = ftMat[0];
	const auto carg = ctMat[0];
	
	const auto tMat2 = farg-tMat1;
	CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat2.cpx());
	for (size_t i=0; i!=tMat1.L(); ++i)
		CPPUNIT_ASSERT_EQUAL(farg-tMat1[i],tMat2[i]);
			
	const auto tMat3 = carg-tMat1;
	CPPUNIT_ASSERT(tMat3.cpx());
	for (size_t i=0; i!=tMat1.L(); ++i)
		CPPUNIT_ASSERT_EQUAL(carg-tMat1[i],tMat3[i]);
}

template<class TT, class FT, class CT>
void test_tMat_el_arithmetic<TT,FT,CT>::test_operator_prod_tlhs() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	const auto tMat1 = rnd<tMat>(M,N);
	const auto ftMat = rnd<fMat>(M,N);
	const auto ctMat = rnd<cMat>(M,N);
	
	const auto farg = ftMat[0];
	const auto carg = ctMat[0];
	
	const auto tMat2 = farg*tMat1;
	CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat2.cpx());
	for (size_t i=0; i!=tMat1.L(); ++i)
		CPPUNIT_ASSERT_EQUAL(farg*tMat1[i],tMat2[i]);
			
	const auto tMat3 = carg*tMat1;
	CPPUNIT_ASSERT(tMat3.cpx());
	for (size_t i=0; i!=tMat1.L(); ++i)
		CPPUNIT_ASSERT_EQUAL(carg*tMat1[i],tMat3[i]);
}

template<class TT, class FT, class CT>
void test_tMat_el_arithmetic<TT,FT,CT>::test_operator_div_tlhs() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	const auto tMat1 = rnd<tMat>(M,N);
	const auto ftMat = rnd<fMat>(M,N);
	const auto ctMat = rnd<cMat>(M,N);
	
	const auto farg = ftMat[0];
	const auto carg = ctMat[0];
	
	const auto tMat2 = farg/tMat1;
	CPPUNIT_ASSERT_EQUAL(tMat1.cpx(),tMat2.cpx());
	for (size_t i=0; i!=tMat1.L(); ++i)
		CPPUNIT_ASSERT_EQUAL(farg/tMat1[i],tMat2[i]);
			
	const auto tMat3 = carg/tMat1;
	CPPUNIT_ASSERT(tMat3.cpx());
	for (size_t i=0; i!=tMat1.L(); ++i)
		CPPUNIT_ASSERT_EQUAL(carg/tMat1[i],tMat3[i]);
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tMat_el_arithmetic<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tMat_el_arithmetic>(
		"test_operator_minus", &test_tMat_el_arithmetic<TT,FT,CT>::test_operator_minus));
	suite->addTest(new CppUnit::TestCaller<test_tMat_el_arithmetic>(
		"test_operator_plus_pluseq", &test_tMat_el_arithmetic<TT,FT,CT>::test_operator_plus_pluseq));
	suite->addTest(new CppUnit::TestCaller<test_tMat_el_arithmetic>(
		"test_operator_minus_minuseq", &test_tMat_el_arithmetic<TT,FT,CT>::test_operator_minus_minuseq));
	suite->addTest(new CppUnit::TestCaller<test_tMat_el_arithmetic>(
		"test_operator_prod_prodeq", &test_tMat_el_arithmetic<TT,FT,CT>::test_operator_prod_prodeq));
	suite->addTest(new CppUnit::TestCaller<test_tMat_el_arithmetic>(
		"test_operator_div_diveq", &test_tMat_el_arithmetic<TT,FT,CT>::test_operator_div_diveq));
	suite->addTest(new CppUnit::TestCaller<test_tMat_el_arithmetic>(
		"test_operator_plus_tlhs", &test_tMat_el_arithmetic<TT,FT,CT>::test_operator_plus_tlhs));
	suite->addTest(new CppUnit::TestCaller<test_tMat_el_arithmetic>(
		"test_operator_minus_tlhs", &test_tMat_el_arithmetic<TT,FT,CT>::test_operator_minus_tlhs));
	suite->addTest(new CppUnit::TestCaller<test_tMat_el_arithmetic>(
		"test_operator_prod_tlhs", &test_tMat_el_arithmetic<TT,FT,CT>::test_operator_prod_tlhs));
	suite->addTest(new CppUnit::TestCaller<test_tMat_el_arithmetic>(
		"test_operator_div_tlhs", &test_tMat_el_arithmetic<TT,FT,CT>::test_operator_div_tlhs));
	
	return suite;
}
