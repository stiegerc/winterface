// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_basic_properties_info.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

#include "blas.h"

using namespace lm__;
using namespace lm__::test;

template<class TT, class FT, class CT>
void test_tMat_basic_properties_info<TT,FT,CT>::test_row_col_cpx() {

	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,N);
		CPPUNIT_ASSERT(!tMat1.row());
		CPPUNIT_ASSERT(!tMat1.col());
	
		const auto m = genRndST(0,M-1);
		auto jr = tMat1.crBegin()+m;
		CPPUNIT_ASSERT(jr->row());
		CPPUNIT_ASSERT(!jr->col());
		
		const auto n = genRndST(0,N-1);
		auto jc = tMat1.ccBegin()+n;
		CPPUNIT_ASSERT(!jc->row());
		CPPUNIT_ASSERT(jc->col());
	}
	{
		const auto tMat1 = rnd<tMat>(1,N);
		CPPUNIT_ASSERT(tMat1.row());
		CPPUNIT_ASSERT(!tMat1.col());
	}
	{
		const auto tMat1 = rnd<tMat>(M,1);
		CPPUNIT_ASSERT(!tMat1.row());
		CPPUNIT_ASSERT(tMat1.col());
	}

	{
		const auto fMat1 = rnd<fMat>(M,N);
		CPPUNIT_ASSERT(!fMat1.cpx());

		const auto cMat1 = rnd<cMat>(M,N);
		CPPUNIT_ASSERT(cMat1.cpx());
	}
}

template<class TT, class FT, class CT>
void test_tMat_basic_properties_info<TT,FT,CT>::test_M_N_L_size_lcap_ccap_empty_incr() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	const auto tMat1 = rnd<tMat>(M,N);
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
	CPPUNIT_ASSERT_EQUAL(M*N,tMat1.L());
	CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat1.size());
	CPPUNIT_ASSERT_EQUAL(M*N,tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.ccap());
	CPPUNIT_ASSERT(!tMat1.empty());
	CPPUNIT_ASSERT_EQUAL(size_t(1),tMat1.incr());
}

template<class TT, class FT, class CT>
void test_tMat_basic_properties_info<TT,FT,CT>::test_square_hermitian_diag() {

	// square
	{
		const auto S = genRndMEST();

		const auto tMat1 = rnd<tMat>(S.M,S.M);
		CPPUNIT_ASSERT(tMat1.square());
		
		const auto tMat2 = rnd<tMat>(S.M,S.N);
		CPPUNIT_ASSERT(!tMat2.square());
	}

	// hermitian
	{
		const size_t M = genRndST(3,3);

		auto tMat1 = rnd<tMat>(M,M);
		tMat1(2,1) = tMat1(1,2)*TT(2.0);

		CPPUNIT_ASSERT(!tMat1.hermitian());

		for (size_t m=0; m!=M; ++m)
			tMat1(m,m) = std::real(tMat1(m,m));
		for (size_t n=1; n!=M; ++n)
			for (size_t m=0; m!=n; ++m)
				lm__::ops::assign(tMat1(m,n),std::conj(tMat1(n,m)));

		CPPUNIT_ASSERT(tMat1.hermitian());
	}
	
	// diag
	{
		const size_t M = genRndST();

		auto tMat1 = rnd<tMat>(M,M)+TT(1.0);
		tMat1(2,1) = tMat1(2,1)+TT(10.0);
		tMat1(1,1) = tMat1(1,1)+TT(10.0);

		CPPUNIT_ASSERT(!tMat1.diag());

		for (size_t n=0; n!=tMat1.N(); ++n)
			for(size_t m=0; m!=tMat1.M(); ++m)
				if (m!=n) tMat1(m,n) = TT(0.0);
		CPPUNIT_ASSERT(tMat1.diag());
	}
}

template<class TT, class FT, class CT>
void test_tMat_basic_properties_info<TT,FT,CT>::test_ob_onb() {
	
	const size_t M = genRndST();
//	const size_t M = 3;
	CPPUNIT_ASSERT(!rnd<tMat>(1,M).ob());
	CPPUNIT_ASSERT(!rnd<tMat>(1,M).onb());

#ifndef NTOLERANT__
	auto tMat1 = rnd_b<tMat>(M,M);
//	std::cout << "\n\n" << tMat1 << "\n\n";
//	std::cout << dot(tMat1.cAt(0),tMat1.cAt(1)) << "\n\n";
//	const bool affe = tMat1.ob();
//	std::cout << "\n\n" << affe << "\n\n";
//	assert(false);

/*	
	{
		cMat tMat1 = rnd_b<cMat>(3,3);
		std::cout << "\n\n" << tMat1 << "\n\n";
		std::cout << "BLAS: ";
		CT res = c_xdotc(3,tMat1.data(),1,tMat1.data()+3,1);
		std::cout << "res: " << res << "\n";


		std::complex<double> rr;
		const int n = 3; const int inc = 1;
		FORTRAN_NAME(zdotc,ZDOTC)(&rr,&n,tMat1.data(),&inc,tMat1.data()+3,&inc);
		std::cout << "rr: " << rr << "\n";

		assert(false);
	}
*/
	CPPUNIT_ASSERT(!tMat1.ob());
	CPPUNIT_ASSERT(!tMat1.onb());

	gsorth(tMat1);
	CPPUNIT_ASSERT(tMat1.ob());
	CPPUNIT_ASSERT(tMat1.onb());

	const auto rnd = rand<tMat>(1,M,FT(2.0),FT(3.0));
	for (size_t n=0; n!=M; ++n)
		for (size_t m=0; m!=M; ++m)
			tMat1(m,n)*=rnd[n];

	CPPUNIT_ASSERT(tMat1.ob());
	CPPUNIT_ASSERT(!tMat1.onb());
#else
	tMat tMat1 = eye(M,M);
	CPPUNIT_ASSERT(tMat1.ob());
	CPPUNIT_ASSERT(tMat1.onb());
	
	const auto r = rnd<tMat>(1,M,1.0,2.0);
	for (size_t n=0; n!=M; ++n)
		for (size_t m=0; m!=M; ++m)
			if (m==n) tMat1(m,n) = r[m];
	CPPUNIT_ASSERT(tMat1.ob());
	CPPUNIT_ASSERT(!tMat1.onb());
#endif
}

template<class TT, class FT, class CT>
CppUnit::Test* test_tMat_basic_properties_info<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tMat_basic_properties_info>(
		"test_row_col_cpx", &test_tMat_basic_properties_info<TT,FT,CT>::test_row_col_cpx));
	suite->addTest(new CppUnit::TestCaller<test_tMat_basic_properties_info>(
		"test_M_N_L_size_lcap_ccap_empty_incr", &test_tMat_basic_properties_info<TT,FT,CT>::test_M_N_L_size_lcap_ccap_empty_incr));
	suite->addTest(new CppUnit::TestCaller<test_tMat_basic_properties_info>(
		"test_square_hermitian_diag", &test_tMat_basic_properties_info<TT,FT,CT>::test_square_hermitian_diag));
	suite->addTest(new CppUnit::TestCaller<test_tMat_basic_properties_info>(
		"test_ob_onb", &test_tMat_basic_properties_info<TT,FT,CT>::test_ob_onb));
	
	return suite;
}
