// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_fn_vector_sums.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


// tests
void test_fn_vector_sums::test_sum_r() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);
		const RE__ f = rnd<fMat>(1,1)[0];

		auto r1 = tMat1.rBegin();
		auto r2 = tMat2.crBegin();
		for (size_t m=0; m!=M; ++m, ++r1, ++r2) {
			const auto ref = tMat1.rGet(m) + tMat2.rGet(m)*f;
			sum(*r1,*r2,f);
			
			for (size_t n=0; n!=N; ++n)
				CPPUNIT_ASSERT_DELTA(ref[n],tMat1(m,n),delta);
		}
	}
	{
		auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);
		const RE__ f = rnd<fMat>(1,1)[0];

		auto r1 = tMat1.rBegin();
		auto r2 = tMat2.crBegin();
		for (size_t m=0; m!=M; ++m, ++r1, ++r2) {
			const auto ref = tMat1.rGet(m) + real(tMat2.rGet(m))*f;
			sum(*r1,*r2,f);
			
			for (size_t n=0; n!=N; ++n)
				CPPUNIT_ASSERT_DELTA(ref[n],tMat1(m,n),delta);
		}
	}
	{
		auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);
		const CPX__ f = rnd<cMat>(1,1)[0];

		auto r1 = tMat1.rBegin();
		auto r2 = tMat2.crBegin();
		for (size_t m=0; m!=M; ++m, ++r1, ++r2) {
			const auto ref = tMat1.rGet(m) + real(tMat2.rGet(m)*f);
			sum(*r1,*r2,f);
			
			for (size_t n=0; n!=N; ++n)
				CPPUNIT_ASSERT_DELTA(ref[n],tMat1(m,n),delta);
		}
	}
	{
		auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);
		const RE__ f = rnd<fMat>(1,1)[0];

		auto r1 = tMat1.rBegin();
		auto r2 = tMat2.crBegin();
		for (size_t m=0; m!=M; ++m, ++r1, ++r2) {
			const auto ref = tMat1.rGet(m) + tMat2.rGet(m)*f;
			sum(*r1,*r2,f);
			
			for (size_t n=0; n!=N; ++n)
				CPPUNIT_ASSERT_DELTA(ref[n],tMat1(m,n),delta);
		}
	}
	{
		auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);
		const CPX__ f = rnd<cMat>(1,1)[0];

		auto r1 = tMat1.rBegin();
		auto r2 = tMat2.crBegin();
		for (size_t m=0; m!=M; ++m, ++r1, ++r2) {
			const auto ref = tMat1.rGet(m) + tMat2.rGet(m)*f;
			sum(*r1,*r2,f);
			
			for (size_t n=0; n!=N; ++n)
				CPPUNIT_ASSERT_DELTA(ref[n],tMat1(m,n),delta);
		}
	}
	{
		auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);
		const CPX__ f = rnd<cMat>(1,1)[0];

		auto r1 = tMat1.rBegin();
		auto r2 = tMat2.crBegin();
		for (size_t m=0; m!=M; ++m, ++r1, ++r2) {
			const auto ref = tMat1.rGet(m) + tMat2.rGet(m)*f;
			sum(*r1,*r2,f);
			
			for (size_t n=0; n!=N; ++n)
				CPPUNIT_ASSERT_DELTA(ref[n],tMat1(m,n),delta);
		}
	}
}

void test_fn_vector_sums::test_sum_c() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);
		const RE__ f = rnd<fMat>(1,1)[0];

		auto c1 = tMat1.cBegin();
		auto c2 = tMat2.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c1, ++c2) {
			const auto ref = tMat1.cGet(n) + tMat2.cGet(n)*f;
			sum(*c1,*c2,f);
			
			for (size_t m=0; m!=M; ++m)
				CPPUNIT_ASSERT_DELTA(ref[m],tMat1(m,n),delta);
		}
	}
	{
		auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);
		const RE__ f = rnd<fMat>(1,1)[0];

		auto c1 = tMat1.cBegin();
		auto c2 = tMat2.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c1, ++c2) {
			const auto ref = tMat1.cGet(n) + real(tMat2.cGet(n))*f;
			sum(*c1,*c2,f);
			
			for (size_t m=0; m!=M; ++m)
				CPPUNIT_ASSERT_DELTA(ref[m],tMat1(m,n),delta);
		}
	}
	{
		auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);
		const CPX__ f = rnd<cMat>(1,1)[0];

		auto c1 = tMat1.cBegin();
		auto c2 = tMat2.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c1, ++c2) {
			const auto ref = tMat1.cGet(n) + real(tMat2.cGet(n)*f);
			sum(*c1,*c2,f);
			
			for (size_t m=0; m!=M; ++m)
				CPPUNIT_ASSERT_DELTA(ref[m],tMat1(m,n),delta);
		}
	}
	{
		auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);
		const RE__ f = rnd<fMat>(1,1)[0];

		auto c1 = tMat1.cBegin();
		auto c2 = tMat2.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c1, ++c2) {
			const auto ref = tMat1.cGet(n) + tMat2.cGet(n)*f;
			sum(*c1,*c2,f);
			
			for (size_t m=0; m!=M; ++m)
				CPPUNIT_ASSERT_DELTA(ref[m],tMat1(m,n),delta);
		}
	}
	{
		auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);
		const CPX__ f = rnd<cMat>(1,1)[0];

		auto c1 = tMat1.cBegin();
		auto c2 = tMat2.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c1, ++c2) {
			const auto ref = tMat1.cGet(n) + tMat2.cGet(n)*f;
			sum(*c1,*c2,f);
			
			for (size_t m=0; m!=M; ++m)
				CPPUNIT_ASSERT_DELTA(ref[m],tMat1(m,n),delta);
		}
	}
	{
		auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);
		const CPX__ f = rnd<cMat>(1,1)[0];

		auto c1 = tMat1.cBegin();
		auto c2 = tMat2.ccBegin();
		for (size_t n=0; n!=N; ++n, ++c1, ++c2) {
			const auto ref = tMat1.cGet(n) + tMat2.cGet(n)*f;
			sum(*c1,*c2,f);
			
			for (size_t m=0; m!=M; ++m)
				CPPUNIT_ASSERT_DELTA(ref[m],tMat1(m,n),delta);
		}
	}
}


// test id
const char* test_fn_vector_sums::test_id() noexcept {
	return "test_fn_vector_sums";
}

CppUnit::Test* test_fn_vector_sums::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_fn_vector_sums>(
		"test_sum_r", &test_fn_vector_sums::test_sum_r));
	suite->addTest(new CppUnit::TestCaller<test_fn_vector_sums>(
		"test_sum_c", &test_fn_vector_sums::test_sum_c));
	
	return suite;
}
