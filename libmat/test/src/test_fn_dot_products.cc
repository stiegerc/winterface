// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_fn_dot_products.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


// tests
void test_fn_dot_products::test_dot() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		RE__ ref(0.0);
		for (size_t i=0; i!=tMat1.L(); ++i) ref+=tMat1[i]*tMat2[i];
		const auto act = dot(tMat1,tMat2);

		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}
	{
		const auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		CPX__ ref(0.0);
		for (size_t i=0; i!=tMat1.L(); ++i) ref+=tMat1[i]*tMat2[i];
		const auto act = dot(tMat1,tMat2);

		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}
	{
		const auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		CPX__ ref(0.0);
		for (size_t i=0; i!=tMat1.L(); ++i) ref+=std::conj(tMat1[i])*tMat2[i];
		const auto act = dot(tMat1,tMat2);

		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}
	{
		const auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		CPX__ ref(0.0);
		for (size_t i=0; i!=tMat1.L(); ++i) ref+=std::conj(tMat1[i])*tMat2[i];
		const auto act = dot(tMat1,tMat2);

		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}

	{
		const auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		auto r1 = tMat1.rBegin();
		for (size_t m1=0; m1!=M; ++m1, ++r1) {

			auto r2 = tMat2.rBegin();
			for (size_t m2=0; m2!=M; ++m2, ++r2) {

				RE__ ref(0.0);
				for (size_t n=0; n!=N; ++n) ref+=tMat1(m1,n)*tMat2(m2,n);
				const auto act = dot(*r1,*r2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
	{
		const auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		auto r1 = tMat1.rBegin();
		for (size_t m1=0; m1!=M; ++m1, ++r1) {

			auto r2 = tMat2.rBegin();
			for (size_t m2=0; m2!=M; ++m2, ++r2) {

				CPX__ ref(0.0);
				for (size_t n=0; n!=N; ++n) ref+=tMat1(m1,n)*tMat2(m2,n);
				const auto act = dot(*r1,*r2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
	{
		const auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		auto r1 = tMat1.rBegin();
		for (size_t m1=0; m1!=M; ++m1, ++r1) {

			auto r2 = tMat2.rBegin();
			for (size_t m2=0; m2!=M; ++m2, ++r2) {

				CPX__ ref(0.0);
				for (size_t n=0; n!=N; ++n) ref+=std::conj(tMat1(m1,n))*tMat2(m2,n);
				const auto act = dot(*r1,*r2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
	{
		const auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		auto r1 = tMat1.rBegin();
		for (size_t m1=0; m1!=M; ++m1, ++r1) {

			auto r2 = tMat2.rBegin();
			for (size_t m2=0; m2!=M; ++m2, ++r2) {

				CPX__ ref(0.0);
				for (size_t n=0; n!=N; ++n) ref+=std::conj(tMat1(m1,n))*tMat2(m2,n);
				const auto act = dot(*r1,*r2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}

	{
		const auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		auto c1 = tMat1.cBegin();
		for (size_t n1=0; n1!=N; ++n1, ++c1) {

			auto c2 = tMat2.cBegin();
			for (size_t n2=0; n2!=N; ++n2, ++c2) {

				RE__ ref(0.0);
				for (size_t m=0; m!=M; ++m) ref+=tMat1(m,n1)*tMat2(m,n2);
				const auto act = dot(*c1,*c2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
	{
		const auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		auto c1 = tMat1.cBegin();
		for (size_t n1=0; n1!=N; ++n1, ++c1) {

			auto c2 = tMat2.cBegin();
			for (size_t n2=0; n2!=N; ++n2, ++c2) {

				CPX__ ref(0.0);
				for (size_t m=0; m!=M; ++m) ref+=tMat1(m,n1)*tMat2(m,n2);
				const auto act = dot(*c1,*c2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
	{
		const auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		auto c1 = tMat1.cBegin();
		for (size_t n1=0; n1!=N; ++n1, ++c1) {

			auto c2 = tMat2.cBegin();
			for (size_t n2=0; n2!=N; ++n2, ++c2) {

				CPX__ ref(0.0);
				for (size_t m=0; m!=M; ++m) ref+=std::conj(tMat1(m,n1))*tMat2(m,n2);
				const auto act = dot(*c1,*c2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
	{
		const auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		auto c1 = tMat1.cBegin();
		for (size_t n1=0; n1!=N; ++n1, ++c1) {

			auto c2 = tMat2.cBegin();
			for (size_t n2=0; n2!=N; ++n2, ++c2) {

				CPX__ ref(0.0);
				for (size_t m=0; m!=M; ++m) ref+=std::conj(tMat1(m,n1))*tMat2(m,n2);
				const auto act = dot(*c1,*c2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
}

void test_fn_dot_products::test_dotu() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		RE__ ref(0.0);
		for (size_t i=0; i!=tMat1.L(); ++i) ref+=tMat1[i]*tMat2[i];
		const auto act = dotu(tMat1,tMat2);

		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}
	{
		const auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		CPX__ ref(0.0);
		for (size_t i=0; i!=tMat1.L(); ++i) ref+=tMat1[i]*tMat2[i];
		const auto act = dotu(tMat1,tMat2);

		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}
	{
		const auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		CPX__ ref(0.0);
		for (size_t i=0; i!=tMat1.L(); ++i) ref+=tMat1[i]*tMat2[i];
		const auto act = dotu(tMat1,tMat2);

		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}
	{
		const auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		CPX__ ref(0.0);
		for (size_t i=0; i!=tMat1.L(); ++i) ref+=tMat1[i]*tMat2[i];
		const auto act = dotu(tMat1,tMat2);

		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}

	{
		const auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		auto r1 = tMat1.rBegin();
		for (size_t m1=0; m1!=M; ++m1, ++r1) {

			auto r2 = tMat2.rBegin();
			for (size_t m2=0; m2!=M; ++m2, ++r2) {

				RE__ ref(0.0);
				for (size_t n=0; n!=N; ++n) ref+=tMat1(m1,n)*tMat2(m2,n);
				const auto act = dotu(*r1,*r2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
	{
		const auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		auto r1 = tMat1.rBegin();
		for (size_t m1=0; m1!=M; ++m1, ++r1) {

			auto r2 = tMat2.rBegin();
			for (size_t m2=0; m2!=M; ++m2, ++r2) {

				CPX__ ref(0.0);
				for (size_t n=0; n!=N; ++n) ref+=tMat1(m1,n)*tMat2(m2,n);
				const auto act = dotu(*r1,*r2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
	{
		const auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		auto r1 = tMat1.rBegin();
		for (size_t m1=0; m1!=M; ++m1, ++r1) {

			auto r2 = tMat2.rBegin();
			for (size_t m2=0; m2!=M; ++m2, ++r2) {

				CPX__ ref(0.0);
				for (size_t n=0; n!=N; ++n) ref+=tMat1(m1,n)*tMat2(m2,n);
				const auto act = dotu(*r1,*r2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
	{
		const auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		auto r1 = tMat1.rBegin();
		for (size_t m1=0; m1!=M; ++m1, ++r1) {

			auto r2 = tMat2.rBegin();
			for (size_t m2=0; m2!=M; ++m2, ++r2) {

				CPX__ ref(0.0);
				for (size_t n=0; n!=N; ++n) ref+=tMat1(m1,n)*tMat2(m2,n);
				const auto act = dotu(*r1,*r2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}

	{
		const auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		auto c1 = tMat1.cBegin();
		for (size_t n1=0; n1!=N; ++n1, ++c1) {

			auto c2 = tMat2.cBegin();
			for (size_t n2=0; n2!=N; ++n2, ++c2) {

				RE__ ref(0.0);
				for (size_t m=0; m!=M; ++m) ref+=tMat1(m,n1)*tMat2(m,n2);
				const auto act = dotu(*c1,*c2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
	{
		const auto tMat1 = rnd<fMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		auto c1 = tMat1.cBegin();
		for (size_t n1=0; n1!=N; ++n1, ++c1) {

			auto c2 = tMat2.cBegin();
			for (size_t n2=0; n2!=N; ++n2, ++c2) {

				CPX__ ref(0.0);
				for (size_t m=0; m!=M; ++m) ref+=tMat1(m,n1)*tMat2(m,n2);
				const auto act = dotu(*c1,*c2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
	{
		const auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<fMat>(M,N);

		auto c1 = tMat1.cBegin();
		for (size_t n1=0; n1!=N; ++n1, ++c1) {

			auto c2 = tMat2.cBegin();
			for (size_t n2=0; n2!=N; ++n2, ++c2) {

				CPX__ ref(0.0);
				for (size_t m=0; m!=M; ++m) ref+=tMat1(m,n1)*tMat2(m,n2);
				const auto act = dotu(*c1,*c2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
	{
		const auto tMat1 = rnd<cMat>(M,N);
		const auto tMat2 = rnd<cMat>(M,N);

		auto c1 = tMat1.cBegin();
		for (size_t n1=0; n1!=N; ++n1, ++c1) {

			auto c2 = tMat2.cBegin();
			for (size_t n2=0; n2!=N; ++n2, ++c2) {

				CPX__ ref(0.0);
				for (size_t m=0; m!=M; ++m) ref+=tMat1(m,n1)*tMat2(m,n2);
				const auto act = dotu(*c1,*c2);

				CPPUNIT_ASSERT_DELTA(ref,act,delta);
			}
		}
	}
}


// test id
const char* test_fn_dot_products::test_id() noexcept {
	return "test_fn_dot_products";
}

CppUnit::Test* test_fn_dot_products::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_fn_dot_products>(
		"test_dot", &test_fn_dot_products::test_dot));
	suite->addTest(new CppUnit::TestCaller<test_fn_dot_products>(
		"test_dotu", &test_fn_dot_products::test_dotu));
	
	return suite;
}
