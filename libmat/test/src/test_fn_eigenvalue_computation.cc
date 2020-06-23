// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_fn_eigenvalue_computation.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


// tests
void test_fn_eigenvalue_computation::test_eig() {
	// fMat
	{
		{
			const fMat tMat1({0.5,0.8660,-0.8660,0.5},2,2);
			const auto tMat2 = eig(tMat1);
			const cMat sol = {CPX__(0.5,0.8660),CPX__(0.5,-0.8660)};
	
			CPPUNIT_ASSERT_DELTA(tMat2[0],sol[0],delta);
			CPPUNIT_ASSERT_DELTA(tMat2[1],sol[1],delta);
		}
		{
			const size_t N = genRndST();
			const auto tMat1 = rnd_b<fMat>(N,N);

			const auto tr1 = trace(tMat1);
			const auto tr2 = sum(eig(tMat1));

			CPPUNIT_ASSERT_DELTA(std::imag(tr2),RE__(0.0),delta);
			CPPUNIT_ASSERT_DELTA(std::real(tr2),tr1,delta);
		}
	}

	// cMat
	{
		{
			const cMat tMat1({CPX__(1.0,4.0),CPX__(2.0,3.0),CPX__(3.0,2.0),CPX__(4.0,1.0)},2,2);
			auto tMat2 = eig(tMat1);
			std::sort(tMat2.begin(),tMat2.end(),[](const CPX__& a, const CPX__& b){return ops::lt(a,b);});
			const cMat sol = {CPX__(0.438447187,0.438447187),CPX__(4.56155281,4.56155281)};

			CPPUNIT_ASSERT_DELTA(tMat2[0],sol[0],delta);
			CPPUNIT_ASSERT_DELTA(tMat2[1],sol[1],delta);
		}
		{
			const size_t N = genRndST();
			const auto tMat1 = rnd_b<cMat>(N,N);

			const auto tr1 = trace(tMat1);
			const auto tr2 = sum(eig(tMat1));

			CPPUNIT_ASSERT_DELTA(tr2,tr1,delta);
		}
	}
}

void test_fn_eigenvalue_computation::test_eigr() {
	
	const size_t N = genRndST();
	
	// fMat	
	{
		const auto tMat1 = rnd_b<fMat>(N,N);
			
		const auto ck = eigr(tMat1);
		
		cMat Vr(N,0); Vr.reserve(N);
		cMat D = zeros<cMat>(N);
		for (size_t i=0; i!=tMat1.M(); ++i) {
			const auto tmp = eGet(ck.Er,ck.Ei,ck.Vr,i);
			D(i,i) = tmp.e;
			Vr.push_back(tmp.v);
		}
		
		const auto m1 = Vr.prod(D);
		const auto m2 = tMat1.prod(Vr);
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_DELTA(m1[i],m2[i],delta);
	}

	// cMat	
	{
		const auto tMat1 = rnd_b<cMat>(N,N);
			
		const auto ck = eigr(tMat1);
		
		const auto m1 = ck.Vr.prod(diag(ck.E));
		const auto m2 = tMat1.prod(ck.Vr);
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_DELTA(m1[i],m2[i],delta);
	}
}

void test_fn_eigenvalue_computation::test_eigl() {
	
	const size_t N = genRndST();
	
	// fMat	
	{
		const auto tMat1 = rnd_b<fMat>(N,N);
			
		const auto ck = eigl(tMat1);
		
		cMat Vl(N,0); Vl.reserve(N);
		cMat D = zeros<cMat>(N);
		for (size_t i=0; i!=tMat1.M(); ++i) {
			const auto tmp = eGet(ck.Er,ck.Ei,ck.Vl,i);
			D(i,i) = tmp.e;
			Vl.push_back(tmp.v);
		}
		
		const auto m1 = D.prod(T(Vl));
		const auto m2 = T(Vl).prod(tMat1);
		
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_DELTA(m1[i],m2[i],delta);
	}

	// cMat	
	{
		const auto tMat1 = rnd_b<cMat>(N,N);
			
		const auto ck = eigl(tMat1);
		
		const auto m1 = diag(ck.E).prod(T(ck.Vl));
		const auto m2 = T(ck.Vl).prod(tMat1);
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_DELTA(m1[i],m2[i],delta);
	}
}

void test_fn_eigenvalue_computation::test_eigrl() {
	
	const size_t N = genRndST();
	
	// fMat	
	{
		const auto tMat1 = rnd_b<fMat>(N,N);
			
		const auto ck = eigrl(tMat1);
		
		cMat Vr(N,0); Vr.reserve(N);
		cMat Vl(N,0); Vl.reserve(N);
		cMat D = zeros<cMat>(N);
		for (size_t i=0; i!=tMat1.M(); ++i) {
			const auto tmpr = eGet(ck.Er,ck.Ei,ck.Vr,i);
			D(i,i) = tmpr.e;
			Vr.push_back(tmpr.v);
			
			const auto tmpl = eGet(ck.Er,ck.Ei,ck.Vl,i);
			Vl.push_back(tmpl.v);
		}
		
		{
			const auto m1 = Vr.prod(D);
			const auto m2 = tMat1.prod(Vr);
			for (size_t i=0; i!=tMat1.L(); ++i)
				CPPUNIT_ASSERT_DELTA(m1[i],m2[i],delta);
		}
		{
			const auto m1 = D.prod(T(Vl));
			const auto m2 = T(Vl).prod(tMat1);
			for (size_t i=0; i!=tMat1.L(); ++i)
				CPPUNIT_ASSERT_DELTA(m1[i],m2[i],delta);
		}
	}

	// cMat	
	{
		const auto tMat1 = rnd_b<cMat>(N,N);
		
		const auto ck = eigrl(tMat1);
		
		{
			const auto m1 = ck.Vr.prod(diag(ck.E));
			const auto m2 = tMat1.prod(ck.Vr);
			for (size_t i=0; i!=tMat1.L(); ++i)
				CPPUNIT_ASSERT_DELTA(m1[i],m2[i],delta);
		}
		{
			const auto m1 = diag(ck.E).prod(T(ck.Vl));
			const auto m2 = T(ck.Vl).prod(tMat1);
			for (size_t i=0; i!=tMat1.L(); ++i)
				CPPUNIT_ASSERT_DELTA(m1[i],m2[i],delta);
		}
	}
}

void test_fn_eigenvalue_computation::test_eigh() {
	// fMat
	{
		{
			const fMat tMat1({RE__(1.0),RE__(2.0),RE__(2.0),RE__(1.0)},2,2);
			const auto tMat2 = eigh(tMat1);
			const fMat sol = {-1.0,3.0};
		
			CPPUNIT_ASSERT(!tMat2.cpx());
			CPPUNIT_ASSERT_DELTA(tMat2[0],sol[0],delta);
			CPPUNIT_ASSERT_DELTA(tMat2[1],sol[1],delta);
		}
		{
			const size_t N = genRndST();
			fMat tMat1;	
			do { tMat1 = rnd<fMat>(N,N); } while(std::abs(det(tMat1))>RE__(1.0));
			tMat1 = tMat1.prod(T(tMat1));

			const auto tr1 = trace(tMat1);
			const auto tr2 = sum(eigh(tMat1));

			CPPUNIT_ASSERT_DELTA(tr2,tr1,delta);
		}
	}	
	
	// cMat
	{
		{
			const cMat tMat1({CPX__(1.0,0.0),CPX__(2.0,-1.0),CPX__(2.0,1.0),CPX__(1.0,0.0)},2,2);
			const auto tMat2 = eigh(tMat1);
			const fMat sol = {-1.2360679775,3.2360679775};
		
			CPPUNIT_ASSERT(!tMat2.cpx());
			CPPUNIT_ASSERT_DELTA(tMat2[0],sol[0],delta);
			CPPUNIT_ASSERT_DELTA(tMat2[1],sol[1],delta);
		}
		{
			const size_t N = genRndST();
			cMat tMat1;	
			do { tMat1 = rnd<cMat>(N,N); } while(std::abs(det(tMat1))>RE__(1.0));
			tMat1 = tMat1.prod(T(tMat1));

			const auto tr1 = trace(tMat1);
			const auto tr2 = sum(eigh(tMat1));

			CPPUNIT_ASSERT_DELTA(tr2,tr1,delta);
		}
	}
}

void test_fn_eigenvalue_computation::test_eighv() {
	
	const size_t N = genRndST();
	
	// fMat	
	{
		auto tMat1 = rnd_b<fMat>(N,N);
		tMat1 = tMat1.prod(T(tMat1));
			
		const auto ck = eighv(tMat1);
		CPPUNIT_ASSERT(!ck.V.cpx());
		for (size_t i=0; i!=N; ++i)
			for (size_t j=0; j!=N; ++j) {
				const auto tmp = dot(ck.V.cAt(i),ck.V.cAt(j));

				if (i==j) {
					CPPUNIT_ASSERT_DELTA(tmp,RE__(1.0),delta);
				} else {
					CPPUNIT_ASSERT_DELTA(tmp,RE__(0.0),delta);
				}
			}

		const auto m1 = ck.V.prod(diag(ck.E));
		const auto m2 = tMat1.prod(ck.V);
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_DELTA(m1[i],m2[i],delta);
	}

	// cMat	
	{
		auto tMat1 = rnd_b<cMat>(N,N);
		tMat1 = tMat1.prod(T(tMat1));
			
		const auto ck = eighv(tMat1);
		for (size_t i=0; i!=N; ++i)
			for (size_t j=0; j!=N; ++j) {
				const auto tmp = dot(ck.V.cAt(i),ck.V.cAt(j));

				if (i==j) {
					CPPUNIT_ASSERT_DELTA(tmp,RE__(1.0),delta);
				} else {
					CPPUNIT_ASSERT_DELTA(tmp,RE__(0.0),delta);
				}
			}

		const auto m1 = ck.V.prod(diag(ck.E));
		const auto m2 = tMat1.prod(ck.V);
		
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_DELTA(m1[i],m2[i],delta);
	}
}


// test id
const char* test_fn_eigenvalue_computation::test_id() noexcept {
	return "test_fn_eigenvalue_computation";
}

CppUnit::Test* test_fn_eigenvalue_computation::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_fn_eigenvalue_computation>(
		"test_eig", &test_fn_eigenvalue_computation::test_eig));
	suite->addTest(new CppUnit::TestCaller<test_fn_eigenvalue_computation>(
		"test_eigr", &test_fn_eigenvalue_computation::test_eigr));
	suite->addTest(new CppUnit::TestCaller<test_fn_eigenvalue_computation>(
		"test_eigl", &test_fn_eigenvalue_computation::test_eigl));
	suite->addTest(new CppUnit::TestCaller<test_fn_eigenvalue_computation>(
		"test_eigrl", &test_fn_eigenvalue_computation::test_eigrl));
	suite->addTest(new CppUnit::TestCaller<test_fn_eigenvalue_computation>(
		"test_eigh", &test_fn_eigenvalue_computation::test_eigh));
	suite->addTest(new CppUnit::TestCaller<test_fn_eigenvalue_computation>(
		"test_eighv", &test_fn_eigenvalue_computation::test_eighv));
	
	return suite;
}
