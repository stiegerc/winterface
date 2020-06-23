// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tFn_functionals.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


template<class TT, class FT, class CT>
void test_tFn_functionals<TT,FT,CT>::test_sum() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,N);
		const TT ref = std::accumulate(tMat1.cbegin(),tMat1.cend(),TT(0.0));
		const TT act = sum(tMat1);
		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		for (auto i=tMat1.crBegin(), e=tMat1.crEnd(); i!=e; ++i) {
			const TT ref = std::accumulate(i->cbegin(),i->cend(),TT(0.0));
			const TT act = sum(*i);
			CPPUNIT_ASSERT_DELTA(ref,act,delta);
		}
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		for (auto i=tMat1.ccBegin(), e=tMat1.ccEnd(); i!=e; ++i) {
			const TT ref = std::accumulate(i->cbegin(),i->cend(),TT(0.0));
			const TT act = sum(*i);
			CPPUNIT_ASSERT_DELTA(ref,act,delta);
		}
	}
}

template<class TT, class FT, class CT>
void test_tFn_functionals<TT,FT,CT>::test_prod() {
	
	const size_t M = genRndST(2,4);
	const size_t N = genRndST(2,4);

	{
		const auto tMat1 = rnd<tMat>(M,N,1.0,1.001);
		const TT ref = std::accumulate(tMat1.cbegin(),tMat1.cend(),TT(1.0),std::multiplies<TT>());
		const TT act = prod(tMat1);
		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}
	{
		const auto tMat1 = rnd<tMat>(M,N,1.0,1.001);
		for (auto i=tMat1.crBegin(), e=tMat1.crEnd(); i!=e; ++i) {
			const TT ref = std::accumulate(i->cbegin(),i->cend(),TT(1.0),std::multiplies<TT>());
			const TT act = prod(*i);
			CPPUNIT_ASSERT_DELTA(ref,act,delta);
		}
	}
	{
		const auto tMat1 = rnd<tMat>(M,N,1.0,1.001);
		for (auto i=tMat1.ccBegin(), e=tMat1.ccEnd(); i!=e; ++i) {
			const TT ref = std::accumulate(i->cbegin(),i->cend(),TT(1.0),std::multiplies<TT>());
			const TT act = prod(*i);
			CPPUNIT_ASSERT_DELTA(ref,act,delta);
		}
	}
}

template<class TT, class FT, class CT>
void test_tFn_functionals<TT,FT,CT>::test_min() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		const auto tMat1 = rnd<tMat>(M,N);
		CPPUNIT_ASSERT_EQUAL(*std::min_element(tMat1.cbegin(),tMat1.cend(),
					[](const TT a, const TT b){return ops::lt_s(a,b);}),
					min(tMat1));
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		for (auto i=tMat1.crBegin(), e=tMat1.crEnd(); i!=e; ++i)
			CPPUNIT_ASSERT_EQUAL(*std::min_element(i->cbegin(),i->cend(),
					[](const TT a, const TT b){return ops::lt_s(a,b);}),
					min(*i));
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		for (auto i=tMat1.ccBegin(), e=tMat1.ccEnd(); i!=e; ++i)
			CPPUNIT_ASSERT_EQUAL(*std::min_element(i->cbegin(),i->cend(),
					[](const TT a, const TT b){return ops::lt_s(a,b);}),
					min(*i));
	}
}

template<class TT, class FT, class CT>
void test_tFn_functionals<TT,FT,CT>::test_max() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		const auto tMat1 = rnd<tMat>(M,N);
		CPPUNIT_ASSERT_EQUAL(*std::max_element(tMat1.cbegin(),tMat1.cend(),
					[](const TT a, const TT b){return ops::lt_s(a,b);}),
					max(tMat1));
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		for (auto i=tMat1.crBegin(), e=tMat1.crEnd(); i!=e; ++i)
			CPPUNIT_ASSERT_EQUAL(*std::max_element(i->cbegin(),i->cend(),
					[](const TT a, const TT b){return ops::lt_s(a,b);}),
					max(*i));
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		for (auto i=tMat1.cBegin(), e=tMat1.cEnd(); i!=e; ++i)
			CPPUNIT_ASSERT_EQUAL(*std::max_element(i->cbegin(),i->cend(),
					[](const TT a, const TT b){return ops::lt_s(a,b);}),
					max(*i));
	}
}

template<class TT, class FT, class CT>
void test_tFn_functionals<TT,FT,CT>::test_mean() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,N);
		const TT ref = std::accumulate(tMat1.cbegin(),tMat1.cend(),TT(0.0))/RE__(M*N);
		const TT act = mean(tMat1);
		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		for (auto i=tMat1.crBegin(), e=tMat1.crEnd(); i!=e; ++i) {
			const TT ref = std::accumulate(i->cbegin(),i->cend(),TT(0.0))/RE__(N);
			const TT act = mean(*i);
			CPPUNIT_ASSERT_DELTA(ref,act,delta);
		}
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		for (auto i=tMat1.ccBegin(), e=tMat1.ccEnd(); i!=e; ++i) {
			const TT ref = std::accumulate(i->cbegin(),i->cend(),TT(0.0))/RE__(M);
			const TT act = mean(*i);
			CPPUNIT_ASSERT_DELTA(ref,act,delta);
		}
	}
}

template<class TT, class FT, class CT>
void test_tFn_functionals<TT,FT,CT>::test_normsq() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,N);
		FT ref(0.0);
		for (auto j: tMat1) ref += std::abs(j)*std::abs(j);
		const FT act = normsq(tMat1);
		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		for (auto i=tMat1.crBegin(), e=tMat1.crEnd(); i!=e; ++i) {
			FT ref(0.0);
			for (auto j: *i) ref += std::abs(j)*std::abs(j);
			const FT act = normsq(*i);
			CPPUNIT_ASSERT_DELTA(ref,act,delta);
		}
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		for (auto i=tMat1.ccBegin(), e=tMat1.ccEnd(); i!=e; ++i) {
			FT ref(0.0);
			for (auto j: *i) ref += std::abs(j)*std::abs(j);
			const FT act = normsq(*i);
			CPPUNIT_ASSERT_DELTA(ref,act,delta);
		}
	}	
}

template<class TT, class FT, class CT>
void test_tFn_functionals<TT,FT,CT>::test_norm() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = rnd<tMat>(M,N);
		FT ref(0.0);
		for (auto j: tMat1) ref += std::abs(j)*std::abs(j);
		ref =  std::sqrt(ref);
		const FT act = norm(tMat1);
		CPPUNIT_ASSERT_DELTA(ref,act,delta);
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		for (auto i=tMat1.crBegin(), e=tMat1.crEnd(); i!=e; ++i) {
			FT ref(0.0);
			for (auto j: *i) ref += std::abs(j)*std::abs(j);
			ref =  std::sqrt(ref);
			const FT act = norm(*i);
			CPPUNIT_ASSERT_DELTA(ref,act,delta);
		}
	}
	{
		const auto tMat1 = rnd<tMat>(M,N);
		for (auto i=tMat1.ccBegin(), e=tMat1.ccEnd(); i!=e; ++i) {
			FT ref(0.0);
			for (auto j: *i) ref += std::abs(j)*std::abs(j);
			ref =  std::sqrt(ref);
			const FT act = norm(*i);
			CPPUNIT_ASSERT_DELTA(ref,act,delta);
		}
	}
}

template<class TT, class FT, class CT>
void test_tFn_functionals<TT,FT,CT>::test_trace() {
	const size_t M = genRndST();
	const auto tMat1 = rnd<tMat>(M,M);
	
	TT ref(0.0);
	for (size_t m=0; m!=M; ++m) ref+=tMat1(m,m);
	const auto act = trace(tMat1);

	CPPUNIT_ASSERT_DELTA(ref,act,delta);
}

template<class TT, class FT, class CT>
void test_tFn_functionals<TT,FT,CT>::test_det() {
	const size_t M = genRndST();
	auto tMat1 = rnd_b<tMat>(M,M);
	gsorth(tMat1);

	TT ref(1.0);
	const auto act = std::abs(det(tMat1));

	CPPUNIT_ASSERT_DELTA(ref,act,delta);
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tFn_functionals<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tFn_functionals>(
		"test_sum", &test_tFn_functionals<TT,FT,CT>::test_sum));
	suite->addTest(new CppUnit::TestCaller<test_tFn_functionals>(
		"test_prod", &test_tFn_functionals<TT,FT,CT>::test_prod));
	suite->addTest(new CppUnit::TestCaller<test_tFn_functionals>(
		"test_min", &test_tFn_functionals<TT,FT,CT>::test_min));
	suite->addTest(new CppUnit::TestCaller<test_tFn_functionals>(
		"test_max", &test_tFn_functionals<TT,FT,CT>::test_max));
	suite->addTest(new CppUnit::TestCaller<test_tFn_functionals>(
		"test_mean", &test_tFn_functionals<TT,FT,CT>::test_mean));
	suite->addTest(new CppUnit::TestCaller<test_tFn_functionals>(
		"test_normsq", &test_tFn_functionals<TT,FT,CT>::test_normsq));
	suite->addTest(new CppUnit::TestCaller<test_tFn_functionals>(
		"test_norm", &test_tFn_functionals<TT,FT,CT>::test_norm));
	suite->addTest(new CppUnit::TestCaller<test_tFn_functionals>(
		"test_trace", &test_tFn_functionals<TT,FT,CT>::test_trace));
	suite->addTest(new CppUnit::TestCaller<test_tFn_functionals>(
		"test_det", &test_tFn_functionals<TT,FT,CT>::test_det));
	
	return suite;
}
