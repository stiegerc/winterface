// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tFn_mat_gen.h"
#include "test_tFn_mat_gen.cc"


// tests
template<>
void test_tFn_mat_gen<CPX__,RE__,CPX__>::test_rand() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	{
		const auto tMat1 = rand<tMat>(M);
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(M,tMat1.N());
	}
	{
		const RE__ l = 0.0;
		const RE__ u = 1.0;
		const auto tMat1 = rand<tMat>(M,N,l,u);
		
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
		for (const auto i: tMat1) {
			CPPUNIT_ASSERT(std::real(i)>=l && std::real(i)<=u);
			CPPUNIT_ASSERT(std::imag(i)>=l && std::imag(i)<=u);
		}
	}
	{
		const RE__ ll_ = -1.0;
		const RE__ ul_ = 0.0;
		const RE__ lu_ = 1.0;
		const RE__ uu_ = 2.0;
	
		const auto l = rand<fMat>(M,N,ll_,ul_);
		const auto u = rand<fMat>(M,N,lu_,uu_);
		
		const auto tMat1 = rand<tMat>(l,u);
		
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
		for (size_t i=0; i!=tMat1.L(); ++i) {
			CPPUNIT_ASSERT(std::real(tMat1[i])>=l[i] && std::real(tMat1[i])<=u[i]);
			CPPUNIT_ASSERT(std::imag(tMat1[i])>=l[i] && std::imag(tMat1[i])<=u[i]);
		}
	}
}

template<>
void test_tFn_mat_gen<CPX__,RE__,CPX__>::test_randi() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	{
		const auto tMat1 = randi<tMat>(M);
		
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(M,tMat1.N());
		for (const auto i: tMat1) {
			CPPUNIT_ASSERT_EQUAL(std::round(std::real(i)),std::real(i));
			CPPUNIT_ASSERT_EQUAL(std::round(std::imag(i)),std::imag(i));
		}
	}
	{
		const long l=0;
		const long u=100;
		const auto tMat1 = randi<tMat>(M,N);
		
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
		for (const auto i: tMat1) {
			CPPUNIT_ASSERT_EQUAL(std::round(std::real(i)),std::real(i));
			CPPUNIT_ASSERT(std::real(i)>=l && std::real(i)<=u);
			CPPUNIT_ASSERT_EQUAL(std::round(std::imag(i)),std::imag(i));
			CPPUNIT_ASSERT(std::imag(i)>=l && std::imag(i)<=u);
		}
	}
	{
		const long ll_ = -100;
		const long ul_ = 0;
		const long lu_ = 1;
		const long uu_ = 100;
	
		const auto l = randi<fMat>(M,N,ll_,ul_);
		const auto u = randi<fMat>(M,N,lu_,uu_);
		
		const auto tMat1 = randi<tMat>(l,u);
		
		CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
		CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
		for (size_t i=0; i!=tMat1.L(); ++i) {
			CPPUNIT_ASSERT_EQUAL(std::round(std::real(tMat1[i])),std::real(tMat1[i]));
			CPPUNIT_ASSERT(std::real(tMat1[i])>=l[i] && std::real(tMat1[i])<=u[i]);
			CPPUNIT_ASSERT_EQUAL(std::round(std::imag(tMat1[i])),std::imag(tMat1[i]));
			CPPUNIT_ASSERT(std::imag(tMat1[i])>=l[i] && std::imag(tMat1[i])<=u[i]);
		}
	}
}


// test id
template<>
const char* test_tFn_mat_gen<CPX__,RE__,CPX__>::test_id() noexcept {
	return "test_cFn_mat_gen";
}

// instantiation
template class test_tFn_mat_gen<CPX__,RE__,CPX__>;
