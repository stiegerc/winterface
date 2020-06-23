// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tFn_math_functions.h"
#include "test_tFn_math_functions.cc"


// tests
template<>
void test_tFn_math_functions<CPX__,RE__,CPX__>::test_ceilEq_ceil() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// ceilEq on tMat
	{
		auto tMat1 = rnd<tMat>(M,N,.1,.9);
		const auto tMat2 = tMat1;
		ceilEq(tMat1);

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::ceil(std::real(tMat2[i])),
					std::ceil(std::imag(tMat2[i]))),tMat1[i]);

	}
	
	// ceilEq on tRow
	{
		auto tMat1 = rnd<tMat>(1,N,.1,.9);
		const auto tMat2 = tMat1;
		ceilEq(*tMat1.rBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::ceil(std::real(tMat2[i])),
					std::ceil(std::imag(tMat2[i]))),tMat1[i]);
	}
	
	// ceilEq on tCol
	{
		auto tMat1 = rnd<tMat>(M,1,.1,.9);
		const auto tMat2 = tMat1;
		ceilEq(*tMat1.cBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::ceil(std::real(tMat2[i])),
					std::ceil(std::imag(tMat2[i]))),tMat1[i]);
	}
	
	// ceil on tMat
	{
		auto tMat1 = rnd<tMat>(M,N,.1,.9);
		const auto tMat2 = ceil(tMat1);

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::ceil(std::real(tMat1[i])),
					std::ceil(std::imag(tMat1[i]))),tMat2[i]);
	}
	
	// ceil on tRow
	{
		auto tMat1 = rnd<tMat>(1,N,.1,.9);
		const auto tMat2 = ceil(*tMat1.crBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::ceil(std::real(tMat1[i])),
					std::ceil(std::imag(tMat1[i]))),tMat2[i]);
	}

	// ceil on tRow
	{
		auto tMat1 = rnd<tMat>(M,1,.1,.9);
		const auto tMat2 = ceil(*tMat1.ccBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::ceil(std::real(tMat1[i])),
					std::ceil(std::imag(tMat1[i]))),tMat2[i]);
	}

#ifndef NTOLERANT__
	{	
		auto tMat1 = rnd<tMat>(2,2);
		tMat1[0] = {-RE__(2.0)*mtol(),RE__(2.0)*mtol()};
		tMat1[1] = {RE__(2.0)*mtol(),-RE__(2.0)*mtol()};
		tMat1[2] = {-RE__(.5)*mtol(),RE__(.5)*mtol()};
		tMat1[3] = {RE__(.5)*mtol(),-RE__(.5)*mtol()};

		const auto tMat2 = ceil(tMat1);
		CPPUNIT_ASSERT_EQUAL(CPX__(0.0,1.0),tMat2[0]);
		CPPUNIT_ASSERT_EQUAL(CPX__(1.0,0.0),tMat2[1]);
		CPPUNIT_ASSERT_EQUAL(CPX__(0.0,0.0),tMat2[2]);
		CPPUNIT_ASSERT_EQUAL(CPX__(0.0,0.0),tMat2[3]);
	}
#endif
}

template<>
void test_tFn_math_functions<CPX__,RE__,CPX__>::test_floorEq_floor() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// floorEq on tMat
	{
		auto tMat1 = rnd<tMat>(M,N,.1,.9);
		const auto tMat2 = tMat1;
		floorEq(tMat1);

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::floor(std::real(tMat2[i])),
					std::floor(std::imag(tMat2[i]))),tMat1[i]);
	}
	
	// floorEq on tRow
	{
		auto tMat1 = rnd<tMat>(1,N,.1,.9);
		const auto tMat2 = tMat1;
		floorEq(*tMat1.rBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::floor(std::real(tMat2[i])),
					std::floor(std::imag(tMat2[i]))),tMat1[i]);
	}
	
	// floorEq on tCol
	{
		auto tMat1 = rnd<tMat>(M,1,.1,.9);
		const auto tMat2 = tMat1;
		floorEq(*tMat1.cBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::floor(std::real(tMat2[i])),
					std::floor(std::imag(tMat2[i]))),tMat1[i]);
	}
	
	// floor on tMat
	{
		auto tMat1 = rnd<tMat>(M,N,.1,.9);
		const auto tMat2 = floor(tMat1);

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::floor(std::real(tMat1[i])),
					std::floor(std::imag(tMat1[i]))),tMat2[i]);
	}
	
	// floor on tRow
	{
		auto tMat1 = rnd<tMat>(1,N,.1,.9);
		const auto tMat2 = floor(*tMat1.crBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::floor(std::real(tMat1[i])),
					std::floor(std::imag(tMat1[i]))),tMat2[i]);
	}

	// floor on tRow
	{
		auto tMat1 = rnd<tMat>(M,1,.1,.9);
		const auto tMat2 = floor(*tMat1.ccBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::floor(std::real(tMat1[i])),
					std::floor(std::imag(tMat1[i]))),tMat2[i]);
	}

#ifndef NTOLERANT__
	{
		auto tMat1 = rnd<tMat>(2,2);
		tMat1[0] = {-RE__(2.0)*mtol(),RE__(1.0)-RE__(2.0)*mtol()};
		tMat1[1] = {RE__(1.0)-RE__(2.0)*mtol(),-RE__(2.0)*mtol()};
		tMat1[2] = {-RE__(.5)*mtol(),RE__(.5)*mtol()};
		tMat1[3] = {RE__(.5)*mtol(),-RE__(.5)*mtol()};

		const auto tMat2 = floor(tMat1);
		CPPUNIT_ASSERT_EQUAL(CPX__(-1.0,0.0),tMat2[0]);
		CPPUNIT_ASSERT_EQUAL(CPX__(0.0,-1.0),tMat2[1]);
		CPPUNIT_ASSERT_EQUAL(CPX__(0.0,0.0),tMat2[2]);
		CPPUNIT_ASSERT_EQUAL(CPX__(0.0,0.0),tMat2[3]);
	}
#endif
}

template<>
void test_tFn_math_functions<CPX__,RE__,CPX__>::test_roundEq_round() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// roundEq on tMat
	{
		auto tMat1 = rnd<tMat>(M,N,.1,.9);
		const auto tMat2 = tMat1;
		roundEq(tMat1);

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::round(std::real(tMat2[i])),
					std::round(std::imag(tMat2[i]))),tMat1[i]);
	}
	
	// roundEq on tRow
	{
		auto tMat1 = rnd<tMat>(1,N,.1,.9);
		const auto tMat2 = tMat1;
		roundEq(*tMat1.rBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::round(std::real(tMat2[i])),
					std::round(std::imag(tMat2[i]))),tMat1[i]);
	}
	
	// roundEq on tCol
	{
		auto tMat1 = rnd<tMat>(M,1,.1,.9);
		const auto tMat2 = tMat1;
		roundEq(*tMat1.cBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::round(std::real(tMat2[i])),
					std::round(std::imag(tMat2[i]))),tMat1[i]);
	}
	
	// round on tMat
	{
		auto tMat1 = rnd<tMat>(M,N,.1,.9);
		const auto tMat2 = round(tMat1);

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::round(std::real(tMat1[i])),
					std::round(std::imag(tMat1[i]))),tMat2[i]);
	}
	
	// round on tRow
	{
		auto tMat1 = rnd<tMat>(1,N,.1,.9);
		const auto tMat2 = round(*tMat1.crBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::round(std::real(tMat1[i])),
					std::round(std::imag(tMat1[i]))),tMat2[i]);
	}

	// round on tRow
	{
		auto tMat1 = rnd<tMat>(M,1,.1,.9);
		const auto tMat2 = round(*tMat1.ccBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(std::round(std::real(tMat1[i])),
					std::round(std::imag(tMat1[i]))),tMat2[i]);
	}
}

template<>
void test_tFn_math_functions<CPX__,RE__,CPX__>::test_signEq_sign() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();
	
	const auto P = RE__(2.0)*round(rnd<tMat>(M,N,RE__(-.4999),RE__(1.4999)))-RE__(1.0);

	const auto sgn_ = [](const CPX__& i) -> CPX__ {
			return std::abs(i)<=RE__(0.0) ? CPX__(0.0): i/std::abs(i);
	};

	// signEq on tMat
	{
		auto tMat1 = rnd<tMat>(M,N,RE__(.1),RE__(.9))*=P;
		const auto tMat2 = tMat1;
		signEq(tMat1);

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(sgn_(tMat2[i]),tMat1[i]);

	}
	
	// signEq on tRow
	{
		auto tMat1 = rnd<tMat>(1,N,RE__(.1),RE__(.9))*=P.rGet(0);
		const auto tMat2 = tMat1;
		signEq(*tMat1.rBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(sgn_(tMat2[i]),tMat1[i]);

	}
	
	// signEq on tCol
	{
		auto tMat1 = rnd<tMat>(M,1,RE__(.1),RE__(.9))*=P.cGet(0);
		const auto tMat2 = tMat1;
		signEq(*tMat1.cBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(sgn_(tMat2[i]),tMat1[i]);

	}
	
	// sign on tMat
	{
		auto tMat1 = rnd<tMat>(M,N,RE__(.1),RE__(.9))*=P;
		const auto tMat2 = sign(tMat1);

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(sgn_(tMat1[i]),tMat2[i]);

	}
	
	// sign on tRow
	{
		auto tMat1 = rnd<tMat>(1,N,RE__(.1),RE__(.9))*=P.rGet(0);
		const auto tMat2 = sign(*tMat1.crBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(sgn_(tMat1[i]),tMat2[i]);

	}

	// sign on tCol
	{
		auto tMat1 = rnd<tMat>(M,1,.1,.9)*=P.cGet(0);
		const auto tMat2 = sign(*tMat1.ccBegin());

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(sgn_(tMat1[i]),tMat2[i]);

	}
#ifndef NTOLERANT__
	{
		auto tMat1 = rnd<tMat>(2,2);
		tMat1[0] = {-RE__(2.0)*mtol(),RE__(2.0)*mtol()};
		tMat1[1] = {RE__(2.0)*mtol(),-RE__(2.0)*mtol()};
		tMat1[2] = {-RE__(.5)*mtol(),RE__(.5)*mtol()};
		tMat1[3] = {RE__(.5)*mtol(),-RE__(.5)*mtol()};

		const auto tMat2 = sign(tMat1);
		CPPUNIT_ASSERT(std::abs(tMat2[0]/std::abs(tMat2[0])-tMat2[0])<delta);
		CPPUNIT_ASSERT(std::abs(tMat2[1]/std::abs(tMat2[1])-tMat2[1])<delta);
		CPPUNIT_ASSERT_EQUAL(CPX__(0.0,0.0),tMat2[2]);
		CPPUNIT_ASSERT_EQUAL(CPX__(0.0,0.0),tMat2[3]);
	}
#endif
}


// test id
template<>
const char* test_tFn_math_functions<CPX__,RE__,CPX__>::test_id() noexcept {
	return "test_cFn_math_functions";
}

// instantiation
template class test_tFn_math_functions<CPX__,RE__,CPX__>;
