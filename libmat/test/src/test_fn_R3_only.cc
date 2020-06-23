// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_fn_R3_only.h"
#include "testTools.h"
#include "lm_fn.h"
#include <iostream>

using namespace lm__;
using namespace lm__::test;


// tests
void test_fn_R3_only:: test_cross_crossn() {
	const auto a = rand<fMat>(3,1);
	const auto b = rand<fMat>(3,1);
	const auto c = cross(a,b);

	CPPUNIT_ASSERT_DELTA(dot(a,c),0.0,mtol());
	CPPUNIT_ASSERT_DELTA(dot(b,c),0.0,mtol());

	const double phi = std::acos(dot(a,b)/norm(a)/norm(b));
	CPPUNIT_ASSERT_DELTA(std::sin(phi)*norm(a)*norm(b),norm(c),mtol());

	const auto d = crossn(a,b);
	CPPUNIT_ASSERT(c==(d*norm(c)));
	CPPUNIT_ASSERT_DELTA(1.0,norm(d),mtol());
}
void test_fn_R3_only::test_getR() {
	const auto a = rand<fMat>(3,1);
	const auto b = rand<fMat>(3,1);
	
	// compare version 1 and 2
	{
		const auto n = cross(a,b)/norm(cross(a,b));
		const double phi = std::acos(dot(a,b)/norm(a)/norm(b));

		const auto R1 = getR(phi,n);
		CPPUNIT_ASSERT_DELTA(det(R1),1.0,mtol());

		const auto c = R1.prod(a);
		CPPUNIT_ASSERT_DELTA(norm(a),norm(c),delta);
		
		const auto cn = c/norm(c);
		const auto bn = b/norm(b);
		for (size_t i=0; i!=3; ++i)
			CPPUNIT_ASSERT_DELTA(bn[i],cn[i],delta);

		const auto R2 = getR(a,b);
		CPPUNIT_ASSERT(R1==R2);
	}

	// compare version 1 to explicit definitions around z-axis
	{
		const RE__ phi = rnd<fMat>(1,1,0,2*M_PI)[0];
		const fMat Rz({std::cos(phi),std::sin(phi),0.0,
			      -std::sin(phi),std::cos(phi),0.0,
			      0.0,0.0,1.0},3);
		const fMat R = getR(phi,fMat({0.0,0.0,1.0},1));
		CPPUNIT_ASSERT(R==Rz);
	}
}


// test id
const char* test_fn_R3_only::test_id() noexcept {
	return "test_fn_R3_only";
}

CppUnit::Test* test_fn_R3_only::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_fn_R3_only>(
		"test_cross_crossn", &test_fn_R3_only::test_cross_crossn));
	suite->addTest(new CppUnit::TestCaller<test_fn_R3_only>(
		"test_getR", &test_fn_R3_only::test_getR));
	
	return suite;
}
