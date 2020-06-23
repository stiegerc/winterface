// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_COMPOUND_ALL_
#define _TEST_COMPOUND_ALL_

#include "ll_types.h"
#include "ll_compound.h"
#include <cppunit/extensions/HelperMacros.h>

class test_compound_all: public CppUnit::TestFixture {
public:
	// types
	typedef ll__::mat_cb<lm__::fMat,ll__::fArray> mat_cb;
	typedef ll__::mat_b<lm__::fMat,lm__::fArray> mat_b;
	typedef ll__::vec_cb<std::vector<size_t>,size_t> vec_cb;
	typedef ll__::vec_b<std::vector<size_t>,size_t> vec_b;
	typedef ll__::mat_vec_cb<ll__::fMat,ll__::fArray,std::vector<size_t>,size_t> mat_vec_cb;
	typedef ll__::mat_vec_b<ll__::fMat,ll__::fArray,std::vector<size_t>,size_t> mat_vec_b;

	void test_mat_cb();
	void test_mat_b();
	void test_vec_cb();
	void test_vec_b();
	void test_mat_vec_cb();
	void test_mat_vec_b();
	void test_compound_classes();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_COMPOUND_ALL_
