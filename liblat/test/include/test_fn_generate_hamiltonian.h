// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_FN_GENERATE_HAMILTONIAN_
#define _TEST_FN_GENERATE_HAMILTONIAN_

#include <cppunit/extensions/HelperMacros.h>

class test_fn_generate_hamiltonian: public CppUnit::TestFixture {
public:
	void test_genHam_from_trans();
	void test_genHam_from_wbh();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_FN_GENERATE_HAMILTONIAN_
