// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_MESH_CTOR_
#define _TEST_MESH_CTOR_

#include "ll_mesh.h"
#include <cppunit/extensions/HelperMacros.h>

class test_mesh_all: public CppUnit::TestFixture {
public:
	void test_itr();
	void test_c_itr();
	void test_mesh();
	void test_generators();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_MESH_CTOR_
