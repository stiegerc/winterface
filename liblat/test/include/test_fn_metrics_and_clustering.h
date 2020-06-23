// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_FN_METRICS_AND_CLUSTERING_
#define _TEST_FN_METRICS_AND_CLUSTERING_

#include <cppunit/extensions/HelperMacros.h>

class test_fn_metrics_and_clustering: public CppUnit::TestFixture {
public:
	void test_genNNmat();
	void test_dist_distb();
	void test_genDmat();
	void test_com();
	void test_dbscan();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;
};

#endif // _TEST_FN_METRICS_AND_CLUSTERING_
