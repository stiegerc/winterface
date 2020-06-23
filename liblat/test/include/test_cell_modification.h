// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_CELL_MODIFICATION_
#define _TEST_CELL_MODIFICATION_

#include "ll_cell.h"
#include <algorithm>
#include <cppunit/extensions/HelperMacros.h>

class test_cell_modification: public CppUnit::TestFixture {
public:
	void test_stress();
	void test_swapDim();
	void test_invDim();
	void test_orient();
	void test_permute();
	void test_rotate();
	void test_scale();
	void test_changeBasis();
	void test_makePrimitive();
	void test_shift();
	void test_diversify();
	void test_collectivize();
	void test_merge();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;

protected:
	// helper function to check whether Ap is sorted
	bool Ap_sorted(const ll_cell& tCell) const noexcept {
		for (const auto t: tCell.types())
			if (!std::is_sorted(tCell.ccBegin(t),tCell.ccEnd(t),lm__::vcmp))
				return false;
		return true;
	}
};

#endif // _TEST_CELL_MODIFICATION_
