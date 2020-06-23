// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_CELL_COMPARISON_
#define _TEST_CELL_COMPARISON_

#include "ll_cell.h"
#include "ll_testTools.h"
#include <cppunit/extensions/HelperMacros.h>

class test_cell_comparison: public CppUnit::TestFixture {
public:
	// types
	typedef ll_cell cell;
	
public:
	void test_getPmat();
	void test_getAvec();
	void test_sameLattice();
	void test_operator_equal_unequal();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;

protected:
	struct cells {
		ll_cell cell1;
		ll_cell cell2;
	};
	inline cells mut_excl_cells(const size_t D) const noexcept {
		cells res = {ll__::test::genRandom(D), ll__::test::genRandom(D)};
		while (res.cell1==res.cell2)
			res.cell2 = ll__::test::genRandom(D);
		return res;
	}
};

#endif // _TEST_CELL_COMPARISON_
