// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_cpxItr_all.h"
#include "test_cpxItr_all.cc"


using namespace lm__;

// tests
template<>
void test_cpxItr_all<1>::test_all_dereference() {
	
	// c_imItr
	{
		// lambda to generate const random range
		auto gen = [](const size_t L) -> const CPX__* {
			CPX__* res = new CPX__ [L];
			rnd(res,L);
			return res;
		};

		const size_t L = genRndST();
		const CPX__* rg = gen(L);

		c_imItr tItr(rg,1);
		for (size_t i=0; i!=L; ++i,++tItr)
			CPPUNIT_ASSERT_EQUAL(std::imag(rg[i]),*tItr);

		delete[] rg;
	}

	// imItr
	{
		// lambda to generate const random range
		auto gen = [](const size_t L) -> CPX__* {
			CPX__* res = new CPX__ [L];
			rnd(res,L);
			return res;
		};

		const size_t L = genRndST();
		CPX__* rg = gen(L);

		imItr tItr(rg,1);
		for (size_t i=0; i!=L; ++i,++tItr) {
			CPPUNIT_ASSERT_EQUAL(std::imag(rg[i]),*tItr);
			
			const auto tmp = rg[i]; *tItr = RE__(1.0);
			CPPUNIT_ASSERT_EQUAL(std::imag(rg[i]),RE__(1.0));
			CPPUNIT_ASSERT_EQUAL(std::real(rg[i]),std::real(tmp));
		}

		delete[] rg;
	}
}

// test id
template<>
const char* test_cpxItr_all<1>::test_id() noexcept {
	return "test_imItr_all";
}

// instantiation
template class test_cpxItr_all<1>;
