// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TEST_IO_
#define _TEST_IO_

#include <cppunit/extensions/HelperMacros.h>
#include "ll_types.h"

class test_io: public CppUnit::TestFixture {
public:
	void test_writePOSCAR();
	void test_readEf();
	void test_genr();
	void test_readLayerMatrix();
	void test_printOmf();
	void test_printOlf();
	void test_readOmf();
	void test_readOlf();
	void test_hrDim();
	void test_readB();
	void test_readAp();
	void test_readWp();
	void test_readHr();
	void test_readAMN();
	void test_readXYZ();
	void test_readEig();
	void test_readWannierTransf();
	void test_readChk();

	static CppUnit::Test* suite();
	static const char* test_id() noexcept;

protected:
	ll__::fMat genRndP_(const size_t d) const noexcept;
};

#endif // _TEST_IO_
