// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tMat_ctor.h"
#include "test_tMat_ctor.cc"
#include "lm_fn.h"


using namespace lm__::test;

// tests
template<>
void test_tMat_ctor<CPX__,RE__,CPX__>::test_file() {
	CPPUNIT_ASSERT_THROW(cMat("data/bla.mat"),std::invalid_argument);	
	CPPUNIT_ASSERT_THROW(cMat("data/var_wc.mat"),std::invalid_argument);	
	CPPUNIT_ASSERT_THROW(cMat("data/alien_sym.mat"),std::invalid_argument);
	CPPUNIT_ASSERT_THROW(cMat("data/bad_i.mat"),std::invalid_argument);
	CPPUNIT_ASSERT_THROW(cMat("data/bad_format1.mat"),std::invalid_argument);		
	CPPUNIT_ASSERT_THROW(cMat("data/bad_format2.mat"),std::invalid_argument);
	CPPUNIT_ASSERT_NO_THROW(cMat("data/generic.mat"));	
	CPPUNIT_ASSERT_NO_THROW(cMat("data/sup_del.mat"));

	{
		const cMat tMat1("data/generic.mat");
		const cMat tMat2({CPX__(1.2,0.0),CPX__(0.0,1.0),CPX__(9.8,0.0),
				 CPX__(3.4,1.2),CPX__(0.0,-1.0),CPX__(-3.4,0.0),
				 CPX__(0.0,1.0),CPX__(0.0,1.4),CPX__(1.2,1.0)},3);
		CPPUNIT_ASSERT(tMat1==tMat2);
	}
	{
		const cMat tMat1("data/sup_del.mat");
		const cMat tMat2({CPX__(1.2,0.0),CPX__(0.0,1.0),CPX__(9.8,0.0),
				 CPX__(3.4,1.2),CPX__(0.0,-1.0),CPX__(-3.4,0.0),
				 CPX__(0.0,1.0),CPX__(0.0,1.4),CPX__(1.2,1.0)},3);
		CPPUNIT_ASSERT(tMat1==tMat2);
	}
	{
		const size_t M = genRndST();
		const size_t N = genRndST();
		const auto tMat1 = rnd<cMat>(M,N);

		tMat1.printToFile("data/rndtest.mat",8);
		const cMat tMat2("data/rndtest.mat");
		
		CPPUNIT_ASSERT(size(tMat1)==size(tMat2));
		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_DELTA(tMat1[i],tMat2[i],RE__(1e-6));
	}
	
	// binary
	{
		const size_t M = genRndST();
		const size_t N = genRndST();
		
		// write non 2d file and check exception
		{
			// open file
			std::ofstream file;
			file.open("data/bad.bin",std::ios::binary);
			if (!file.good())
				throw(std::invalid_argument("open file \'data/bad.bin\' failed"));

			// write header
			const std::string hdr ="cMat";
			file.write(hdr.c_str(),5);

			// write dimensions
			uint64_t dim;
			do dim=genRndST();
			while (dim==2);
			file.write((char*) &dim, sizeof(uint64_t));
			file.write((char*) &M, sizeof(uint64_t));
			file.write((char*) &N, sizeof(uint64_t));
			
			file.close();

			// check exception
			CPPUNIT_ASSERT_THROW(fMat("data/bad.bin"),std::invalid_argument);
		}
		{
			const fMat tMat1 = rnd<fMat>(M,N);
			tMat1.writeToFile("data/rndtest.bin");
			const cMat tMat2("data/rndtest.bin");
			CPPUNIT_ASSERT(tMat1==lm__::real(tMat2));
			CPPUNIT_ASSERT(lm__::zeros<fMat>(M,N)==lm__::imag(tMat2));
		}
		{
			const cMat tMat1 = rnd<cMat>(M,N);
			tMat1.writeToFile("data/rndtest.bin");
			const cMat tMat2("data/rndtest.bin");
			CPPUNIT_ASSERT(tMat1==tMat2);
		}
	}
}
template<>
void test_tMat_ctor<CPX__,RE__,CPX__>::test_fstream() {
	{
		std::ifstream file;
		file.open("data/generic.mat");
		CPPUNIT_ASSERT_THROW(cMat(file,4,3),std::invalid_argument);
		file.close();
	}
	{
		std::ifstream file;
		file.open("data/generic.mat");
		CPPUNIT_ASSERT_THROW(cMat(file,3,4),std::invalid_argument);
		file.close();
	}
	{
		std::ifstream file;
		file.open("data/ifstream_wc.mat");

		CPPUNIT_ASSERT_THROW(cMat(file,3,3),std::invalid_argument);
		file.close();
	}
}
template<>
void test_tMat_ctor<CPX__,RE__,CPX__>::test_re_im() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	const fMat re = rnd<fMat>(M,N);
	const fMat im = rnd<fMat>(M,N);

	// from Row
	{
		const auto ri = re.crBegin()+genRndST(0,M-1);
		const auto ii = im.crBegin()+genRndST(0,M-1);
		const cMat tMat(*ri,*ii);
		for (size_t i=0; i!=tMat.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__((*ri)[i],(*ii)[i]),tMat[i]);
	}

	// from Col
	{
		const auto ri = re.ccBegin()+genRndST(0,N-1);
		const auto ii = im.ccBegin()+genRndST(0,N-1);
		const cMat tMat(*ri,*ii);
		for (size_t i=0; i!=tMat.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__((*ri)[i],(*ii)[i]),tMat[i]);
	}

	// from Mat
	{
		const cMat tMat(re,im);
		for (size_t i=0; i!=tMat.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(CPX__(re[i],im[i]),tMat[i]);
	}
}
template<>
void test_tMat_ctor<CPX__,RE__,CPX__>::test_copy_fArray() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	const auto tMat1 = rnd<fMat>(M,N);
	
	// from Row
	{
		const auto i = tMat1.crBegin()+genRndST(0,M-1);
		const cMat tMat2(*i);
		
		CPPUNIT_ASSERT(msize(*i)==msize(tMat2));
		for (size_t j=0; j!=tMat2.L(); ++j)
			CPPUNIT_ASSERT_EQUAL((*i)[j],std::real(tMat2[j]));
	}
	
	// from Col
	{
		const auto i = tMat1.ccBegin()+genRndST(0,N-1);
		const cMat tMat2(*i);
		
		CPPUNIT_ASSERT(size(*i)==size(tMat2));
		for (size_t j=0; j!=tMat2.L(); ++j)
			CPPUNIT_ASSERT_EQUAL((*i)[j],std::real(tMat2[j]));
	}

	// from Mat
	{
		const cMat tMat2(tMat1);
		
		CPPUNIT_ASSERT(msize(tMat2)==msize(tMat1));
		for (size_t j=0; j!=tMat2.L(); ++j)
			CPPUNIT_ASSERT_EQUAL(tMat1[j],std::real(tMat2[j]));
	}
}
template<>
void test_tMat_ctor<CPX__,RE__,CPX__>::test_copy_cArray() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	const auto tMat1 = rnd<cMat>(M,N);
	
	// from Row
	{
		const auto i = tMat1.crBegin()+genRndST(0,M-1);
		const cMat tMat2(*i);
		CPPUNIT_ASSERT((*i)==tMat2);
	}

	// from Col
	{
		const auto i = tMat1.ccBegin()+genRndST(0,N-1);
		const cMat tMat2(*i);
		CPPUNIT_ASSERT((*i)==tMat2);
	}

	// from Mat
	{
		const cMat tMat2(tMat1);
		CPPUNIT_ASSERT(tMat1==tMat2);
	}
}


// test id
template<>
const char* test_tMat_ctor<CPX__,RE__,CPX__>::test_id() noexcept {
	return "test_cMat_ctor";
}


// instantiation
template class test_tMat_ctor<CPX__,RE__,CPX__>;
