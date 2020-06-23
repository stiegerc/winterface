// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_ctor.h"
#include "testTools.h"
#include <iostream>


using namespace lm__::test;

template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_mn() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	const tMat tMat1(M,N);
	CPPUNIT_ASSERT(!tMat1.empty());
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
	CPPUNIT_ASSERT_EQUAL(M*N,tMat1.L());
	CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat1.ccap());
	CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat1.lcap());

	const tMat tMat2(M,0);
	CPPUNIT_ASSERT(tMat2.empty());
	CPPUNIT_ASSERT_EQUAL(M,tMat2.M());
	CPPUNIT_ASSERT_EQUAL(size_t(0),tMat2.N());
	CPPUNIT_ASSERT_EQUAL(size_t(0),tMat2.L());
	CPPUNIT_ASSERT_EQUAL(size_t(0),tMat2.ccap());
	CPPUNIT_ASSERT_EQUAL(size_t(0),tMat2.lcap());
	
	const tMat tMat3(0,N);
	CPPUNIT_ASSERT(tMat3.empty());
	CPPUNIT_ASSERT_EQUAL(size_t(0),tMat3.M());
	CPPUNIT_ASSERT_EQUAL(N,tMat3.N());
	CPPUNIT_ASSERT_EQUAL(size_t(0),tMat3.L());
	CPPUNIT_ASSERT_EQUAL(size_t(0),tMat3.ccap());
	CPPUNIT_ASSERT_EQUAL(size_t(0),tMat3.lcap());
}

template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_m() {
	const size_t M = genRndST();
	const tMat tMat1(M);
		
	CPPUNIT_ASSERT(!tMat1.empty());
	CPPUNIT_ASSERT(tMat1.square());
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(M,tMat1.N());
	CPPUNIT_ASSERT_EQUAL(M*M,tMat1.L());
	CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat1.ccap());
	CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat1.lcap());
}


template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_default() {
	const tMat tMat1;
		
	CPPUNIT_ASSERT(tMat1.empty());
	CPPUNIT_ASSERT(!tMat1.M() && !tMat1.N() && !tMat1.L());
	CPPUNIT_ASSERT(!tMat1.lcap());
	CPPUNIT_ASSERT(!tMat1.ccap());
}

template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_prealloc() {
	const size_t M = genRndST();
	const size_t N = genRndST();
	TT* data_ = new TT [M*N];
		
	const tMat tMat1(data_,M,N);
	
	CPPUNIT_ASSERT(!tMat1.empty());
	CPPUNIT_ASSERT_EQUAL(M,tMat1.M());
	CPPUNIT_ASSERT_EQUAL(N,tMat1.N());
	CPPUNIT_ASSERT_EQUAL(M*N,tMat1.L());
	CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat1.ccap());
	CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(const_cast<const TT*>(data_),tMat1.data());
}

template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_from_itrs() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	const auto tMat1 = rnd<tMat>(M,N);
	
	// row itrs
	{
		const size_t m1 = genRndST(0,M-1);
		const size_t m2 = genRndST(m1,M-1);
		
		const tMat tMat2(tMat1.rBegin()+m1,tMat1.rBegin()+m2);
		CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat2.N());
		CPPUNIT_ASSERT_EQUAL(m2-m1,tMat2.M());
		for (auto i=tMat2.rBegin(),e=tMat2.rEnd(),j=tMat1.rBegin()+m1; i!=e; ++i,++j)
			CPPUNIT_ASSERT_EQUAL(*i,*j);
		
		const tMat tMat3(tMat1.crBegin()+m1,tMat1.crBegin()+m2);
		CPPUNIT_ASSERT(tMat2==tMat3);
	}
	
	// col itrs
	{
		const size_t n1 = genRndST(0,N-1);
		const size_t n2 = genRndST(n1,N-1);
		
		const tMat tMat2(tMat1.cBegin()+n1,tMat1.cBegin()+n2);
		CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.M());
		CPPUNIT_ASSERT_EQUAL(n2-n1,tMat2.N());
		for (auto i=tMat2.cBegin(),e=tMat2.cEnd(),j=tMat1.cBegin()+n1; i!=e; ++i,++j)
			CPPUNIT_ASSERT_EQUAL(*i,*j);
		
		const tMat tMat3(tMat1.ccBegin()+n1,tMat1.ccBegin()+n2);
		CPPUNIT_ASSERT(tMat2==tMat3);
	}
}

template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_initl_row_col() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	// from rows
	{
		const auto tMat1 = rnd<tMat>(6,N);
		std::vector<size_t> I = {0,1,2,3,4,5};
		I.resize(3);
		std::random_shuffle(I.begin(),I.end());

		const tMat tMat2 = {tMat1.rAt(I[0]),tMat1.rAt(I[1]),tMat1.rAt(I[2])};

		CPPUNIT_ASSERT_EQUAL(size_t(3),tMat2.M());
		CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat2.N());
		CPPUNIT_ASSERT_EQUAL(tMat2.M()*tMat2.N(),tMat2.lcap());
		CPPUNIT_ASSERT(tMat1.rAt(I[0])==tMat2.rAt(0));
		CPPUNIT_ASSERT(tMat1.rAt(I[1])==tMat2.rAt(1));
		CPPUNIT_ASSERT(tMat1.rAt(I[2])==tMat2.rAt(2));
	}
	
	// from cols
	{
		const auto tMat1 = rnd<tMat>(M,6);
		std::vector<size_t> I = {0,1,2,3,4,5};
		I.resize(3);
		std::random_shuffle(I.begin(),I.end());

		const tMat tMat2 = {tMat1.cAt(I[0]),tMat1.cAt(I[1]),tMat1.cAt(I[2])};

		CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.M());
		CPPUNIT_ASSERT_EQUAL(size_t(3),tMat2.N());
		CPPUNIT_ASSERT_EQUAL(tMat2.M()*tMat2.N(),tMat2.lcap());
		CPPUNIT_ASSERT(tMat1.cAt(I[0])==tMat2.cAt(0));
		CPPUNIT_ASSERT(tMat1.cAt(I[1])==tMat2.cAt(1));
		CPPUNIT_ASSERT(tMat1.cAt(I[2])==tMat2.cAt(2));
	}
}

template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_initl() {
	const tMat tMat1({TT(1.0),TT(2.0),TT(3.0),TT(4.0)});

	CPPUNIT_ASSERT(!tMat1.empty());
	CPPUNIT_ASSERT_EQUAL(size_t(4),tMat1.M());
	CPPUNIT_ASSERT_EQUAL(size_t(1),tMat1.N());
	CPPUNIT_ASSERT_EQUAL(size_t(4),tMat1.L());
	CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat1.ccap());
	CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(TT(1.0),tMat1[0]);
	CPPUNIT_ASSERT_EQUAL(TT(2.0),tMat1[1]);
	CPPUNIT_ASSERT_EQUAL(TT(3.0),tMat1[2]);
	CPPUNIT_ASSERT_EQUAL(TT(4.0),tMat1[3]);
}

template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_initl_m() {
	const tMat tMat1({TT(1.0),TT(2.0),TT(3.0),TT(4.0)},2);
		
	CPPUNIT_ASSERT(!tMat1.empty());
	CPPUNIT_ASSERT_EQUAL(size_t(2),tMat1.M());
	CPPUNIT_ASSERT_EQUAL(size_t(2),tMat1.N());
	CPPUNIT_ASSERT_EQUAL(size_t(4),tMat1.L());
	CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat1.ccap());
	CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(TT(1.0),tMat1[0]);
	CPPUNIT_ASSERT_EQUAL(TT(2.0),tMat1[1]);
	CPPUNIT_ASSERT_EQUAL(TT(3.0),tMat1[2]);
	CPPUNIT_ASSERT_EQUAL(TT(4.0),tMat1[3]);
}

template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_initl_mn() {
	const tMat tMat1({TT(1.0),TT(2.0),TT(3.0),TT(4.0)},2,2);
		
	CPPUNIT_ASSERT(!tMat1.empty());
	CPPUNIT_ASSERT_EQUAL(size_t(2),tMat1.M());
	CPPUNIT_ASSERT_EQUAL(size_t(2),tMat1.N());
	CPPUNIT_ASSERT_EQUAL(size_t(4),tMat1.L());
	CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat1.ccap());
	CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat1.lcap());
	CPPUNIT_ASSERT_EQUAL(TT(1.0),tMat1[0]);
	CPPUNIT_ASSERT_EQUAL(TT(2.0),tMat1[1]);
	CPPUNIT_ASSERT_EQUAL(TT(3.0),tMat1[2]);
	CPPUNIT_ASSERT_EQUAL(TT(4.0),tMat1[3]);
}

template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_vec_tt() {
	const size_t N = genRndST();
		
	std::vector<TT> inp(N);
	std::iota(inp.begin(),inp.end(),1);
	
	// to col
	{	
		const tMat tMat1(inp,true);
		CPPUNIT_ASSERT_EQUAL(inp.size(),tMat1.L());
		CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat1.M());
		for (size_t i=0; i<inp.size(); ++i)
			CPPUNIT_ASSERT_EQUAL(inp[i],tMat1[i]);
	}

	// to row
	{
		const tMat tMat1(inp,false);
		CPPUNIT_ASSERT_EQUAL(inp.size(),tMat1.L());
		CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat1.N());
		for (size_t i=0; i<inp.size(); ++i)
			CPPUNIT_ASSERT_EQUAL(inp[i],tMat1[i]);
	}
}

template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_vec_size_t() {
	const auto inp = genRndI(1,20,0,40);
	
	// to col
	{
		const tMat tMat1(inp,true);
		CPPUNIT_ASSERT_EQUAL(inp.size(),tMat1.L());
		CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat1.M());
		for (size_t i=0; i<inp.size(); ++i)
			CPPUNIT_ASSERT_EQUAL(TT(inp[i]),tMat1[i]);
	}
	
	// to row
	{
		const tMat tMat1(inp,false);
		CPPUNIT_ASSERT_EQUAL(inp.size(),tMat1.L());
		CPPUNIT_ASSERT_EQUAL(tMat1.L(),tMat1.N());
		for (size_t i=0; i<inp.size(); ++i)
			CPPUNIT_ASSERT_EQUAL(TT(inp[i]),tMat1[i]);
	}
}

template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_copy() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	const auto tMat1 = rnd<tMat>(M,N);
	const auto tMat2(tMat1);

	CPPUNIT_ASSERT_EQUAL(tMat1.M(),tMat2.M());
	CPPUNIT_ASSERT_EQUAL(tMat1.N(),tMat2.N());
	CPPUNIT_ASSERT_EQUAL(tMat2.L(),tMat2.lcap());
	CPPUNIT_ASSERT_EQUAL(tMat2.N(),tMat2.ccap());
	for (size_t i=0; i!=M*N; ++i)
		CPPUNIT_ASSERT_EQUAL(tMat1[i],tMat2[i]);
}

template<class TT, class FT, class CT>
void test_tMat_ctor<TT,FT,CT>::test_move() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	auto tMat1 = rnd<tMat>(M,N);
	const auto ptr = tMat1.data();

	const auto tMat2 = std::move(tMat1);
	CPPUNIT_ASSERT_EQUAL(M,tMat2.M());
	CPPUNIT_ASSERT_EQUAL(N,tMat2.N());
	CPPUNIT_ASSERT_EQUAL(const_cast<const TT*>(ptr),tMat2.data());
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tMat_ctor<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_mn", &test_tMat_ctor<TT,FT,CT>::test_mn));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_m", &test_tMat_ctor<TT,FT,CT>::test_m));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_default", &test_tMat_ctor<TT,FT,CT>::test_default));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_prealloc", &test_tMat_ctor<TT,FT,CT>::test_prealloc));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_from_itrs", &test_tMat_ctor<TT,FT,CT>::test_from_itrs));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_initl_row_col", &test_tMat_ctor<TT,FT,CT>::test_initl_row_col));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_initl", &test_tMat_ctor<TT,FT,CT>::test_initl));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_initl_m", &test_tMat_ctor<TT,FT,CT>::test_initl_m));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_initl_mn", &test_tMat_ctor<TT,FT,CT>::test_initl_mn));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_file", &test_tMat_ctor<TT,FT,CT>::test_file));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_fstream", &test_tMat_ctor<TT,FT,CT>::test_fstream));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_vec_tt", &test_tMat_ctor<TT,FT,CT>::test_vec_tt));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_vec_size_t", &test_tMat_ctor<TT,FT,CT>::test_vec_size_t));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_re_im", &test_tMat_ctor<TT,FT,CT>::test_re_im));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_copy_fArray", &test_tMat_ctor<TT,FT,CT>::test_copy_fArray));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_copy_cArray", &test_tMat_ctor<TT,FT,CT>::test_copy_cArray));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_copy", &test_tMat_ctor<TT,FT,CT>::test_copy));
	suite->addTest(new CppUnit::TestCaller<test_tMat_ctor>(
		"test_move", &test_tMat_ctor<TT,FT,CT>::test_move));
	
	return suite;
}
