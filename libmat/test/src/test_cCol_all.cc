// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "testTools.h"

#include "test_tCol_all.h"
#include "test_tCol_all.cc"


template<>
void test_tCol_all<CPX__,RE__,CPX__>::test_assign() {

	const size_t M = genRndST();
	const size_t N = genRndST();
	
	// from real scalar
	{
		auto tMat1 = rnd<cMat>(M,N);
		auto f = rnd<fMat>(1,1)[0];

		const size_t n1 = genRndST(0,N-1);
		auto j = tMat1.cBegin()+n1;

		*j = f;
		CPPUNIT_ASSERT(std::all_of(j->begin(),j->end(),[&f](const CPX__ i){return i==f;}));
	}
	
	// from complex scalar
	{
		auto tMat1 = rnd<cMat>(M,N);
		auto f = rnd<cMat>(1,1)[0];

		const size_t n1 = genRndST(0,N-1);
		auto j = tMat1.cBegin()+n1;

		*j = f;
		CPPUNIT_ASSERT(std::all_of(j->begin(),j->end(),[&f](const CPX__ i){return i==f;}));
	}
	
	// from real row
	{
		auto tMat1 = rnd<cMat>(M,N);
		auto tMat2 = rnd<fMat>(M,M);

		const size_t n1 = genRndST(0,N-1);
		auto j1 = tMat1.cBegin()+n1;

		const size_t m2 = genRndST(0,M-1);
		auto j2 = tMat2.crBegin()+m2;

		*j1 = *j2;
		auto i2=j2->cbegin();
		CPPUNIT_ASSERT(std::all_of(j1->begin(),j1->end(),[&i2](const CPX__ i){return i==*i2++;}));
	}
	
	// from complex row
	{
		auto tMat1 = rnd<cMat>(M,N);
		auto tMat2 = rnd<cMat>(M,M);

		const size_t n1 = genRndST(0,N-1);
		auto j1 = tMat1.cBegin()+n1;

		const size_t m2 = genRndST(0,M-1);
		auto j2 = tMat2.crBegin()+m2;

		*j1 = *j2;
		auto i2=j2->cbegin();
		CPPUNIT_ASSERT(std::all_of(j1->begin(),j1->end(),[&i2](const CPX__ i){return i==*i2++;}));
	}

	// from real column
	{
		auto tMat1 = rnd<cMat>(M,N);
		auto tMat2 = rnd<fMat>(M,N);

		const size_t n1 = genRndST(0,N-1);
		auto j1 = tMat1.cBegin()+n1;

		const size_t n2 = genRndST(0,N-1);
		auto j2 = tMat2.ccBegin()+n2;

		*j1 = *j2;
		auto i2=j2->cbegin();
		CPPUNIT_ASSERT(std::all_of(j1->begin(),j1->end(),[&i2](const CPX__ i){return i==*i2++;}));
	}
	
	// from complex column
	{
		auto tMat1 = rnd<cMat>(M,N);
		auto tMat2 = rnd<cMat>(M,N);

		const size_t n1 = genRndST(0,N-1);
		auto j1 = tMat1.cBegin()+n1;

		const size_t n2 = genRndST(0,N-1);
		auto j2 = tMat2.ccBegin()+n2;

		*j1 = *j2;
		auto i2=j2->cbegin();
		CPPUNIT_ASSERT(std::all_of(j1->begin(),j1->end(),[&i2](const CPX__ i){return i==*i2++;}));
	}
	
	// from real matrix
	{
		auto tMat1 = rnd<cMat>(M,N);
		auto tMat2 = rnd<fMat>(M,1);

		const size_t n1 = genRndST(0,N-1);
		auto j1 = tMat1.cBegin()+n1;

		*j1 = tMat2;;
		auto i2=tMat2.cbegin();
		CPPUNIT_ASSERT(std::all_of(j1->begin(),j1->end(),[&i2](const CPX__ i){return i==*i2++;}));
	}
	
	// from complex matrix
	{
		auto tMat1 = rnd<cMat>(M,N);
		auto tMat2 = rnd<cMat>(M,1);

		const size_t n1 = genRndST(0,N-1);
		auto j1 = tMat1.cBegin()+n1;

		*j1 = tMat2;;
		auto i2=tMat2.cbegin();
		CPPUNIT_ASSERT(std::all_of(j1->begin(),j1->end(),[&i2](const CPX__ i){return i==*i2++;}));
	}
}

template<>
void test_tCol_all<CPX__,RE__,CPX__>::test_matrix_arithmetic() {
	
	const size_t M1 = genRndST();
	const size_t N1 = genRndST();
	const size_t M2 = genRndST();
	const size_t N2 = genRndST();

	// with real row
	{
		auto tMat1 = rnd<cMat>(M1,N1);
		auto tMat2 = rnd<fMat>(M2,N2);

		const size_t n_ = genRndST(0,N1-1);
		const size_t m_ = genRndST(0,M2-1);

		auto jc = tMat1.ccBegin()+n_;
		auto jr = tMat2.crBegin()+m_;
		
		const auto tMat3 = jc->prod(*jr);
		CPPUNIT_ASSERT(tMat3.cpx());
		CPPUNIT_ASSERT_EQUAL(M1,tMat3.M());
		CPPUNIT_ASSERT_EQUAL(N2,tMat3.N());

		for (size_t n=0; n!=N2; ++n)
			for (size_t m=0; m!=M1; ++m)
				CPPUNIT_ASSERT_DELTA((*jc)[m]*(*jr)[n],tMat3(m,n),delta);
	}
	
	// with complex row
	{
		auto tMat1 = rnd<cMat>(M1,N1);
		auto tMat2 = rnd<cMat>(M2,N2);

		const size_t n_ = genRndST(0,N1-1);
		const size_t m_ = genRndST(0,M2-1);

		auto jc = tMat1.ccBegin()+n_;
		auto jr = tMat2.crBegin()+m_;
		
		const auto tMat3 = jc->prod(*jr);
		CPPUNIT_ASSERT(tMat3.cpx());
		CPPUNIT_ASSERT_EQUAL(M1,tMat3.M());
		CPPUNIT_ASSERT_EQUAL(N2,tMat3.N());

		for (size_t n=0; n!=N2; ++n)
			for (size_t m=0; m!=M1; ++m)
				CPPUNIT_ASSERT_DELTA((*jc)[m]*(*jr)[n],tMat3(m,n),delta);
	}

	// with real matrix
	{
		auto tMat1 = rnd<cMat>(M1,N1);
		auto tMat2 = rnd<fMat>(1,N2);

		const size_t n_ = genRndST(0,N1-1);

		auto j = tMat1.ccBegin()+n_;
		
		const auto tMat3 = j->prod(tMat2);
		CPPUNIT_ASSERT(tMat3.cpx());
		CPPUNIT_ASSERT_EQUAL(M1,tMat3.M());
		CPPUNIT_ASSERT_EQUAL(N2,tMat3.N());

		for (size_t n=0; n!=N2; ++n)
			for (size_t m=0; m!=M1; ++m)
				CPPUNIT_ASSERT_DELTA((*j)[m]*tMat2[n],tMat3(m,n),delta);
	}
	
	// with complex matrix
	{
		auto tMat1 = rnd<cMat>(M1,N1);
		auto tMat2 = rnd<cMat>(1,N2);

		const size_t n_ = genRndST(0,N1-1);

		auto j = tMat1.ccBegin()+n_;
		
		const auto tMat3 = j->prod(tMat2);
		CPPUNIT_ASSERT(tMat3.cpx());
		CPPUNIT_ASSERT_EQUAL(M1,tMat3.M());
		CPPUNIT_ASSERT_EQUAL(N2,tMat3.N());

		for (size_t n=0; n!=N2; ++n)
			for (size_t m=0; m!=M1; ++m)
				CPPUNIT_ASSERT_DELTA((*j)[m]*tMat2[n],tMat3(m,n),delta);
	}
}

// test id
template<>
const char* test_tCol_all<CPX__,RE__,CPX__>::test_id() noexcept {
	return "test_cCol_all";
}

// instantiation
template class test_tCol_all<CPX__,RE__,CPX__>;
