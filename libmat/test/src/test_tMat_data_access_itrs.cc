// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_tMat_data_access_itrs.h"
#include "testTools.h"
#include <iostream>

using namespace lm__::test;

template<class TT, class FT, class CT>
void test_tMat_data_access_itrs<TT,FT,CT>::test_data_access() {

	const size_t M = genRndST();
	const size_t N = genRndST();

	// const
	{
		const auto tMat1 = rnd<tMat>(M,N);
		const auto d = tMat1.data();

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(i,size_t(&tMat1[i]-d));

		size_t i=0;
		for (size_t n=0; n!=N; ++n)
			for (size_t m=0; m!=M; ++m,++i)
				CPPUNIT_ASSERT_EQUAL(i,size_t(&tMat1(m,n)-d));

		CPPUNIT_ASSERT_EQUAL(d,&tMat1.front());
		CPPUNIT_ASSERT_EQUAL(d+tMat1.L()-1,&tMat1.back());
	}
	
	// non const
	{
		auto tMat1 = rnd<tMat>(M,N);
		auto d = tMat1.data();

		for (size_t i=0; i!=tMat1.L(); ++i)
			CPPUNIT_ASSERT_EQUAL(size_t(&tMat1[i]-d),i);

		size_t i=0;
		for (size_t n=0; n!=N; ++n)
			for (size_t m=0; m!=M; ++m,++i)
				CPPUNIT_ASSERT_EQUAL(size_t(&tMat1(m,n)-d),i);
		
		CPPUNIT_ASSERT_EQUAL(d,&tMat1.front());
		CPPUNIT_ASSERT_EQUAL(d+tMat1.L()-1,&tMat1.back());
	}
}

template<class TT, class FT, class CT>
void test_tMat_data_access_itrs<TT,FT,CT>::test_el_itr() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// const
	{
		const auto tMat1 = rnd<tMat>(M,N);

		{
			const auto tMat1 = rnd<tMat>(0,0);
			auto i = tMat1.begin();
			auto ie = tMat1.end();
			CPPUNIT_ASSERT_EQUAL(i,ie);
		}
		{
			auto i = tMat1.begin();
			auto ie = tMat1.end();

			CPPUNIT_ASSERT_EQUAL(tMat1.front(),*i);
			CPPUNIT_ASSERT_EQUAL(tMat1.L(),size_t(ie-i));
		
			for (size_t j=0; j!=tMat1.L(); ++j,++i)
				CPPUNIT_ASSERT_EQUAL(tMat1[j],*i);
		}
		{
			auto i = tMat1.cbegin();
			auto ie = tMat1.cend();

			CPPUNIT_ASSERT_EQUAL(tMat1.front(),*i);
			CPPUNIT_ASSERT_EQUAL(tMat1.L(),size_t(ie-i));
			
			for (size_t j=0; j!=tMat1.L(); ++j,++i)
				CPPUNIT_ASSERT_EQUAL(tMat1[j],*i);
		}
		{
			auto i = tMat1.rbegin();
			auto ie = tMat1.rend();

			CPPUNIT_ASSERT_EQUAL(tMat1.back(),*i);
			CPPUNIT_ASSERT_EQUAL(tMat1.L(),size_t(ie-i));
		
			for (size_t j=0; j!=tMat1.L(); ++j,++i)
				CPPUNIT_ASSERT_EQUAL(tMat1[tMat1.L()-1-j],*i);
		}
		{
			auto i = tMat1.crbegin();
			auto ie = tMat1.crend();

			CPPUNIT_ASSERT_EQUAL(tMat1.back(),*i);
			CPPUNIT_ASSERT_EQUAL(tMat1.L(),size_t(ie-i));
		
			for (size_t j=0; j!=tMat1.L(); ++j,++i)
				CPPUNIT_ASSERT_EQUAL(tMat1[tMat1.L()-1-j],*i);
		}
	}
	
	// non const
	{
		auto tMat1 = rnd<tMat>(M,N);
		
		{
			auto tMat1 = rnd<tMat>(0,0);
			auto i = tMat1.begin();
			auto ie = tMat1.end();
			CPPUNIT_ASSERT_EQUAL(i,ie);
		}
		{
			auto i = tMat1.begin();
			auto ie = tMat1.end();

			CPPUNIT_ASSERT_EQUAL(tMat1.front(),*i);
			CPPUNIT_ASSERT_EQUAL(tMat1.L(),size_t(ie-i));
		
			for (size_t j=0; j!=tMat1.L(); ++j,++i)
				CPPUNIT_ASSERT_EQUAL(tMat1[j],*i);
		}
		{
			auto i = tMat1.rbegin();
			auto ie = tMat1.rend();

			CPPUNIT_ASSERT_EQUAL(tMat1.back(),*i);
			CPPUNIT_ASSERT_EQUAL(tMat1.L(),size_t(ie-i));
		
			for (size_t j=0; j!=tMat1.L(); ++j,++i)
				CPPUNIT_ASSERT_EQUAL(tMat1[tMat1.L()-1-j],*i);
		}
	}
}

template<class TT, class FT, class CT>
void test_tMat_data_access_itrs<TT,FT,CT>::test_diag_itr() {

	const size_t M = genRndST();
	const size_t N = genRndST();

	// const
	{
		const auto tMat1 = rnd<tMat>(M,N);
		
		{
			const auto tMat1 = rnd<tMat>(0,0);
			auto i = tMat1.dbegin();
			auto ie = tMat1.dend();
			CPPUNIT_ASSERT_EQUAL(i,ie);
		}
		{
			// along m
			{
				{
					const size_t m_ = genRndST(0,M-1);

					auto i = tMat1.dbegin(m_,true);
					auto ie = tMat1.dend(m_,true);

					size_t m=m_, n=0;
					while (i!=ie)
						CPPUNIT_ASSERT_EQUAL(tMat1(m++,n++),*i++);
					CPPUNIT_ASSERT(m==M || n==N);
				}
				{
					const size_t m_ = genRndST(0,M-1);

					auto i = tMat1.cdbegin(m_,true);
					auto ie = tMat1.cdend(m_,true);

					size_t m=m_, n=0;
					while (i!=ie)
						CPPUNIT_ASSERT_EQUAL(tMat1(m++,n++),*i++);
					CPPUNIT_ASSERT(m==M || n==N);
				}
				{
					const size_t m_ = genRndST(0,M-1);

					auto i = tMat1.dbegin(m_,true);
					auto ie = tMat1.dend(m_,true);
					auto ri = tMat1.rdbegin(m_,true);
					auto rie = tMat1.rdend(m_,true);

					CPPUNIT_ASSERT_EQUAL((ie-i),(rie-ri));
					CPPUNIT_ASSERT_EQUAL(ie,ri.base());
					CPPUNIT_ASSERT_EQUAL(i,rie.base());
				}
				{
					const size_t m_ = genRndST(0,M-1);

					auto i = tMat1.cdbegin(m_,true);
					auto ie = tMat1.cdend(m_,true);
					auto ri = tMat1.crdbegin(m_,true);
					auto rie = tMat1.crdend(m_,true);

					CPPUNIT_ASSERT_EQUAL((ie-i),(rie-ri));
					CPPUNIT_ASSERT_EQUAL(ie,ri.base());
					CPPUNIT_ASSERT_EQUAL(i,rie.base());
				}
			}

			// along n
			{
				{
					const size_t n_ = genRndST(0,N-1);

					auto i = tMat1.dbegin(n_,false);
					auto ie = tMat1.dend(n_,false);

					size_t m=0, n=n_;
					while (i!=ie)
						CPPUNIT_ASSERT_EQUAL(tMat1(m++,n++),*i++);
					CPPUNIT_ASSERT(m==M || n==N);
				}
				{
					const size_t n_ = genRndST(0,N-1);

					auto i = tMat1.cdbegin(n_,false);
					auto ie = tMat1.cdend(n_,false);

					size_t m=0, n=n_;
					while (i!=ie)
						CPPUNIT_ASSERT_EQUAL(tMat1(m++,n++),*i++);
					CPPUNIT_ASSERT(m==M || n==N);
				}
				{
					const size_t n_ = genRndST(0,N-1);

					auto i = tMat1.dbegin(n_,false);
					auto ie = tMat1.dend(n_,false);
					auto ri = tMat1.rdbegin(n_,false);
					auto rie = tMat1.rdend(n_,false);

					CPPUNIT_ASSERT_EQUAL((ie-i),(rie-ri));
					CPPUNIT_ASSERT_EQUAL(ie,ri.base());
					CPPUNIT_ASSERT_EQUAL(i,rie.base());
				}
				{
					const size_t n_ = genRndST(0,N-1);

					auto i = tMat1.cdbegin(n_,false);
					auto ie = tMat1.cdend(n_,false);
					auto ri = tMat1.crdbegin(n_,false);
					auto rie = tMat1.crdend(n_,false);

					CPPUNIT_ASSERT_EQUAL((ie-i),(rie-ri));
					CPPUNIT_ASSERT_EQUAL(ie,ri.base());
					CPPUNIT_ASSERT_EQUAL(i,rie.base());
				}
			}
		}
	}

	// non const
	{
		auto tMat1 = rnd<tMat>(M,N);
		
		{
			auto tMat1 = rnd<tMat>(0,0);
			auto i = tMat1.dbegin();
			auto ie = tMat1.dend();
			CPPUNIT_ASSERT_EQUAL(i,ie);
		}
		{
			// along m
			{
				{
					const size_t m_ = genRndST(0,M-1);

					auto i = tMat1.dbegin(m_,true);
					auto ie = tMat1.dend(m_,true);

					size_t m=m_, n=0;
					while (i!=ie)
						CPPUNIT_ASSERT_EQUAL(tMat1(m++,n++),*i++);
					CPPUNIT_ASSERT(m==M || n==N);
				}
				{
					const size_t m_ = genRndST(0,M-1);

					auto i = tMat1.cdbegin(m_,true);
					auto ie = tMat1.cdend(m_,true);

					size_t m=m_, n=0;
					while (i!=ie)
						CPPUNIT_ASSERT_EQUAL(tMat1(m++,n++),*i++);
					CPPUNIT_ASSERT(m==M || n==N);
				}
				{
					const size_t m_ = genRndST(0,M-1);

					auto i = tMat1.dbegin(m_,true);
					auto ie = tMat1.dend(m_,true);
					auto ri = tMat1.rdbegin(m_,true);
					auto rie = tMat1.rdend(m_,true);

					CPPUNIT_ASSERT_EQUAL((ie-i),(rie-ri));
					CPPUNIT_ASSERT_EQUAL(ie,ri.base());
					CPPUNIT_ASSERT_EQUAL(i,rie.base());
				}
				{
					const size_t m_ = genRndST(0,M-1);

					auto i = tMat1.cdbegin(m_,true);
					auto ie = tMat1.cdend(m_,true);
					auto ri = tMat1.crdbegin(m_,true);
					auto rie = tMat1.crdend(m_,true);

					CPPUNIT_ASSERT_EQUAL((ie-i),(rie-ri));
					CPPUNIT_ASSERT_EQUAL(ie,ri.base());
					CPPUNIT_ASSERT_EQUAL(i,rie.base());
				}
			}

			// along n
			{
				{
					const size_t n_ = genRndST(0,N-1);

					auto i = tMat1.dbegin(n_,false);
					auto ie = tMat1.dend(n_,false);

					size_t m=0, n=n_;
					while (i!=ie)
						CPPUNIT_ASSERT_EQUAL(tMat1(m++,n++),*i++);
					CPPUNIT_ASSERT(m==M || n==N);
				}
				{
					const size_t n_ = genRndST(0,N-1);

					auto i = tMat1.cdbegin(n_,false);
					auto ie = tMat1.cdend(n_,false);

					size_t m=0, n=n_;
					while (i!=ie)
						CPPUNIT_ASSERT_EQUAL(tMat1(m++,n++),*i++);
					CPPUNIT_ASSERT(m==M || n==N);
				}
				{
					const size_t n_ = genRndST(0,N-1);

					auto i = tMat1.dbegin(n_,false);
					auto ie = tMat1.dend(n_,false);
					auto ri = tMat1.rdbegin(n_,false);
					auto rie = tMat1.rdend(n_,false);

					CPPUNIT_ASSERT_EQUAL((ie-i),(rie-ri));
					CPPUNIT_ASSERT_EQUAL(ie,ri.base());
					CPPUNIT_ASSERT_EQUAL(i,rie.base());
				}
				{
					const size_t n_ = genRndST(0,N-1);

					auto i = tMat1.cdbegin(n_,false);
					auto ie = tMat1.cdend(n_,false);
					auto ri = tMat1.crdbegin(n_,false);
					auto rie = tMat1.crdend(n_,false);

					CPPUNIT_ASSERT_EQUAL((ie-i),(rie-ri));
					CPPUNIT_ASSERT_EQUAL(ie,ri.base());
					CPPUNIT_ASSERT_EQUAL(i,rie.base());
				}
			}
		}
	}
}

template<class TT, class FT, class CT>
void test_tMat_data_access_itrs<TT,FT,CT>::test_row_col_access() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// const
	{
		const auto tMat1 = rnd<tMat>(M,N);
		
		// rAt, rFront, rBack
		{
			const size_t m = genRndST(0,M-1);
			auto j = tMat1.rBegin()+m;

			CPPUNIT_ASSERT(*j==tMat1.rAt(m));
			CPPUNIT_ASSERT(*tMat1.rBegin()==tMat1.rFront());
			CPPUNIT_ASSERT(*(tMat1.rBegin()+tMat1.M()-1)==tMat1.rBack());
		}
		
		// cAt, cFront, cBack
		{
			const size_t n = genRndST(0,N-1);
			auto j = tMat1.cBegin()+n;

			CPPUNIT_ASSERT(*j==tMat1.cAt(n));
			CPPUNIT_ASSERT(*tMat1.cBegin()==tMat1.cFront());
			CPPUNIT_ASSERT(*(tMat1.cBegin()+tMat1.N()-1)==tMat1.cBack());
		}
	}
	
	// non const
	{
		auto tMat1 = rnd<tMat>(M,N);
		
		// rAt, rFront, rBack
		{
			const size_t m = genRndST(0,M-1);
			auto j = tMat1.rBegin()+m;

			CPPUNIT_ASSERT(*j==tMat1.rAt(m));
			CPPUNIT_ASSERT(*tMat1.rBegin()==tMat1.rFront());
			CPPUNIT_ASSERT(*(tMat1.rBegin()+tMat1.M()-1)==tMat1.rBack());
		}
		
		// cAt, cFront, cBack
		{
			const size_t n = genRndST(0,N-1);
			auto j = tMat1.cBegin()+n;

			CPPUNIT_ASSERT(*j==tMat1.cAt(n));
			CPPUNIT_ASSERT(*tMat1.cBegin()==tMat1.cFront());
			CPPUNIT_ASSERT(*(tMat1.cBegin()+tMat1.N()-1)==tMat1.cBack());
		}
	}
}

template<class TT, class FT, class CT>
void test_tMat_data_access_itrs<TT,FT,CT>::test_row_itr() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// const
	{
		const auto tMat1 = rnd<tMat>(M,N);
		
		{
			const auto tMat1 = rnd<tMat>(0,0);
			auto i = tMat1.rBegin();
			auto ie = tMat1.rEnd();
			CPPUNIT_ASSERT_EQUAL(i,ie);
		}
		{
			auto i = tMat1.rBegin();
			auto ie = tMat1.rEnd();

			CPPUNIT_ASSERT_EQUAL(tMat1.M(),size_t(ie-i));
			for (size_t j=0; j!=tMat1.M(); ++j,++i)
				CPPUNIT_ASSERT(tMat1.rAt(j)==*i);
		}
		{
			auto i = tMat1.crBegin();
			auto ie = tMat1.crEnd();

			CPPUNIT_ASSERT_EQUAL(tMat1.M(),size_t(ie-i));
			for (size_t j=0; j!=tMat1.M(); ++j,++i)
				CPPUNIT_ASSERT(tMat1.rAt(j)==*i);
		}
	}
	
	// non const
	{
		auto tMat1 = rnd<tMat>(M,N);
		{
			auto tMat1 = rnd<tMat>(0,0);
			auto i = tMat1.rBegin();
			auto ie = tMat1.rEnd();
			CPPUNIT_ASSERT_EQUAL(i,ie);
		}
		{
			auto i = tMat1.rBegin();
			auto ie = tMat1.rEnd();

			CPPUNIT_ASSERT_EQUAL(tMat1.M(),size_t(ie-i));
			for (size_t j=0; j!=tMat1.M(); ++j,++i)
				CPPUNIT_ASSERT(tMat1.rAt(j)==*i);
		}
	}
}

template<class TT, class FT, class CT>
void test_tMat_data_access_itrs<TT,FT,CT>::test_col_itr() {
	
	const size_t M = genRndST();
	const size_t N = genRndST();

	// const
	{
		const auto tMat1 = rnd<tMat>(M,N);
		
		{
			const auto tMat1 = rnd<tMat>(0,0);
			auto i = tMat1.cBegin();
			auto ie = tMat1.cEnd();
			CPPUNIT_ASSERT_EQUAL(i,ie);
		}
		{
			auto i = tMat1.cBegin();
			auto ie = tMat1.cEnd();

			CPPUNIT_ASSERT_EQUAL(tMat1.N(),size_t(ie-i));
			for (size_t j=0; j!=tMat1.N(); ++j,++i)
				CPPUNIT_ASSERT(tMat1.cAt(j)==*i);
		}
		{
			auto i = tMat1.ccBegin();
			auto ie = tMat1.ccEnd();

			CPPUNIT_ASSERT_EQUAL(tMat1.N(),size_t(ie-i));
			for (size_t j=0; j!=tMat1.N(); ++j,++i)
				CPPUNIT_ASSERT(tMat1.cAt(j)==*i);
		}
	}
	
	// non const
	{
		auto tMat1 = rnd<tMat>(M,N);
		
		{
			auto tMat1 = rnd<tMat>(0,0);
			auto i = tMat1.cBegin();
			auto ie = tMat1.cEnd();
			CPPUNIT_ASSERT_EQUAL(i,ie);
		}
		{
			auto i = tMat1.cBegin();
			auto ie = tMat1.cEnd();

			CPPUNIT_ASSERT_EQUAL(tMat1.N(),size_t(ie-i));
			for (size_t j=0; j!=tMat1.N(); ++j,++i)
				CPPUNIT_ASSERT(tMat1.cAt(j)==*i);
		}
	}
}


template<class TT, class FT, class CT>
CppUnit::Test* test_tMat_data_access_itrs<TT,FT,CT>::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_tMat_data_access_itrs>(
		"test_data_access", &test_tMat_data_access_itrs<TT,FT,CT>::test_data_access));
	suite->addTest(new CppUnit::TestCaller<test_tMat_data_access_itrs>(
		"test_el_itr", &test_tMat_data_access_itrs<TT,FT,CT>::test_el_itr));
	suite->addTest(new CppUnit::TestCaller<test_tMat_data_access_itrs>(
		"test_diag_itr", &test_tMat_data_access_itrs<TT,FT,CT>::test_diag_itr));
	suite->addTest(new CppUnit::TestCaller<test_tMat_data_access_itrs>(
		"test_row_col_access", &test_tMat_data_access_itrs<TT,FT,CT>::test_row_col_access));
	suite->addTest(new CppUnit::TestCaller<test_tMat_data_access_itrs>(
		"test_row_itr", &test_tMat_data_access_itrs<TT,FT,CT>::test_row_itr));
	suite->addTest(new CppUnit::TestCaller<test_tMat_data_access_itrs>(
		"test_col_itr", &test_tMat_data_access_itrs<TT,FT,CT>::test_col_itr));
	
	return suite;
}
