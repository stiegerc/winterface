// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_fn_metrics_and_clustering.h"
#include "ll_fn.h"
#include "ll_cell.h"
#include "ll_testTools.h"
#include "lm_testTools.h"


using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;


void test_fn_metrics_and_clustering::test_genNNmat() {

	// check no restriction
	{
		const fMat NN({-1.0,-1.0,-1.0,
			       -1.0,-1.0, 0.0,
			       -1.0,-1.0, 1.0,
			       -1.0, 0.0,-1.0,
			       -1.0, 0.0, 0.0,
			       -1.0, 0.0, 1.0,
			       -1.0, 1.0,-1.0,
			       -1.0, 1.0, 0.0,
			       -1.0, 1.0, 1.0,
				0.0,-1.0,-1.0,
				0.0,-1.0, 0.0,
				0.0,-1.0, 1.0,
				0.0, 0.0,-1.0,
				0.0, 0.0, 0.0,
				0.0, 0.0, 1.0,
				0.0, 1.0,-1.0,
				0.0, 1.0, 0.0,
				0.0, 1.0, 1.0,
				1.0,-1.0,-1.0,
				1.0,-1.0, 0.0,
				1.0,-1.0, 1.0,
				1.0, 0.0,-1.0,
				1.0, 0.0, 0.0,
				1.0, 0.0, 1.0,
				1.0, 1.0,-1.0,
				1.0, 1.0, 0.0,
				1.0, 1.0, 1.0},DIM__);

		const auto nn = genNNmat({false,false,false});
		CPPUNIT_ASSERT(nn==NN);
	}
	
	// check default to be no restrictions
	{
		CPPUNIT_ASSERT(genNNmat()==genNNmat({false,false,false}));
	}
	
	// check y axis restricted
	{
		const fMat NN({-1.0, 0.0,-1.0,
			       -1.0, 0.0, 0.0,
			       -1.0, 0.0, 1.0,
				0.0, 0.0,-1.0,
				0.0, 0.0, 0.0,
				0.0, 0.0, 1.0,
				1.0, 0.0,-1.0,
				1.0, 0.0, 0.0,
				1.0, 0.0, 1.0},DIM__);
		
		const auto nn = genNNmat({false,true,false});
		CPPUNIT_ASSERT(nn==NN);
	}
	
	// check z axis restricted
	{
		const fMat NN({-1.0,-1.0, 0.0,
			       -1.0, 0.0, 0.0,
			       -1.0, 1.0, 0.0,
				0.0,-1.0, 0.0,
				0.0, 0.0, 0.0,
				0.0, 1.0, 0.0,
				1.0,-1.0, 0.0,
				1.0, 0.0, 0.0,
				1.0, 1.0, 0.0},DIM__);
		
		const auto nn = genNNmat({false,false,true});
		CPPUNIT_ASSERT(nn==NN);
	}

	// check y and z axis restricted
	{
		const fMat NN({-1.0, 0.0, 0.0,
				0.0, 0.0, 0.0,
				1.0, 0.0, 0.0,},DIM__);
		
		const auto nn = genNNmat({false,true,true});
		CPPUNIT_ASSERT(nn==NN);
	}
	
	// check x,y,z axis restricted
	{
		const auto nn = genNNmat({true,true,true});
		CPPUNIT_ASSERT(nn==zeros<fMat>(DIM__,1));
	}
}

void test_fn_metrics_and_clustering::test_dist_distb() {
	
	// check dist, no f
        {		
		const auto NN = genNNmat();

		const auto B = ones<fMat>(3)-eye<fMat>(3);
		{
			const auto d = dist(B,fMat({0.0,0.0,0.0}),fMat({0.5,0.0,0.0}),NN);
			CPPUNIT_ASSERT_DELTA(norm(B.cAt(0))/2.0,d,mtol());
                }
                {
			const auto d = dist(B,fMat({0.0,0.0,0.0}),fMat({0.25,0.0,0.0}),NN);
			CPPUNIT_ASSERT_DELTA(norm(B.cAt(0))/4.0,d,mtol());
                }
                {
			const auto d = dist(B,fMat({0.0,0.0,0.0}),fMat({0.75,0.0,0.0}),NN);
			CPPUNIT_ASSERT_DELTA(norm(B.cAt(0))/4.0,d,mtol());
                }
	}

	// check dist
	{
		const auto NN = genNNmat();

		// check using zincblende
		{
			const auto B = .5*(ones<fMat>(3)-eye<fMat>(3));
			
			fMat Ap(3,0); Ap.reserve(2);
			Ap.push_back(zeros<fMat>(3,1));
			Ap.push_back(rand<fMat>(3,1,.23,.27));
			
			// check dynamic f
			const double fdyn = -1.1;

			const auto tb01 = dist(B,Ap.cAt(0),Ap.cAt(1),fdyn,NN);
			const auto tb10 = dist(B,Ap.cAt(1),Ap.cAt(0),fdyn,NN);

			const auto mb = norm(B.prod(ones<fMat>(3,1)*.27001));
			CPPUNIT_ASSERT(tb01>0.0 && tb01<mb);
			CPPUNIT_ASSERT_DELTA(tb01,tb10,mtol());
			{
				const auto tb01_o = dist(B,Ap.cAt(0),Ap.cAt(1),-1.0,NN);
				const auto tb01_s = dist(B,Ap.cAt(0),Ap.cAt(1),NN);
				CPPUNIT_ASSERT_DELTA(tb01_o,tb01_s,mtol());
			}

			// check hard cutoff f
			{
				const auto tb01_hc = dist(B,Ap.cAt(0),Ap.cAt(1),mb,NN);
				CPPUNIT_ASSERT(tb01_hc<=mb);
			}

			// check without f versus f==0.0
			{
				const auto tb01_nf = dist(B,Ap.cAt(0),Ap.cAt(1),NN);
				const auto tb01_z = dist(B,Ap.cAt(0),Ap.cAt(1),0.0,NN);
				CPPUNIT_ASSERT_DELTA(tb01_nf,tb01_z,mtol());
			}
		}
		
	}

	// check distb, no f
	{
		const auto NN = genNNmat();
		
		// check using zincblende
		{
			const auto B = .5*(ones<fMat>(3)-eye<fMat>(3));
			const fMat Ap({0.0,0.0,0.0,.25,.25,.25},3);

			const auto tb01 = distb(B,Ap.cAt(0),Ap.cAt(1),NN);
			const auto tb10 = distb(B,Ap.cAt(1),Ap.cAt(0),NN);

			CPPUNIT_ASSERT_DELTA(tb01.val,sqrt(3.0)/4.0,mtol());
			CPPUNIT_ASSERT_DELTA(tb01.val,tb10.val,mtol());
			CPPUNIT_ASSERT_EQUAL(tb01.bonds.N(),tb10.bonds.N());
			CPPUNIT_ASSERT_EQUAL(size_t(0),tb01.bonds.N()%2);

			for (auto i=tb01.bonds.ccBegin(),e=tb01.bonds.ccEnd(); i!=e; ++i)
				CPPUNIT_ASSERT(std::find(tb10.bonds.ccBegin(),tb10.bonds.ccEnd(),-*i).i()
					!= ptrdiff_t(tb01.bonds.N()));
		}
	}
	
	// check distb
	{
		const auto NN = genNNmat();

		// check using slightly randomized zincblende
		{
			const auto B = .5*(ones<fMat>(3)-eye<fMat>(3));
			
			fMat Ap(3,0); Ap.reserve(2);
			Ap.push_back(zeros<fMat>(3,1));
			Ap.push_back(rand<fMat>(3,1,.23,.27));
			
			
			// check dynamic f
			const double fdyn = -1.1;

			const auto tb01 = distb(B,Ap.cAt(0),Ap.cAt(1),fdyn,NN);
			const auto tb10 = distb(B,Ap.cAt(1),Ap.cAt(0),fdyn,NN);

			const auto mb = norm(B.prod(ones<fMat>(3,1)*.27001));
			CPPUNIT_ASSERT(tb01.val>0.0 && tb01.val<mb);
			CPPUNIT_ASSERT_DELTA(tb01.val,tb10.val,mtol());
			CPPUNIT_ASSERT_EQUAL(tb01.bonds.N(),tb10.bonds.N());
			
			for (auto i=tb01.bonds.ccBegin(),e=tb01.bonds.ccEnd(); i!=e; ++i)
				CPPUNIT_ASSERT(std::find(tb10.bonds.ccBegin(),tb10.bonds.ccEnd(),-*i).i()
						!= ptrdiff_t(tb01.bonds.N()));
			
		
			auto n = mnorm(B.prod(tb01.bonds)); std::sort(n.begin(),n.end());
			for (size_t i=1; i<n.N(); ++i)
				CPPUNIT_ASSERT(-n[i-1]*fdyn > n[i]);
			{
				auto tb01_o = distb(B,Ap.cAt(0),Ap.cAt(1),-1.0,NN);
				auto tb01_s = distb(B,Ap.cAt(0),Ap.cAt(1),n[0],NN);
				CPPUNIT_ASSERT_DELTA(tb01_o.val,tb01_s.val,mtol());
				CPPUNIT_ASSERT(tb01_o.bonds==tb01_s.bonds);
			}

			
			// check hard cutoff f
			for (size_t i=0; i<n.N(); ++i) {
				const double fhc = n[i];
				const auto tb01_hc = distb(B,Ap.cAt(0),Ap.cAt(1),fhc,NN);
				CPPUNIT_ASSERT(all(mnorm(B.prod(tb01_hc.bonds)).leq(fhc)));
			}
			{
				const auto tb01_z = distb(B,Ap.cAt(0),Ap.cAt(1),0.0,NN);
				const auto tb01_s = distb(B,Ap.cAt(0),Ap.cAt(1),n[0],NN);
				CPPUNIT_ASSERT_DELTA(tb01_z.val,tb01_s.val,mtol());
			}
		}
	}
}

void test_fn_metrics_and_clustering::test_genDmat() {

	// same positions
	{	
		const size_t D = DIM__;
		const auto tCell1 = genRandom(D);
		
		const auto r = genRndrv(D);
		const auto NN = genNNmat(r);
		
		// no f
		const auto Dmat = genDmat(tCell1.B(),tCell1.Ap(),NN);
		CPPUNIT_ASSERT(Dmat.square());
		CPPUNIT_ASSERT_EQUAL(tCell1.N(),Dmat.N());
		CPPUNIT_ASSERT(T(Dmat)==Dmat);
		CPPUNIT_ASSERT(std::all_of(Dmat.cdbegin(),Dmat.cdend(),
			[](const double i){return lm__::ops::z(i);}));
		
		// with f as direct cutoff
		{
			const double f = std::accumulate(Dmat.cbegin(),Dmat.cend(),0.0)/Dmat.L();
			const auto Dmat_p = genDmat(tCell1.B(),tCell1.Ap(),f,NN);
			CPPUNIT_ASSERT(size(Dmat)==size(Dmat_p));
			CPPUNIT_ASSERT(T(Dmat_p)==Dmat_p);
			CPPUNIT_ASSERT(std::all_of(Dmat_p.cdbegin(),Dmat_p.cdend(),
				[](const double i){return lm__::ops::z(i);}));
			
			CPPUNIT_ASSERT(all(Dmat.leq(Dmat_p)));
		}
		
		// with f as fractional
		{
			const double f = -1.1;
			const auto Dmat_m = genDmat(tCell1.B(),tCell1.Ap(),f,NN);
			CPPUNIT_ASSERT(size(Dmat)==size(Dmat_m));
			CPPUNIT_ASSERT(T(Dmat_m)==Dmat_m);
			CPPUNIT_ASSERT(std::all_of(Dmat_m.cdbegin(),Dmat_m.cdend(),
				[](const double i){return lm__::ops::z(i);}));
			
			CPPUNIT_ASSERT(all(Dmat.leq(Dmat_m)));
		}

		// primitive Zincblende, no f
		{
			const double a0 = genRndDouble(1.1,1.9);
			const auto tCell = genZincblende(D,a0);

			const auto D = genDmat(tCell.B(),tCell.Ap(),NN);
			
			CPPUNIT_ASSERT_EQUAL(size_t(2),D.M());
			CPPUNIT_ASSERT_EQUAL(size_t(2),D.N());

			const fMat Dref({0.0,sqrt(3.0)*a0/4.0,sqrt(3.0)*a0/4.0,0.0},2,2);
			CPPUNIT_ASSERT(Dref==D);
		}
	}

	// random different positions
	{
		const size_t M = 3;
		const size_t N1 = genRndST();
		const size_t N2 = genRndST();
		const double fdyn = -1.1;
		
		const auto B = rand<fMat>(M,M)%1.0;
		const auto P1 = rand<fMat>(M,N1)%1.0;
		const auto P2 = rand<fMat>(M,N2)%1.0;
		
		const auto r = genRndrv(M);
		const auto NN = genNNmat(r);

		// no f
		{
			const auto D = genDmat(B,P1,P2,NN);
			for (auto i=P1.ccBegin(),ie=P1.ccEnd(); i!=ie; ++i)
				for (auto j=P2.ccBegin(),je=P2.ccEnd(); j!=je; ++j)
					CPPUNIT_ASSERT_DELTA(D(i.i(),j.i()),dist(B,*i,*j,NN),mtol());
		}
		// with f
		{
			const auto D = genDmat(B,P1,P2,fdyn,NN);
			for (auto i=P1.ccBegin(),ie=P1.ccEnd(); i!=ie; ++i)
				for (auto j=P2.ccBegin(),je=P2.ccEnd(); j!=je; ++j)
					CPPUNIT_ASSERT_DELTA(D(i.i(),j.i()),dist(B,*i,*j,fdyn,NN),mtol());
		}
	}
}

void test_fn_metrics_and_clustering::test_com() {
	
	// check com
	{
		{
			const fMat Ap({0.1,0.1,0.1,
					0.9,0.9,0.9},3);
			CPPUNIT_ASSERT(com(Ap)==fMat({0.0,0.0,0.0}));
		}	
		{
			const fMat Ap({0.4,0.4,0.4,
					0.6,0.6,0.6},3);
			CPPUNIT_ASSERT(com(Ap)==fMat({0.5,0.5,0.5}));
		}	
		{
			const fMat Ap({0.1,0.6,0.7,
					0.9,0.4,0.1},3);
			CPPUNIT_ASSERT(com(Ap)==fMat({0.0,0.5,0.9}));
		}

		{
			const fMat Ap({0.1,0.6,0.7,
					0.9,0.4,0.1},3);
			const fMat w({1,2},1,2);
			CPPUNIT_ASSERT(com(Ap,w)==fMat({0.962183826553324,
							0.462183826553324,
							0.027034170688061}));
		}
	}
}

void test_fn_metrics_and_clustering::test_dbscan() {
	// check dbscan
	{
		// check exceptions
		{
			CPPUNIT_ASSERT_THROW(dbscan(fMat(2,3),2,0.1), std::invalid_argument);
			CPPUNIT_ASSERT_THROW(dbscan(-ones<fMat>(3),2,0.1), std::invalid_argument);
			CPPUNIT_ASSERT_THROW(dbscan(ones<fMat>(3),2,-0.1), std::invalid_argument);
		}

		// check simple case
		{
			const double s=.2;
			fMat Ap(3,0); Ap.reserve(21);
			Ap.push_back(rand<fMat>(3,9,1.0-s,1.0));
			Ap.push_back(rand<fMat>(3,11,0.0,s));
			Ap.push_back(fMat({.5,.5,.5}));
		
			fMat D(Ap.N());
			for (size_t i=0; i<D.M(); ++i)
				for (size_t j=0; j<D.N(); ++j)
					D(i,j) = norm(Ap.cAt(i)-Ap.cAt(j));

			auto C = dbscan(D,2,s);

			wi ck = {{0,1,2,3,4,5,6,7,8},
			       {9,10,11,12,13,14,15,16,17,18,19},
			       {20}};


			CPPUNIT_ASSERT(std::unique(C.begin(),C.end())-C.begin()==3);
			CPPUNIT_ASSERT(C[0].size()+C[1].size()+C[2].size()==21);

			for (const auto& i: C) {
				auto c = ck.begin();
				while(c<ck.end()) {
					if (i==*c) break;
					++c;
				}
				CPPUNIT_ASSERT(c!=ck.end());
			}
		}
	}
}


const char* test_fn_metrics_and_clustering::test_id() noexcept {
	return "test_fn_metrics_and_clustering";
}

CppUnit::Test* test_fn_metrics_and_clustering::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_fn_metrics_and_clustering>(
		"test_genNNmat", &test_fn_metrics_and_clustering::test_genNNmat));
	suite->addTest(new CppUnit::TestCaller<test_fn_metrics_and_clustering>(
		"test_dist_distb", &test_fn_metrics_and_clustering::test_dist_distb));
	suite->addTest(new CppUnit::TestCaller<test_fn_metrics_and_clustering>(
		"test_genDmat", &test_fn_metrics_and_clustering::test_genDmat));
	suite->addTest(new CppUnit::TestCaller<test_fn_metrics_and_clustering>(
		"test_com", &test_fn_metrics_and_clustering::test_com));
	suite->addTest(new CppUnit::TestCaller<test_fn_metrics_and_clustering>(
		"test_dbscan", &test_fn_metrics_and_clustering::test_dbscan));
	
	return suite;
}
