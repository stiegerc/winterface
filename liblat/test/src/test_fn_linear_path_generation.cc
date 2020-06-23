// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_fn_linear_path_generation.h"
#include "ll_fn.h"
#include "ll_cell.h"
#include "ll_testTools.h"
#include "lm_testTools.h"


using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;


void test_fn_linear_path_generation::test_genPath() {
	const size_t smax=100;
	
	const auto S = round(rand<fMat>(1,smax,1.501,20.499));
	const auto D = round(rand<fMat>(1,smax,0.501,10.499));
	
	// explicit number of points version
	{
		// check inserting empty P
		CPPUNIT_ASSERT(genPath(fMat(3,0),std::vector<size_t>(2)).empty());
		
		// check inserting P of size dx1
		{
			fMat P(3,1);
			CPPUNIT_ASSERT(genPath(P,std::vector<size_t>(2)).path()==P);
			CPPUNIT_ASSERT(genPath(P,std::vector<size_t>(2)).pos()==0.0);
		}
		
		// check if exceptions are thrown
		{
			CPPUNIT_ASSERT_THROW(genPath(fMat(3,5),{3,4,3}),std::invalid_argument);
			CPPUNIT_ASSERT_THROW(genPath(fMat(3,5),{3,4,3,1}),std::invalid_argument);
		}

		// check randomly generated P, Nps
		{
			for (size_t i=0; i!=smax; ++i) {
				const size_t s = S[i];
				const size_t d = D[i];
				
				const auto P = rand<fMat>(d,s,-10.0,10.0);
				const auto v = round(rand<fMat>(1,s-1,1.501,200.499));
				std::vector<size_t> Nps; Nps.reserve(s-1);
				for (auto j: v) Nps.push_back(size_t(j));
				
				const auto path = genPath(P,Nps);
				const size_t Np = std::accumulate(Nps.begin(),Nps.end(),0);
				
				CPPUNIT_ASSERT_EQUAL(d,path.dim());
				CPPUNIT_ASSERT_EQUAL(Np,path.N());
				
				// check start and ending points
				{
					size_t sind=0;
					for (size_t j=0; j!=Nps.size(); ++j) {
						CPPUNIT_ASSERT(path.cAt(sind)==P.cAt(j));
						CPPUNIT_ASSERT(path.cAt(sind+Nps[j]-1)==P.cAt(j+1));
						sind+=Nps[j];
					}
				}
				
				// check that points are on a straight line
				{
					size_t sind=0;
					for (size_t j=0; j!=Nps.size(); ++j) {
						const auto v = path.cAt(sind+1) - 
							       path.cAt(sind);
						for (size_t ind=sind; ind!=sind+Nps[j]-1; ++ind)
							CPPUNIT_ASSERT(v==path.cAt(ind+1)-path.cAt(ind));
						sind+=Nps[j];
					}
				}
			}
		}
	}
	
	// total number of points, equal spacing version
	{
		// check inserting empty P
		CPPUNIT_ASSERT(genPath(fMat(3,0),2).empty());
		
		// check inserting P of size dx1
		{
			const fMat P(3,1);
			CPPUNIT_ASSERT(genPath(P,2).path()==P);
			CPPUNIT_ASSERT(genPath(P,2).pos()==0.0);
		}
		
		// check if exceptions are thrown
		{
			const fMat P({1.0,0.0,1.1,0.0,1.0,1.0},2,3);
			CPPUNIT_ASSERT_THROW(genPath(P,3),std::invalid_argument);
		}

		// check randomly generated P, Np
		{
			for (size_t i=0; i!=smax; ++i) {
				const size_t s = S[i];
				const size_t d = D[i];
				
				fMat P(d,s);
				P.cFront() = rand<fMat>(d,1,-2.0,2.0);
				size_t j=1;
				while (j!=P.N()) {
					const auto t = rand<fMat>(d,1,-2.0,2.0);
					const auto n = norm(t-P.cAt(j-1));
					if (n>1.0) P.cAt(j++) = t;
				}

				const size_t Np = rand<fMat>(1,1,double(10*s)-.499,double(30*s)+0.499)[0];
				const auto path = genPath(P,Np);
				
				CPPUNIT_ASSERT_EQUAL(d,path.dim());
				CPPUNIT_ASSERT_EQUAL(Np,path.N());
				
				// check if all points of P are included
				std::vector<size_t> inds(P.N());
				size_t cind=0;;
				for (size_t j=0; j!=P.N(); ++j) {
					inds[j] = std::find(path.ccBegin()+cind,
							path.ccEnd(),P.cAt(j)).i();
					CPPUNIT_ASSERT(inds[j]<path.N());
					cind=inds[j];
				}
				
				// check if the extracted inds are in ascending order
				size_t cur;
				for (size_t j=0; j!=inds.size()-1; ++j) {
					cur = inds[j];
					CPPUNIT_ASSERT(cur<inds[j+1]);
				}
				
				// create Nps vector for other version of Nps, then compare paths
				std::vector<size_t> Nps(inds.size()-1);
				for (size_t j=0; j!=Nps.size(); ++j)
					Nps[j] = inds[j+1]-inds[j];
				++Nps[0];
				CPPUNIT_ASSERT_EQUAL(path.N(),size_t(std::accumulate(Nps.begin(),Nps.end(),0)));
				
				const auto path_ = genPath(P,Nps);
				CPPUNIT_ASSERT(path.path()==path_.path());
				
				// check if entries in Nps correspond to distances in P
				double l=0.0;
				for (size_t j=0; j!=Nps.size(); ++j)
					l+=norm(P.cAt(j+1)-P.cAt(j));
				for (size_t j=0; j!=Nps.size(); ++j)
					CPPUNIT_ASSERT(std::abs(norm(P.cAt(j+1)
					-P.cAt(j))/l-double(Nps[j])/double(Np))<1e-1);
					// allow 10% missmatch in comparing ratios of steps
					// and lengths, gets better for more steps
			}
		}
	}
}


const char* test_fn_linear_path_generation::test_id() noexcept {
	return "test_fn_linear_path_generation";
}

CppUnit::Test* test_fn_linear_path_generation::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_fn_linear_path_generation>(
		"test_genPath", &test_fn_linear_path_generation::test_genPath));
	
	return suite;
}
