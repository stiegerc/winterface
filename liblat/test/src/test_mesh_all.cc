// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_mesh_all.h"
#include "lm_testTools.h"
#include "ll_testTools.h"
#include "aux_io.h"
#include "ll_io.h"
#include <iostream>


using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;
using namespace aux;


void test_mesh_all::test_itr() {

	const size_t M = genRndST(1,10);
	const size_t D1 = genRndST(2,10);
	const size_t D2 = genRndST(2,10);

	auto tMat = rand<fMat>(M,D1*D2);

	// ctor, cast
	{
		ll_mesh<>::itr itr(tMat.cBegin(),D2);
		CPPUNIT_ASSERT(tMat.cBegin()==itr.vitr());
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)D2,itr.s());

		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),(ptrdiff_t)itr);
		CPPUNIT_ASSERT_EQUAL(size_t(0),(size_t)itr);
	}

	// increment, decrement, comparison
	{
		ll_mesh<>::itr itr(tMat.cBegin(),D2);
		
		++itr;
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)itr.vitr(),itr.s());

		auto jtr1=itr++;
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)jtr1.vitr(),jtr1.s());
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)itr.vitr(),2*itr.s());
		CPPUNIT_ASSERT(jtr1!=itr);
		CPPUNIT_ASSERT(jtr1<itr);
		CPPUNIT_ASSERT(jtr1<=itr);

		--itr;
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)itr.vitr(),itr.s());
		CPPUNIT_ASSERT(jtr1==itr);
		CPPUNIT_ASSERT(jtr1<=itr);
		CPPUNIT_ASSERT(jtr1>=itr);
		
		auto jtr2=itr--;
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)jtr2.vitr(),jtr2.s());
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)itr.vitr(),0*itr.s());
		CPPUNIT_ASSERT(jtr2!=itr);
		CPPUNIT_ASSERT(jtr2>itr);
		CPPUNIT_ASSERT(jtr2>=itr);
	}

	// arithmetic
	{
		const size_t S = genRndST();
		ll_mesh<>::itr itr(tMat.cBegin(),D2);

		itr+=S;
		CPPUNIT_ASSERT_EQUAL(S*itr.s(),(size_t)itr);
		CPPUNIT_ASSERT_EQUAL((size_t)itr.vitr(),(size_t)itr);
		
		auto jtr1 = itr-S;
		CPPUNIT_ASSERT_EQUAL((size_t)0,(size_t)jtr1);

		itr-=S;
		CPPUNIT_ASSERT_EQUAL((size_t)0,(size_t)itr);
		
		auto jtr2 = itr+S;
		CPPUNIT_ASSERT_EQUAL(S*jtr2.s(),(size_t)jtr2);
		
		auto jtr3 = S+itr;
		CPPUNIT_ASSERT_EQUAL(S*jtr3.s(),(size_t)jtr3);
	}

	// difference
	{
		const ptrdiff_t S = genRndST();
		ll_mesh<>::itr itr(tMat.cBegin(),D2);

		CPPUNIT_ASSERT_EQUAL(S,(itr+S)-itr);
	}

	// assignment
	{
		ll_mesh<>::itr itr(tMat.cBegin(),D2);
		ll_mesh<>::itr jtr(tMat.cBegin(),D1);
		
		jtr = itr;
		CPPUNIT_ASSERT(itr==jtr);
	}

	// swap
	{
		ll_mesh<>::itr itr(tMat.cBegin(),D2);
		ll_mesh<>::itr jtr(tMat.cBegin(),D1);
		
		auto icp = itr;
		auto jcp = jtr;

		swap(itr,jtr);
		CPPUNIT_ASSERT(icp==jtr);
		CPPUNIT_ASSERT(jcp==itr);
	}

	// const dereference
	{
		ll_mesh<>::itr itr(tMat.cBegin(),D2);

		CPPUNIT_ASSERT_EQUAL(*itr.vitr(),*itr);
		CPPUNIT_ASSERT_EQUAL(itr.vitr()[0],itr[0]);
		CPPUNIT_ASSERT_EQUAL(itr.vitr()[itr.s()],itr[1]);
		CPPUNIT_ASSERT_EQUAL(itr.vitr()->M(),itr->M());
		CPPUNIT_ASSERT_EQUAL(itr.vitr()->N(),itr->N());
	}

	// reverse iterator
	{
		ll_mesh<>::itr itr(tMat.cBegin(),D2);
		ll_mesh<>::r_itr r_itr(itr);

		CPPUNIT_ASSERT_EQUAL(itr,r_itr.base());
		CPPUNIT_ASSERT_EQUAL(itr,ll_mesh<>::itr(r_itr));
		CPPUNIT_ASSERT(*(r_itr-1)==*itr);
	}
}

void test_mesh_all::test_c_itr() {

	const size_t M = genRndST(1,10);
	const size_t D1 = genRndST(2,10);
	const size_t D2 = genRndST(2,10);

	auto tMat = rand<fMat>(M,D1*D2);

	// ctor, cast
	{
		ll_mesh<>::c_itr itr(tMat.cBegin(),D2);
		CPPUNIT_ASSERT(tMat.cBegin()==itr.vitr());
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)D2,itr.s());

		CPPUNIT_ASSERT_EQUAL(ptrdiff_t(0),(ptrdiff_t)itr);
		CPPUNIT_ASSERT_EQUAL(size_t(0),(size_t)itr);
	}

	// increment, decrement, comparison
	{
		ll_mesh<>::c_itr itr(tMat.cBegin(),D2);
		
		++itr;
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)itr.vitr(),itr.s());

		auto jtr1=itr++;
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)jtr1.vitr(),jtr1.s());
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)itr.vitr(),2*itr.s());
		CPPUNIT_ASSERT(jtr1!=itr);
		CPPUNIT_ASSERT(jtr1<itr);
		CPPUNIT_ASSERT(jtr1<=itr);

		--itr;
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)itr.vitr(),itr.s());
		CPPUNIT_ASSERT(jtr1==itr);
		CPPUNIT_ASSERT(jtr1<=itr);
		CPPUNIT_ASSERT(jtr1>=itr);
		
		auto jtr2=itr--;
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)jtr2.vitr(),jtr2.s());
		CPPUNIT_ASSERT_EQUAL((ptrdiff_t)itr.vitr(),0*itr.s());
		CPPUNIT_ASSERT(jtr2!=itr);
		CPPUNIT_ASSERT(jtr2>itr);
		CPPUNIT_ASSERT(jtr2>=itr);
	}

	// arithmetic
	{
		const size_t S = genRndST();
		ll_mesh<>::c_itr itr(tMat.cBegin(),D2);

		itr+=S;
		CPPUNIT_ASSERT_EQUAL(S*itr.s(),(size_t)itr);
		CPPUNIT_ASSERT_EQUAL((size_t)itr.vitr(),(size_t)itr);
		
		auto jtr1 = itr-S;
		CPPUNIT_ASSERT_EQUAL((size_t)0,(size_t)jtr1);

		itr-=S;
		CPPUNIT_ASSERT_EQUAL((size_t)0,(size_t)itr);
		
		auto jtr2 = itr+S;
		CPPUNIT_ASSERT_EQUAL(S*jtr2.s(),(size_t)jtr2);
		
		auto jtr3 = S+itr;
		CPPUNIT_ASSERT_EQUAL(S*jtr3.s(),(size_t)jtr3);
	}

	// difference
	{
		const ptrdiff_t S = genRndST();
		ll_mesh<>::c_itr itr(tMat.cBegin(),D2);

		CPPUNIT_ASSERT_EQUAL(S,(itr+S)-itr);
	}

	// assignment
	{
		ll_mesh<>::c_itr itr(tMat.cBegin(),D2);
		ll_mesh<>::c_itr jtr(tMat.cBegin(),D1);
		
		jtr = itr;
		CPPUNIT_ASSERT(itr==jtr);
	}

	// swap
	{
		ll_mesh<>::c_itr itr(tMat.cBegin(),D2);
		ll_mesh<>::c_itr jtr(tMat.cBegin(),D1);
		
		auto icp = itr;
		auto jcp = jtr;

		swap(itr,jtr);
		CPPUNIT_ASSERT(icp==jtr);
		CPPUNIT_ASSERT(jcp==itr);
	}

	// const dereference
	{
		ll_mesh<>::c_itr itr(tMat.cBegin(),D2);

		CPPUNIT_ASSERT_EQUAL(*itr.vitr(),*itr);
		CPPUNIT_ASSERT_EQUAL(itr.vitr()[0],itr[0]);
		CPPUNIT_ASSERT_EQUAL(itr.vitr()[itr.s()],itr[1]);
		CPPUNIT_ASSERT_EQUAL(itr.vitr()->M(),itr->M());
		CPPUNIT_ASSERT_EQUAL(itr.vitr()->N(),itr->N());
	}
	
	// reverse iterator
	{
		ll_mesh<>::c_itr itr(tMat.cBegin(),D2);
		ll_mesh<>::cr_itr r_itr(itr);

		CPPUNIT_ASSERT_EQUAL(itr,r_itr.base());
		CPPUNIT_ASSERT_EQUAL(itr,ll_mesh<>::c_itr(r_itr));
		CPPUNIT_ASSERT(*(r_itr-1)==*itr);
	}
}


void test_mesh_all::test_mesh() {
	
	const size_t M = genRndST(1,5);
	const size_t ND = genRndST(1,5);

	std::vector<size_t> D; D.reserve(ND);
	while (D.capacity()!=D.size())
		D.push_back(genRndST(1,5));

	std::vector<size_t> maj(ND);
	{
		std::iota(maj.begin(),maj.end(),0);

		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(maj.begin(),maj.end(),g);
	}
	
	const auto tMat = rand<fMat>(M,
		std::accumulate(D.cbegin(),D.cend(),size_t(1),
			std::multiplies<size_t>()));
	
	// default ctor
	{
		const ll_mesh<> tMesh;
		CPPUNIT_ASSERT(tMesh.empty());
		CPPUNIT_ASSERT_EQUAL(size_t(0),tMesh.size());
		CPPUNIT_ASSERT_EQUAL(size_t(0),tMesh.meshDim());
		CPPUNIT_ASSERT_EQUAL(size_t(0),tMesh.M());
		CPPUNIT_ASSERT_EQUAL(fMat(),tMesh.base_);
		CPPUNIT_ASSERT(tMesh.D().empty());
		CPPUNIT_ASSERT(tMesh.s().empty());
	}
	
	// from fMat
	{
		const ll_mesh<> tMesh(tMat,maj,D);
		CPPUNIT_ASSERT(!tMesh.empty());
		CPPUNIT_ASSERT_EQUAL(tMat.N(),tMesh.size());
		CPPUNIT_ASSERT_EQUAL(D.size(),tMesh.meshDim());
		CPPUNIT_ASSERT_EQUAL(tMat.M(),tMesh.M());
		for (size_t i=0; i!=D.size(); ++i)
			CPPUNIT_ASSERT_EQUAL(D[i],tMesh.D(i));
		CPPUNIT_ASSERT(D==tMesh.D());

		const size_t ref = genRndST(0,tMesh.N()-1);
		const size_t i = tMesh.posToi(tMesh.iTopos(ref));
		CPPUNIT_ASSERT_EQUAL(ref,i);
	}

	// from baseDim
	{
		const ll_mesh<> tMesh(tMat.M(),maj,D);
		CPPUNIT_ASSERT(!tMesh.empty());
		CPPUNIT_ASSERT_EQUAL(tMat.N(),tMesh.size());
		CPPUNIT_ASSERT_EQUAL(D.size(),tMesh.meshDim());
		CPPUNIT_ASSERT_EQUAL(tMat.M(),tMesh.M());
		for (size_t i=0; i!=D.size(); ++i)
			CPPUNIT_ASSERT_EQUAL(D[i],tMesh.D(i));
		CPPUNIT_ASSERT(D==tMesh.D());

		const size_t ref = genRndST(0,tMesh.N()-1);
		const size_t i = tMesh.posToi(tMesh.iTopos(ref));
		CPPUNIT_ASSERT_EQUAL(ref,i);
	}

	// check begin,end
	{
		const ll_mesh<> tMesh(tMat,maj,D);

		std::vector<size_t> pos(tMesh.meshDim());
		for (size_t d=0; d!=tMesh.meshDim(); ++d)
			pos[d] = genRndST(0,tMesh.D(d)-1);

		for (size_t d=0; d!=tMesh.meshDim(); ++d) {
			
			const auto itr = tMesh.cbegin(d,pos);
			auto pos_ = pos; pos_[d]=0;
			CPPUNIT_ASSERT(pos_==tMesh.iTopos((ptrdiff_t)itr));
			
			const auto etr = tMesh.cend(d,pos);
			CPPUNIT_ASSERT_EQUAL(ptrdiff_t(tMesh.D(d)),etr-itr);
			CPPUNIT_ASSERT_EQUAL(-ptrdiff_t(tMesh.D(d)),itr-etr);
		}
	}
}

void test_mesh_all::test_generators() {
	std::random_device rd;
	std::mt19937 g(rd());
	
	// 1D
	{
		const auto tMesh1 = genMesh(fMat({0.0,1.0},1),{0},{1});
		CPPUNIT_ASSERT_EQUAL(fMat({0.0},1),tMesh1.base_);
		CPPUNIT_ASSERT(std::vector<size_t>{1}==tMesh1.D());
		CPPUNIT_ASSERT(std::vector<size_t>{1}==tMesh1.s());
		
		const auto tMesh2 = genMesh(fMat({0.0,1.0},1),{0},{10});
		CPPUNIT_ASSERT_EQUAL(0.0,tMesh2.base_.front());
		CPPUNIT_ASSERT_EQUAL(1.0,tMesh2.base_.back());
		auto db = tMesh2.base_.get({},{1,2,3,4,5,6,7,8,9}) - 
			  tMesh2.base_.get({},{0,1,2,3,4,5,6,7,8});
		CPPUNIT_ASSERT_EQUAL(db[0]*ones<fMat>(1,9),db);
	}

	// 2D
	{
		fMat bounds({0.0,0.0,1.0,1.0},2);
		std::vector<size_t> N = {genRndST(1,10),genRndST(1,10)};
		std::vector<size_t> maj = {0,1};
		std::shuffle(maj.begin(),maj.end(),g);

		const auto tMesh1 = genMesh(bounds,maj,N);

		{
			fMat ck(2,0); ck.reserve(tMesh1.size());
			for (size_t j=0; j!=tMesh1.D(1); ++j)
			for (auto i=tMesh1.cbegin(0,{0,j}),e=tMesh1.cend(0,{0,j}); i!=e; ++i)
				ck.push_back(*i);
			
			for (auto i=tMesh1.base_.ccBegin(),e=tMesh1.base_.ccEnd(); i!=e; ++i)
				CPPUNIT_ASSERT(std::find(ck.ccBegin(),ck.ccEnd(),*i)!=ck.ccEnd());
		}
		{
			fMat ck(2,0); ck.reserve(tMesh1.size());
			for (size_t j=0; j!=tMesh1.D(0); ++j)
			for (auto i=tMesh1.cbegin(1,{j,0}),e=tMesh1.cend(1,{j,0}); i!=e; ++i)
				ck.push_back(*i);

			for (auto i=tMesh1.base_.ccBegin(),e=tMesh1.base_.ccEnd(); i!=e; ++i)
				CPPUNIT_ASSERT(std::find(ck.ccBegin(),ck.ccEnd(),*i)!=ck.ccEnd());
		}
	}

	// ND
	{
		const size_t N_ = genRndST(3,3);

		auto bounds = ones<fMat>(N_,2);
		bounds.cFront()*=genRndDouble();
		bounds.cBack()*=genRndDouble();

		std::vector<size_t> N; N.reserve(N_);
		while (N.size()!=N.capacity())
			N.push_back(genRndST(1,5));

		std::vector<size_t> maj(N_);
		std::iota(maj.begin(),maj.end(),0);
		std::shuffle(maj.begin(),maj.end(),g);

		const auto tMesh1 = genMesh(bounds,maj,N);

		fMat ck(N_,0); ck.reserve(tMesh1.size());
		if (N_==3)
			for (size_t j1=0; j1!=tMesh1.D(1); ++j1)
			for (size_t j2=0; j2!=tMesh1.D(2); ++j2)
			for (auto i=tMesh1.cbegin(0,{0,j1,j2}),
				  e=tMesh1.cend(0,{0,j1,j2}); i!=e; ++i) {
				ck.push_back(*i);
			}
		if (N_==4)
			for (size_t j1=0; j1!=tMesh1.D(1); ++j1)
			for (size_t j2=0; j2!=tMesh1.D(2); ++j2)
			for (size_t j3=0; j3!=tMesh1.D(3); ++j3)
			for (auto i=tMesh1.cbegin(0,{0,j1,j2,j3}),
				  e=tMesh1.cend(0,{0,j1,j2,j3}); i!=e; ++i)
				ck.push_back(*i);
		if (N_==5)
			for (size_t j1=0; j1!=tMesh1.D(1); ++j1)
			for (size_t j2=0; j2!=tMesh1.D(2); ++j2)
			for (size_t j3=0; j3!=tMesh1.D(3); ++j3)
			for (size_t j4=0; j4!=tMesh1.D(4); ++j4)
			for (auto i=tMesh1.cbegin(0,{0,j1,j2,j3,j4}),
				  e=tMesh1.cend(0,{0,j1,j2,j3,j4}); i!=e; ++i)
				ck.push_back(*i);

		for (auto i=tMesh1.base_.ccBegin(),e=tMesh1.base_.ccEnd(); i!=e; ++i)
			CPPUNIT_ASSERT(std::find(ck.ccBegin(),ck.ccEnd(),*i)!=ck.ccEnd());
	}

	// 2D mesh MATLAB style, write to file
	{
		const auto tMesh = genMesh(0.0,1.0,2,{1,0},{false,false});
		tMesh.writeToFile("outp/tMesh.bin");
	}

	{
		const auto B = ones<fMat>(3)-eye<fMat>(3);
		const double rho = 1000.0;
		const double lb = 0.0;
		const double ub = 1.0;
		std::vector<size_t> maj = {1,0,2};

		{
			const auto tMesh = genMesh_cell(rho,B,lb,ub,maj,{false,false,false});
			tMesh.writeToFile("outp/tMesh_000__.bin");
		}
		{
			const auto tMesh = genMesh_cell(rho,B,lb,ub,maj,{true,false,false});
			tMesh.writeToFile("outp/tMesh_100.bin");
		}
		{
			const auto tMesh = genMesh_cell(rho,B,lb,ub,maj,{false,true,false});
			tMesh.writeToFile("outp/tMesh_010.bin");
		}
		{
			const auto tMesh = genMesh_cell(rho,B,lb,ub,maj,{false,false,true});
			tMesh.writeToFile("outp/tMesh_001.bin");
		}
		{
			const auto tMesh = genMesh_cell(rho,B,lb,ub,maj,{true,true,true});
			tMesh.writeToFile("outp/tMesh_111.bin");
		}
	}
}


const char* test_mesh_all::test_id() noexcept {
	return "test_mesh_all";
}

CppUnit::Test* test_mesh_all::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_mesh_all>(
		"test_itr", &test_mesh_all::test_itr));
	suite->addTest(new CppUnit::TestCaller<test_mesh_all>(
		"test_c_itr", &test_mesh_all::test_c_itr));
	suite->addTest(new CppUnit::TestCaller<test_mesh_all>(
		"test_mesh", &test_mesh_all::test_mesh));
	suite->addTest(new CppUnit::TestCaller<test_mesh_all>(
		"test_generators", &test_mesh_all::test_generators));
	
	return suite;
}
