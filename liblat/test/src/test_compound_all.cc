// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_compound_all.h"
#include "ll_compound.h"
#include "ll_testTools.h"
#include "lm_testTools.h"

using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;

void test_compound_all::test_mat_cb() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	mat_cb te;
	CPPUNIT_ASSERT(te.empty());
	
	fMat rnd = rand<fMat>(M,N);
	fMat cpy = rnd;
	auto ptr = cpy.data();

	mat_cb t(std::move(cpy));
	CPPUNIT_ASSERT(rnd==t.mat_);
	CPPUNIT_ASSERT_EQUAL(ptr,t.mat_.data());
	CPPUNIT_ASSERT_EQUAL(M,t.dim());
	CPPUNIT_ASSERT_EQUAL(N,t.N());
	CPPUNIT_ASSERT_EQUAL(N,t.capacity());

	const size_t n = genRndST(0,N-1);
	CPPUNIT_ASSERT_EQUAL(rnd.cAt(n),t.cAt(n));
	CPPUNIT_ASSERT_EQUAL(rnd.cFront(),t.cFront());
	CPPUNIT_ASSERT_EQUAL(rnd.cBack(),t.cBack());

	CPPUNIT_ASSERT_EQUAL(t.cFront(),*t.cBegin());
	CPPUNIT_ASSERT_EQUAL(t.cBack(),*(t.cEnd()-1));
	CPPUNIT_ASSERT_EQUAL(t.cFront(),*t.ccBegin());
	CPPUNIT_ASSERT_EQUAL(t.cBack(),*(t.ccEnd()-1));

	t.reserve(2*t.N());
	CPPUNIT_ASSERT_EQUAL(2*N,t.capacity());
	for (auto i=rnd.ccBegin(),e=rnd.ccEnd(); i!=e; ++i)
		t.push_back(*i);
	CPPUNIT_ASSERT_EQUAL(2*N,t.N());
	CPPUNIT_ASSERT_EQUAL(2*N,t.capacity());
	for (size_t i=0; i!=N; ++i)
		CPPUNIT_ASSERT_EQUAL(rnd.cAt(i),t.cAt(i+N));

	t.resize(N+1);
	CPPUNIT_ASSERT_EQUAL(N+1,t.N());
	t.pop_back();
	CPPUNIT_ASSERT_EQUAL(N,t.N());
	t.shrink_to_fit();
	CPPUNIT_ASSERT_EQUAL(N,t.capacity());
}
void test_compound_all::test_mat_b() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	mat_b te;
	CPPUNIT_ASSERT(te.empty());
		
	fMat rnd = rand<fMat>(M,N);
	fMat cpy = rnd;
	auto ptr = cpy.data();

	mat_b t(std::move(cpy));
	CPPUNIT_ASSERT(rnd==t.mat_);
	CPPUNIT_ASSERT_EQUAL(ptr,t.data());
	CPPUNIT_ASSERT_EQUAL(M,t.dim());
	CPPUNIT_ASSERT_EQUAL(N,t.N());
	CPPUNIT_ASSERT_EQUAL(N,t.capacity());

	const size_t n = genRndST(0,N-1);
	CPPUNIT_ASSERT_EQUAL(rnd.cAt(n),t.cAt(n));
	CPPUNIT_ASSERT_EQUAL(rnd.cFront(),t.cFront());
	CPPUNIT_ASSERT_EQUAL(rnd.cBack(),t.cBack());

	CPPUNIT_ASSERT_EQUAL(t.cFront(),*t.cBegin());
	CPPUNIT_ASSERT_EQUAL(t.cBack(),*(t.cEnd()-1));
	CPPUNIT_ASSERT_EQUAL(t.cFront(),*t.ccBegin());
	CPPUNIT_ASSERT_EQUAL(t.cBack(),*(t.ccEnd()-1));

	t.reserve(2*t.N());
	CPPUNIT_ASSERT_EQUAL(2*N,t.capacity());
	for (auto i=rnd.ccBegin(),e=rnd.ccEnd(); i!=e; ++i)
		t.push_back(*i);
	CPPUNIT_ASSERT_EQUAL(2*N,t.N());
	CPPUNIT_ASSERT_EQUAL(2*N,t.capacity());
	for (size_t i=0; i!=N; ++i)
		CPPUNIT_ASSERT_EQUAL(rnd.cAt(i),t.cAt(i+N));

	t.resize(N+1);
	CPPUNIT_ASSERT_EQUAL(N+1,t.N());
	t.pop_back();
	CPPUNIT_ASSERT_EQUAL(N,t.N());
	t.shrink_to_fit();
	CPPUNIT_ASSERT_EQUAL(N,t.capacity());
}
void test_compound_all::test_vec_cb() {
	const size_t N = genRndST();

	vec_cb ve;
	CPPUNIT_ASSERT(ve.empty());
		
	std::vector<size_t> rnd(N);
	std::iota(rnd.begin(),rnd.end(),0);
	std::random_shuffle(rnd.begin(),rnd.end());
	std::vector<size_t> cpy = rnd;
	auto ptr = cpy.data();

	vec_cb v(std::move(cpy));
	CPPUNIT_ASSERT(rnd==v.vec_);
	CPPUNIT_ASSERT_EQUAL((size_t)ptr,(size_t)v.data());
	CPPUNIT_ASSERT_EQUAL(N,v.size());
	CPPUNIT_ASSERT_EQUAL(N,v.capacity());

	const size_t n = genRndST(0,N-1);
	CPPUNIT_ASSERT_EQUAL(rnd[n],v[n]);
	CPPUNIT_ASSERT_EQUAL(rnd.front(),v.front());
	CPPUNIT_ASSERT_EQUAL(rnd.back(),v.back());

	CPPUNIT_ASSERT_EQUAL(v.front(),*v.begin());
	CPPUNIT_ASSERT_EQUAL(v.back(),*(v.end()-1));
	CPPUNIT_ASSERT_EQUAL(v.front(),*v.cbegin());
	CPPUNIT_ASSERT_EQUAL(v.back(),*(v.cend()-1));

	v.reserve(2*v.size());
	CPPUNIT_ASSERT_EQUAL(2*N,v.capacity());
	for (auto i=rnd.cbegin(),e=rnd.cend(); i!=e; ++i)
		v.push_back(*i);
	CPPUNIT_ASSERT_EQUAL(2*N,v.size());
	CPPUNIT_ASSERT_EQUAL(2*N,v.capacity());
	for (size_t i=0; i!=N; ++i)
		CPPUNIT_ASSERT_EQUAL(rnd[i],v[i+N]);

	v.resize(N+1);
	CPPUNIT_ASSERT_EQUAL(N+1,v.size());
	v.pop_back();
	CPPUNIT_ASSERT_EQUAL(N,v.size());
	v.shrink_to_fit();
	CPPUNIT_ASSERT_EQUAL(N,v.capacity());
}
void test_compound_all::test_vec_b() {
	const size_t N = genRndST();

	vec_b ve;
	CPPUNIT_ASSERT(ve.empty());
		
	std::vector<size_t> rnd(N);
	std::iota(rnd.begin(),rnd.end(),0);
	std::random_shuffle(rnd.begin(),rnd.end());
	std::vector<size_t> cpy = rnd;
	auto ptr = cpy.data();

	vec_b v(std::move(cpy));
	CPPUNIT_ASSERT(rnd==v.vec_);
	CPPUNIT_ASSERT_EQUAL((size_t)ptr,(size_t)v.data());
	CPPUNIT_ASSERT_EQUAL(N,v.size());
	CPPUNIT_ASSERT_EQUAL(N,v.capacity());

	const size_t n = genRndST(0,N-1);
	CPPUNIT_ASSERT_EQUAL(rnd[n],v[n]);
	CPPUNIT_ASSERT_EQUAL(rnd.front(),v.front());
	CPPUNIT_ASSERT_EQUAL(rnd.back(),v.back());

	CPPUNIT_ASSERT_EQUAL(v.front(),*v.begin());
	CPPUNIT_ASSERT_EQUAL(v.back(),*(v.end()-1));
	CPPUNIT_ASSERT_EQUAL(v.front(),*v.cbegin());
	CPPUNIT_ASSERT_EQUAL(v.back(),*(v.cend()-1));

	v.reserve(2*v.size());
	CPPUNIT_ASSERT_EQUAL(2*N,v.capacity());
	for (auto i=rnd.cbegin(),e=rnd.cend(); i!=e; ++i)
		v.push_back(*i);
	CPPUNIT_ASSERT_EQUAL(2*N,v.size());
	CPPUNIT_ASSERT_EQUAL(2*N,v.capacity());
	for (size_t i=0; i!=N; ++i)
		CPPUNIT_ASSERT_EQUAL(rnd[i],v[i+N]);

	v.resize(N+1);
	CPPUNIT_ASSERT_EQUAL(N+1,v.size());
	v.pop_back();
	CPPUNIT_ASSERT_EQUAL(N,v.size());
	v.shrink_to_fit();
	CPPUNIT_ASSERT_EQUAL(N,v.capacity());
}
void test_compound_all::test_mat_vec_cb() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	mat_vec_cb tve;
	CPPUNIT_ASSERT(tve.empty());
		
	fMat mrnd = rand<fMat>(M,N);
	fMat mcpy = mrnd;
	auto mptr = mcpy.data();
	
	std::vector<size_t> vrnd(N);
	std::iota(vrnd.begin(),vrnd.end(),0);
	std::random_shuffle(vrnd.begin(),vrnd.end());
	std::vector<size_t> vcpy = vrnd;
	auto vptr = vcpy.data();

	mat_vec_cb tv(std::move(mcpy),std::move(vcpy));
	CPPUNIT_ASSERT(mrnd==tv.mat_);
	CPPUNIT_ASSERT(vrnd==tv.vec_);
	CPPUNIT_ASSERT_EQUAL(mptr,tv.mdata());
	CPPUNIT_ASSERT_EQUAL(vptr,tv.vdata());
	CPPUNIT_ASSERT_EQUAL(M,tv.dim());
	CPPUNIT_ASSERT_EQUAL(N,tv.N());
	CPPUNIT_ASSERT_EQUAL(N,tv.size());
	CPPUNIT_ASSERT_EQUAL(N,tv.capacity());

	tv.reserve(2*tv.size());
	CPPUNIT_ASSERT_EQUAL(2*N,tv.capacity());
	auto j = mrnd.ccBegin();
	for (auto i=vrnd.cbegin(),e=vrnd.cend(); i!=e; ++i,++j)
		tv.push_back(*j,*i);
	CPPUNIT_ASSERT_EQUAL(2*N,tv.size());
	CPPUNIT_ASSERT_EQUAL(2*N,tv.capacity());
	for (size_t i=0; i!=N; ++i) {
		CPPUNIT_ASSERT_EQUAL(mrnd.cAt(i),tv.cAt(i+N));
		CPPUNIT_ASSERT_EQUAL(vrnd[i],tv[i+N]);
	}
	
	tv.resize(N+1);
	CPPUNIT_ASSERT_EQUAL(N+1,tv.N());
	CPPUNIT_ASSERT_EQUAL(N+1,tv.size());
	tv.pop_back();
	CPPUNIT_ASSERT_EQUAL(N,tv.N());
	CPPUNIT_ASSERT_EQUAL(N,tv.size());
	tv.shrink_to_fit();
	CPPUNIT_ASSERT_EQUAL(N,tv.capacity());
}
void test_compound_all::test_mat_vec_b() {
	const size_t M = genRndST();
	const size_t N = genRndST();

	mat_vec_cb tve;
	CPPUNIT_ASSERT(tve.empty());
		
	fMat mrnd = rand<fMat>(M,N);
	fMat mcpy = mrnd;
	auto mptr = mcpy.data();
	
	std::vector<size_t> vrnd(N);
	std::iota(vrnd.begin(),vrnd.end(),0);
	std::random_shuffle(vrnd.begin(),vrnd.end());
	std::vector<size_t> vcpy = vrnd;
	auto vptr = vcpy.data();

	mat_vec_cb tv(std::move(mcpy),std::move(vcpy));
	CPPUNIT_ASSERT(mrnd==tv.mat_);
	CPPUNIT_ASSERT(vrnd==tv.vec_);
	CPPUNIT_ASSERT_EQUAL(mptr,tv.mdata());
	CPPUNIT_ASSERT_EQUAL(vptr,tv.vdata());
	CPPUNIT_ASSERT_EQUAL(M,tv.dim());
	CPPUNIT_ASSERT_EQUAL(N,tv.N());
	CPPUNIT_ASSERT_EQUAL(N,tv.size());
	CPPUNIT_ASSERT_EQUAL(N,tv.capacity());

	tv.reserve(2*tv.size());
	CPPUNIT_ASSERT_EQUAL(2*N,tv.capacity());
	auto j = mrnd.ccBegin();
	for (auto i=vrnd.cbegin(),e=vrnd.cend(); i!=e; ++i,++j)
		tv.push_back(*j,*i);
	CPPUNIT_ASSERT_EQUAL(2*N,tv.size());
	CPPUNIT_ASSERT_EQUAL(2*N,tv.capacity());
	for (size_t i=0; i!=N; ++i) {
		CPPUNIT_ASSERT_EQUAL(mrnd.cAt(i),tv.cAt(i+N));
		CPPUNIT_ASSERT_EQUAL(vrnd[i],tv[i+N]);
	}
	
	tv.resize(N+1);
	CPPUNIT_ASSERT_EQUAL(N+1,tv.N());
	CPPUNIT_ASSERT_EQUAL(N+1,tv.size());
	tv.pop_back();
	CPPUNIT_ASSERT_EQUAL(N,tv.N());
	CPPUNIT_ASSERT_EQUAL(N,tv.size());
	tv.shrink_to_fit();
	CPPUNIT_ASSERT_EQUAL(N,tv.capacity());
}


void test_compound_all::test_compound_classes() {
	const size_t D = genRndST(1,10);
	const size_t N = genRndST(20,30);
	
	auto tMat1 = rand<fMat>(D,N);
	auto mVec1 = rand<fMat>(1,N);

	aTv tVec1(N,0);
	std::iota(tVec1.begin(),tVec1.end(),aT(0));
	std::random_shuffle(tVec1.begin(),tVec1.end());

	idv tVec2 = genRndIdv(N);
		
	const std::string tStr = "abcdefghiklmnopqrstuvwxyz1234567890";
	
	// Ap_T
	{
		const Ap_T tComp(tMat1,tVec1);
		CPPUNIT_ASSERT(tComp.mat_==tComp.Ap());
		CPPUNIT_ASSERT(tComp.vec_==tComp.T());
	}

	// Ap_id
	{
		const Ap_id tComp(tMat1,tVec2);
		CPPUNIT_ASSERT(tComp.mat_==tComp.Ap());
		CPPUNIT_ASSERT(tComp.vec_==tComp.id());

	}

	// Wp_s
	{
		const Wp_s tComp(tMat1,mVec1);
		CPPUNIT_ASSERT(tComp.mat_==tComp.Wp());
		CPPUNIT_ASSERT(tComp.vec_==tComp.s());

	}

	// R_H
	{
		const size_t Nw = genRndST();
		cMatv hVec(N);
		for (auto& h: hVec) h = rand<cMat>(Nw,Nw);
		auto tMat2 = randi<fMat>(D,N);
		std::sort(tMat2.cBegin(),tMat2.cEnd(),vcmp);

		R_H<> tComp(std::move(tMat2),std::move(hVec));
		CPPUNIT_ASSERT(tComp.mat_==tComp.R());
		CPPUNIT_ASSERT(tComp.vec_==tComp.H());
		CPPUNIT_ASSERT(tComp[N/2]==tComp.center());
		CPPUNIT_ASSERT_EQUAL(Nw,tComp.Nw());
	}

	// k_U
	{
		const size_t Nw = genRndST();
		const size_t Nb = genRndST();

		std::vector<cMat> uVec(N);
		for (auto& u: uVec) u = rand<cMat>(Nb,Nw);

		k_U tComp(tMat1,std::move(uVec));
		CPPUNIT_ASSERT(tComp.mat_==tComp.k());
		CPPUNIT_ASSERT(tComp.vec_==tComp.U());
		CPPUNIT_ASSERT_EQUAL(Nw,tComp.Nw());
		CPPUNIT_ASSERT_EQUAL(Nb,tComp.Nb());
	}

	// p_p
	{
		std::vector<size_t> Nps;
		while (std::accumulate(Nps.cbegin(),Nps.cend(),size_t(0))<N)
			Nps.push_back(genRndST(2,5));
		if (std::accumulate(Nps.cbegin(),Nps.cend(),size_t(0))!=N) {
			Nps.pop_back();
			Nps.push_back(N-std::accumulate(Nps.cbegin(),Nps.cend(),size_t(0)));
		}

		p_p tComp(tMat1,mVec1,Nps);
		CPPUNIT_ASSERT(tComp.mat_==tComp.path());
		CPPUNIT_ASSERT(tComp.vec_==tComp.pos());
		CPPUNIT_ASSERT(Nps==tComp.Nps());
	}
}


const char* test_compound_all::test_id() noexcept {
	return "test_compound_all";
}

CppUnit::Test* test_compound_all::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_compound_all>(
		"test_mat_cb", &test_compound_all::test_mat_cb));
	suite->addTest(new CppUnit::TestCaller<test_compound_all>(
		"test_mat_b", &test_compound_all::test_mat_b));
	suite->addTest(new CppUnit::TestCaller<test_compound_all>(
		"test_vec_cb", &test_compound_all::test_vec_cb));
	suite->addTest(new CppUnit::TestCaller<test_compound_all>(
		"test_vec_b", &test_compound_all::test_vec_b));
	suite->addTest(new CppUnit::TestCaller<test_compound_all>(
		"test_mat_vec_cb", &test_compound_all::test_mat_vec_cb));
	suite->addTest(new CppUnit::TestCaller<test_compound_all>(
		"test_mat_vec_b", &test_compound_all::test_mat_vec_b));
	suite->addTest(new CppUnit::TestCaller<test_compound_all>(
		"test_compound_classes", &test_compound_all::test_compound_classes));
	
	return suite;
}
