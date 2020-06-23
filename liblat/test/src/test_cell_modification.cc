// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_cell_modification.h"
#include "ll_testTools.h"
#include "aux_io.h"
#include "libmat.h"
#include <iostream>
#include <algorithm>
#include "ll_io.h"
#include "ll_fn.h"

using namespace ll__;
using namespace ll__::test;
using namespace lm__;
using namespace lm__::test;
using namespace aux;


void test_cell_modification::test_stress() {
	const auto D = genRndST(1,5);
	auto tCell1 = genRandom(D);
	const auto tCell2 = tCell1;

	fMat S;
	do S = rand<fMat>(D,D);
	while (det(S)<=mtol());

	tCell1.stress(S);
	CPPUNIT_ASSERT(S.prod(tCell2.B())==tCell1.B());
}

void test_cell_modification::test_swapDim() {
	const auto D = genRndST(1,5);
	auto tCell1 = genRandom(D);
	const auto tCell2 = tCell1;

	const auto n1 = genRndST(0,D-1);
	const auto n2 = genRndST(0,D-1);

	tCell1.swapDim(n1,n2);

	CPPUNIT_ASSERT(tCell2.B().cAt(n1)==tCell1.B().cAt(n2));
	CPPUNIT_ASSERT(tCell2.B().cAt(n2)==tCell1.B().cAt(n1));
	CPPUNIT_ASSERT_DELTA(tCell2.vol(),tCell1.vol(),mtol());
	CPPUNIT_ASSERT_EQUAL(tCell2.N(),tCell1.N());
	CPPUNIT_ASSERT(tCell2.Ntype()==tCell1.Ntype());
	CPPUNIT_ASSERT(tCell1.sameLattice(tCell2));
	CPPUNIT_ASSERT(Ap_sorted(tCell1));
}

void test_cell_modification::test_invDim() {
	const auto D = genRndST(1,5);
	auto tCell1 = genRandom(D);
	const auto tCell2 = tCell1;
	
	const auto n = genRndST(0,D-1);
	tCell1.invDim(n);
	
	CPPUNIT_ASSERT(tCell2.B().cAt(n)==(-tCell1.B().cAt(n)));
	CPPUNIT_ASSERT_DELTA(tCell2.vol(),tCell1.vol(),mtol());
	CPPUNIT_ASSERT_EQUAL(tCell2.N(),tCell1.N());
	CPPUNIT_ASSERT(tCell2.Ntype()==tCell1.Ntype());

	CPPUNIT_ASSERT(tCell1.sameLattice(tCell2));
	CPPUNIT_ASSERT(Ap_sorted(tCell1));
}

void test_cell_modification::test_orient() {
	const auto D = genRndST(1,5);
	const auto r = genRndrv(D);
	
	auto tCell1 = genRandom(D);

	// sgn 1.0
	if (std::all_of(r.begin(),r.end(),[](const bool i){return i;})) {
		const auto tCell2 = tCell1;
		tCell1.orient(1.0,r);
		
		CPPUNIT_ASSERT(tCell2==tCell1);
	} else {
		const auto tCell2 = tCell1;
		tCell1.orient(1.0,r);

		CPPUNIT_ASSERT_EQUAL(double(1.0),tCell1.sign());
		CPPUNIT_ASSERT_DELTA(tCell2.vol(),tCell1.vol(),mtol());
		CPPUNIT_ASSERT_EQUAL(tCell2.N(),tCell1.N());
		CPPUNIT_ASSERT(tCell2.Ntype()==tCell1.Ntype());
		CPPUNIT_ASSERT(tCell2.sameLattice(tCell1));
	}
	CPPUNIT_ASSERT(Ap_sorted(tCell1));

	// sgn -1.0
	if (std::all_of(r.begin(),r.end(),[](const bool i){return i;})) {
		const auto tCell2 = tCell1;
		tCell1.orient(-1.0,r);
		
		CPPUNIT_ASSERT(tCell2==tCell1);
	} else {
		const auto tCell2 = tCell1;
		tCell1.orient(-1.0,r);

		CPPUNIT_ASSERT_EQUAL(double(-1.0),tCell1.sign());
		CPPUNIT_ASSERT_DELTA(tCell2.vol(),tCell1.vol(),mtol());
		CPPUNIT_ASSERT_EQUAL(tCell2.N(),tCell1.N());
		CPPUNIT_ASSERT(tCell2.Ntype()==tCell1.Ntype());
		CPPUNIT_ASSERT(tCell2.sameLattice(tCell1));
	}
	CPPUNIT_ASSERT(Ap_sorted(tCell1));
}

void test_cell_modification::test_permute() {
	const auto D = genRndST(1,5);
	auto tCell1 = genRandom(D);
	const auto tCell2 = tCell1;

	const auto n1 = genRndST(0,D-1);
	const auto n2 = genRndST(0,D-1);

	auto P = eye<fMat>(D);
	swap(P.cAt(n1),P.cAt(n2));

	tCell1.swapDim(n1,n2);
	tCell1.permute(P);

	CPPUNIT_ASSERT(Ap_sorted(tCell1));
	CPPUNIT_ASSERT(tCell1==tCell2);
}

void test_cell_modification::test_rotate() {
	const auto D = genRndST(1,5);
	auto tCell1 = genRandom(D);
	const auto tCell2 = tCell1;

	const auto R = gsorth(rand<fMat>(D,D));
	
	tCell1.rotate(R);
	CPPUNIT_ASSERT(R.prod(tCell2.B())==tCell1.B());
	CPPUNIT_ASSERT(Ap_sorted(tCell1));

	tCell1.rotate(T(R));
	CPPUNIT_ASSERT(tCell2==tCell1);
}

void test_cell_modification::test_scale() {
	const auto D = genRndST(1,5);
	auto tCell1 = genRandom(D);
	const auto tCell2 = tCell1;

	// single double version
	{
		double f;
		do f = genRndDouble(-2.0,2.0);
		while (std::abs(f)<1e-3);

		tCell1.scale(f);
		CPPUNIT_ASSERT(tCell2.B()*f==tCell1.B());
		CPPUNIT_ASSERT(Ap_sorted(tCell1));

		tCell1.scale(1.0/f);
		CPPUNIT_ASSERT(tCell2==tCell1);
	}

	// multiple double version
	{
		fMat f;
		do f = rand<fMat>(D,1,-2.0,2.0);
		while (any(abs(f).lt(1e-3)));

		tCell1.scale(f);
		auto j=f.cbegin();
		for (auto i1=tCell1.B().ccBegin(),i2=tCell2.B().ccBegin(),
				e1=tCell1.B().ccEnd(); i1!=e1; ++i1,++i2,++j)
			CPPUNIT_ASSERT((*i2 * *j) == *i1);

		tCell1.scale(1.0/f);
		CPPUNIT_ASSERT(tCell2==tCell1);
	}
}

void test_cell_modification::test_changeBasis() {
	const auto D = genRndST(1,5);
	auto tCell1 = genRandom(D);
	
	const auto C = rnd_i(D,-2,2);
	const auto NB = tCell1.B().prod(C);

	auto tCell2 = tCell1;
	tCell2.changeBasis(NB);
	CPPUNIT_ASSERT(Ap_sorted(tCell2));

	const size_t f = std::round(std::abs(det(C)));
		
	CPPUNIT_ASSERT(NB==tCell2.B());
	CPPUNIT_ASSERT_EQUAL(size_t(f*tCell1.N()),tCell2.N());

	CPPUNIT_ASSERT_EQUAL(tCell1.Nspecies(),tCell2.Nspecies());
	for (size_t i=0; i!=tCell2.Nspecies(); ++i)
		CPPUNIT_ASSERT_EQUAL(size_t(f*tCell1.Ntype(i)),tCell2.Ntype(i));

	CPPUNIT_ASSERT(cunique(tCell2.Ap()));

	tCell2.changeBasis(tCell1.B());
	CPPUNIT_ASSERT(Ap_sorted(tCell2));
	CPPUNIT_ASSERT(tCell2==tCell1);
}

void test_cell_modification::test_makePrimitive() {

	const auto D = genRndST(1,5);
	const auto tCell1 = genRandom(D);

	// with r
	{
		auto tCell2 = tCell1;
		
		const auto NB = tCell2.B().prod(rnd_i(D,-2,2));
		tCell2.changeBasis(NB);
		const auto sgn = tCell2.sign();
		const auto r = genRndrv(D);
		const auto tCell3 = tCell2;

		tCell2.makePrimitive(r);

		CPPUNIT_ASSERT(Ap_sorted(tCell2));
		CPPUNIT_ASSERT(tCell2.vol()<=tCell3.vol()+mtol());
		CPPUNIT_ASSERT(tCell2.N()<=tCell3.N());
		CPPUNIT_ASSERT(tCell2.validBasis(tCell1.B()));
		CPPUNIT_ASSERT(tCell2.validBasis(tCell3.B()));
		CPPUNIT_ASSERT_EQUAL(sgn,tCell2.sign());

		// check restricted dimensions are conserved
		const auto ri = inds(r);
		for (const auto i: ri)
			CPPUNIT_ASSERT(tCell2.B().cAt(i)==NB.cAt(i));
	
		tCell2.changeBasis(tCell1.B());
		CPPUNIT_ASSERT(tCell1==tCell2);
	}

	// orthogonal complement
	{
		const size_t D = 3;
		const auto tCell1 = genZincblende(D,1.0);
		auto tCell2 = tCell1;
		const auto NB = diag(randi<fMat>(D,1,1.0,5.0));

		tCell2.changeBasis(NB);
		const auto sgn = tCell2.sign();
		const auto r = genRndrv(D);
		const auto tCell3 = tCell2;

		tCell2.makePrimitiveInSubspace(r);
		
		CPPUNIT_ASSERT(Ap_sorted(tCell2));
		CPPUNIT_ASSERT(tCell2.vol()<=tCell3.vol()+mtol());
		CPPUNIT_ASSERT(tCell2.N()<=tCell3.N());
		CPPUNIT_ASSERT(tCell2.validBasis(tCell1.B()));
		CPPUNIT_ASSERT(tCell2.validBasis(tCell3.B()));
		CPPUNIT_ASSERT_EQUAL(sgn,tCell2.sign());
		
		// check restricted dimensions are conserved
		const auto ri = inds(r);
		for (const auto i: ri)
			CPPUNIT_ASSERT(tCell2.B().cAt(i)==NB.cAt(i));

		// check new vectors are orthogonal to old ones
		const auto nri = ninds(r);
		for (const auto ni: nri)
			for (const auto i: ri)
				CPPUNIT_ASSERT_DELTA(0.0,dot(tCell2.B().cAt(ni),tCell2.B().cAt(i)),mtol());
	}

	// full restricted
	{
		auto tCell2 = tCell1;

		tCell2.changeBasis(tCell2.B().prod(rnd_i(D,-2,2)));
		const auto tCell3 = tCell2;
		tCell2.makePrimitive(rv(D,true));

		CPPUNIT_ASSERT(tCell3==tCell2);
	}
}

void test_cell_modification::test_shift() {
	const auto D = genRndST(1,5);
	auto tCell1 = genRandom(D);
	const auto tCell2 = tCell1;

	const auto sh = rand<fMat>(D,1,.2,.8);
	
	tCell1.shift(sh);
	CPPUNIT_ASSERT(Ap_sorted(tCell1));
	CPPUNIT_ASSERT(tCell2.B()==tCell1.B());
	CPPUNIT_ASSERT(tCell2.types()==tCell1.types());
	CPPUNIT_ASSERT(tCell2.Ntype()==tCell1.Ntype());

	tCell1.shift(-sh);
	CPPUNIT_ASSERT(Ap_sorted(tCell1));
	CPPUNIT_ASSERT(tCell2==tCell1);
}

void test_cell_modification::test_diversify() {
	
	const auto D = genRndST(1,5);

	// no ts
	{
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;

		tCell1.diversify();
		CPPUNIT_ASSERT(tCell1.B()==tCell2.B());
		CPPUNIT_ASSERT(tCell1.Ap()==tCell2.Ap());
		CPPUNIT_ASSERT(tCell1.types()==rg(0,1,tCell1.N()));
		CPPUNIT_ASSERT(tCell1.Ntype()==aCv(tCell1.N(),1));
		CPPUNIT_ASSERT(Ap_sorted(tCell1));
	}

	// no ts, enlarged cell
	{
		fMat C = randi<fMat>(D,D,-2,2);
		do C = randi<fMat>(D,D,-2,2);
		while (std::abs(det(C))<mtol());
		
		auto tCell1 = genRandom(D).expand(C);
		const auto tCell2 = tCell1;

		tCell1.diversify();
		CPPUNIT_ASSERT(tCell1.B()==tCell2.B());
		CPPUNIT_ASSERT(tCell1.Ap()==tCell2.Ap());
		CPPUNIT_ASSERT(tCell1.types()==rg(0,1,tCell1.N()));
		CPPUNIT_ASSERT(tCell1.Ntype()==aCv(tCell1.N(),1));
		CPPUNIT_ASSERT(Ap_sorted(tCell1));

		CPPUNIT_ASSERT_EQUAL(tCell2.N(),tCell1.id().size());
		
		const auto N = tCell2.Ntype();
		size_t k=0;
		for (size_t i=0; i!=N.size(); ++i)
			for (size_t j=0; j!=N[i]; ++j,++k)
				CPPUNIT_ASSERT_EQUAL(tCell2.stripId(tCell2.id()[i]),
						     tCell1.stripId(tCell1.id()[k]));

		// check id is properly indexed
		auto uid = tCell2.id();
		std::sort(uid.begin(),uid.end());
		uid.resize(std::distance(uid.begin(),std::unique(uid.begin(),uid.end())));
		CPPUNIT_ASSERT_EQUAL(tCell2.id().size(),uid.size());
	}

	// empty ts
	{
		const aTv ts = {};
		
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;

		tCell1.diversify(ts);
		CPPUNIT_ASSERT(tCell1==tCell2);
		CPPUNIT_ASSERT(Ap_sorted(tCell1));
	}

	// with ts
	{
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;
		
		const auto ts = inds(genRndrv(tCell1.types().size()));

		tCell1.diversify(ts);
		CPPUNIT_ASSERT(tCell1.B()==tCell2.B());
		CPPUNIT_ASSERT(tCell1.Ap()==tCell2.Ap());

		const auto N = tCell1.Ntype();
		CPPUNIT_ASSERT_EQUAL(aC(tCell1.N()),std::accumulate(N.cbegin(),N.cend(),aC(0)));
		for (const auto t: ts) {
			auto i = std::find(tCell1.T_.cbegin(),tCell1.T_.cend(),tCell2.T_[t]);
			while (*i!=tCell2.T_[t+1]) {
				CPPUNIT_ASSERT_EQUAL(aT(1),*(i+1)-*i);
				++i;
			}
		}
		CPPUNIT_ASSERT(Ap_sorted(tCell1));
	}
}

void test_cell_modification::test_collectivize() {

	const auto D =  genRndST(1,5);

	auto tCell1 = genRandom(D);
	const auto tCell2 = tCell1;

	// no ts
	{
		tCell1.collectivize();
		CPPUNIT_ASSERT(tCell1.B()==tCell2.B());
		CPPUNIT_ASSERT(tCell1.types()==aTv(1,0));
		CPPUNIT_ASSERT(tCell1.Ntype()==aCv(1,tCell1.N()));
		CPPUNIT_ASSERT(Ap_sorted(tCell1));
		CPPUNIT_ASSERT_EQUAL(size_t(1),tCell1.id().size());
		CPPUNIT_ASSERT_EQUAL(tCell2.softstripId(tCell2.id().front()),
			tCell1.id().front());
	}

	// no ts, enlarged cell
	{
		fMat C = randi<fMat>(D,D,-2,2);
		do C = randi<fMat>(D,D,-2,2);
		while (std::abs(det(C))<mtol());
		
		auto tCell1 = genRandom(D).expand(C);
		const auto tCell2 = tCell1;
		
		tCell1.collectivize();
		CPPUNIT_ASSERT(tCell1.B()==tCell2.B());
		CPPUNIT_ASSERT(tCell1.types()==aTv(1,0));
		CPPUNIT_ASSERT(tCell1.Ntype()==aCv(1,tCell1.N()));
		CPPUNIT_ASSERT(Ap_sorted(tCell1));
		CPPUNIT_ASSERT_EQUAL(size_t(1),tCell1.id().size());
		CPPUNIT_ASSERT_EQUAL(tCell2.softstripId(tCell2.id().front()),
			tCell1.id().front());
	}
	
	// with ts
	{
		auto tCell1 = genRandom(D);
		const auto tCell2 = tCell1;
		
		const auto ts = inds(genRndrv(tCell1.types().size()));

		tCell1.collectivize(ts);
		CPPUNIT_ASSERT(tCell1.B()==tCell2.B());
		CPPUNIT_ASSERT_EQUAL(tCell2.N(),tCell1.N());

		size_t N=0;
		for (const auto t: ts) N+=tCell2.Ntype(t);
		fMat rAp(tCell2.dim(),0); rAp.reserve(N);
		
		if (ts.size()>1) {
			// check collectivized
			for (const auto t: ts)
				for (auto i=tCell2.ccBegin(t),e=tCell2.ccEnd(t); i!=e; ++i)
					rAp.push_back(*i);
			std::sort(rAp.cBegin(),rAp.cEnd(),vcmp);
			CPPUNIT_ASSERT(std::equal(rAp.ccBegin(),rAp.ccEnd(),tCell1.ccBegin(ts.front())));
			CPPUNIT_ASSERT_EQUAL(tCell2.id(ts.front()),tCell1.id(ts.front()));

			// check other
			const auto t1 = tCell1.types();
			const auto t2 = tCell2.types();
			CPPUNIT_ASSERT_EQUAL(t2.size(),t1.size()+ts.size()-1);

			auto tc = ts.cbegin();
			for (auto i2=t2.cbegin(),i1=t2.cbegin(),e2=t2.cend(); i2!=e2; ++i2) {
				if (*i2 == ts.front()) ++i1;
				if (tc!=ts.cend() && *i2 == *tc) ++tc;
				else {
					CPPUNIT_ASSERT_EQUAL(tCell2.Ntype(*i2),tCell1.Ntype(*i1));
					CPPUNIT_ASSERT(std::equal(tCell2.ccBegin(*i2),tCell2.ccEnd(*i2),
						tCell1.ccBegin(*i1)));
					CPPUNIT_ASSERT_EQUAL(tCell2.softstripId(tCell2.id(*i2)),
							     tCell1.softstripId(tCell1.id(*i1)));
					++i1;
				}
			}

			for (const auto& ET: tCell1.equalTypes()) {
				if (ET.size()==1)
					CPPUNIT_ASSERT_EQUAL(tCell1.softstripId(tCell1.id(ET.front())),
							tCell1.id(ET.front()));
			}
		}
	}
}

void test_cell_modification::test_merge() {
	const auto D =  genRndST(1,5);
		
	auto tCell1 = genRndCell(D);
	const auto tCell2 = tCell1;
	

	// merge with same cell has no effect
	{
		tCell1.merge(tCell1);
		CPPUNIT_ASSERT(tCell2==tCell1);
	}

	// merge with all new positions of same type
	{
		auto nAp = tCell1.Ap();
		do {
			nAp += rand<fMat>(D,tCell1.N());
			nAp %= 1.0;
		} while (std::any_of(nAp.ccBegin(),nAp.ccEnd(),[&tCell1](const auto& i)->bool{
			return std::find(tCell1.ccBegin(),tCell1.ccEnd(),i)!=tCell1.ccEnd();
		}));
		const ll_cell tCell3(tCell1.B(),std::move(nAp),tCell1.Ntype(),tCell1.id());

		tCell1.merge(tCell3);
		CPPUNIT_ASSERT_EQUAL(tCell2.N()*2,tCell1.N());
		
		for (auto i=tCell2.ccBegin(),e=tCell2.ccEnd(); i!=e; ++i)
			CPPUNIT_ASSERT(tCell1.ccEnd()!=std::find(tCell1.ccBegin(),tCell1.ccEnd(),*i));
		for (auto i=tCell3.ccBegin(),e=tCell3.ccEnd(); i!=e; ++i)
			CPPUNIT_ASSERT(tCell1.ccEnd()!=std::find(tCell1.ccBegin(),tCell1.ccEnd(),*i));
	}

	// merge unmerged with merged cell results in same cell
	{
		auto tCell3 = tCell2;
		tCell3.merge(tCell1);
		CPPUNIT_ASSERT(tCell1==tCell3);
	}
	
	// merge with new random cell
	{
		auto tCell1 = tCell2;
		auto tCell3 = genRandom(D);
		tCell3 = ll_cell(tCell1.B(),tCell3.Ap(),tCell3.Ntype(),tCell3.id());

		tCell1.merge(tCell3);
		CPPUNIT_ASSERT(tCell2.N()<=tCell1.N());

		auto id_ = tCell1.id();

		for (const auto& i: tCell2.id())
			if (tCell2.softstripId(i)!=i)
				CPPUNIT_ASSERT(std::find(id_.cbegin(),id_.cend(),i)!=id_.cend());
		for (const auto& i: tCell3.id())
			if (tCell3.softstripId(i)!=i)
				CPPUNIT_ASSERT(std::find(id_.cbegin(),id_.cend(),i)!=id_.cend());

		for (const auto& ET: tCell1.equalTypes()) {
			if (ET.size()<2) continue;
			CPPUNIT_ASSERT(std::all_of(ET.cbegin(),ET.cend(),[&tCell1](const aT t)->bool
				{ return tCell1.softstripId(tCell1.id(t))!=tCell1.id(t); }));
		}
	}
}


const char* test_cell_modification::test_id() noexcept {
	return "test_cell_modification";
}

CppUnit::Test* test_cell_modification::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());

	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_stress", &test_cell_modification::test_stress));
	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_swapDim", &test_cell_modification::test_swapDim));
	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_invDim", &test_cell_modification::test_invDim));
	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_orient", &test_cell_modification::test_orient));
	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_permute", &test_cell_modification::test_permute));
	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_rotate", &test_cell_modification::test_rotate));
	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_scale", &test_cell_modification::test_scale));
	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_changeBasis", &test_cell_modification::test_changeBasis));
	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_makePrimitive", &test_cell_modification::test_makePrimitive));
	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_shift", &test_cell_modification::test_shift));
	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_diversify", &test_cell_modification::test_diversify));
	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_collectivize", &test_cell_modification::test_collectivize));
	suite->addTest(new CppUnit::TestCaller<test_cell_modification>(
		"test_merge", &test_cell_modification::test_merge));
	
	return suite;
}
