// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_io.h"
#include "ll_io.h"
#include "ll_fn.h"
#include "ll_cell.h"
#include "ll_testTools.h"
#include "lm_testTools.h"
#include "aux_io.h"
#include <random>
#include <chrono>


using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;
using namespace aux;

fMat test_io::genRndP_(const size_t d) const noexcept {
	std::vector<size_t> I(d);
	std::iota(I.begin(),I.end(),0);
	std::shuffle(I.begin(), I.end(), std::default_random_engine(
		std::chrono::system_clock::now().time_since_epoch().count()));

	fMat res(d,0); res.reserve(d);
	for (const auto i: I)
		res.push_back(cId<fMat>(d,i));
	return res;
}


void test_io::test_writePOSCAR() {
	const size_t D = genRndST(1,5);
	const std::string prefix = "outp/";

	// check exceptions
	{
		// B not empty
		CPPUNIT_ASSERT_THROW(printPOSCAR("bla",fMat(),fMat(),aTv()),std::invalid_argument);

		// B not square
		CPPUNIT_ASSERT_THROW(printPOSCAR("bla",fMat(3,2),fMat(3,4),aTv()),std::invalid_argument);

		// B.M()!=Ap.M()
		CPPUNIT_ASSERT_THROW(printPOSCAR("bla",fMat(3),fMat(4,5),aTv()),std::invalid_argument);

		// !T.empty() && T.size()!=Ap.N()
		CPPUNIT_ASSERT_THROW(printPOSCAR("bla",fMat(3),fMat(3,4),{0,1},{},1.0,true),
				std::invalid_argument);

		// legal input no throw
		CPPUNIT_ASSERT_NO_THROW(printPOSCAR(prefix+"POSCAR_legal",fMat(3),fMat(3,4),
					{0,0,1,1},{"S","Mo"},1.0,true));
	}
	
	// random cell, direct
	{
		const auto tCell1 = genRandom(D);
		printPOSCAR(prefix+"POSCAR1",tCell1,true,false,16,"random_POSCAR_test_direct");
		ll_cell tCell2(prefix+"POSCAR1",D);

		CPPUNIT_ASSERT(tCell1==tCell2);
	}
	// random cell, cartesian
	{
		const auto tCell1 = genRandom(D);
		printPOSCAR(prefix+"POSCAR2",tCell1,false,false,16,"random_POSCAR_test_cart");
		ll_cell tCell2(prefix+"POSCAR2",D);

		CPPUNIT_ASSERT(tCell1==tCell2);
	}
	// random cell, cell input
	{
		const auto tCell1 = genRandom(D);
		printPOSCAR(prefix+"POSCAR3",tCell1,true,false,16,"random_POSCAR_test_cell");
		ll_cell tCell2(prefix+"POSCAR3",D);

		CPPUNIT_ASSERT(tCell1==tCell2);
	}
	// zincblende
	{
		const auto tCell1 = genZincblende(3,1.0);
		printPOSCAR(prefix+"POSCAR4",tCell1,true,false,16,"GaAs_test");
		ll_cell tCell2(prefix+"POSCAR4",3);

		CPPUNIT_ASSERT(tCell1==tCell2);
	}
}

void test_io::test_readEf() {
	const double ref = -1.5338;
	const double Ef = readEf("data/dft/OUTCAR");
	CPPUNIT_ASSERT_DELTA(ref,Ef,1e-4);

	// full vs minimal
	{
		const std::string path = "data/w90/";
		const auto Ef_full = readEf(path+"mos2/OUTCAR");
		const auto Ef_min = readEf(path+"minimal/OUTCAR");
		CPPUNIT_ASSERT_EQUAL(Ef_full,Ef_min);
	}
}

void test_io::test_genr() {
	// d = 3
	{
		const auto r = genr(3);
		CPPUNIT_ASSERT(rv({false,false,false})==r);
	}
	
	// d = 2
	{
		const auto r = genr(2);
		CPPUNIT_ASSERT(rv({false,true,false})==r);
	}
	
	// d = 1
	{
		const auto r = genr(1);
		CPPUNIT_ASSERT(rv({false,true,true})==r);
	}
}

void test_io::test_readLayerMatrix() {
	const auto LM = readLayerMatrix("data/omen/Layer_Matrix.dat");

	const auto ckmat = fMat("data/omen/lmck").T();
	const auto cktype = fMat("data/omen/lmt");

	CPPUNIT_ASSERT(LM.Ap()==ckmat);
	CPPUNIT_ASSERT(fMat(LM.T())==cktype);
}

void test_io::test_printOmf() {
	// check exceptions
	{
		CPPUNIT_ASSERT_THROW(printOmf("bla",{-1.0,1.0},{},{1.3}),std::invalid_argument);
		CPPUNIT_ASSERT_THROW(printOmf("bla",{-1.0,1.0},{2,3},{1.3}),std::invalid_argument);
	}

	// write test file
	{
		const std::vector<size_t> Norb = {3,4,3,3,4,3};
		const std::vector<double> id = {1.0,2.0,3.0,4.0,5.0,6.0};
		const std::string fileName = "outp/test_mat_par";

		printOmf(fileName,{-1.0,1.0},Norb,id);
	}
}

void test_io::test_printOlf() {
	ll_cell tCell("data/dft/POSCAR1");
	tCell.changeBasis(findBasis(tCell.B(),eye<fMat>(3)));
	const auto NN = genNNmat({false,false,false});

	idv id; id.reserve(tCell.N());
	for (const auto t: tCell.type())
		id.push_back("A"+std::to_string(t+1));
	
	// check exceptions
	{
		CPPUNIT_ASSERT_THROW(printOlf("bla",{},tCell.B().prod(tCell.Ap()),
					id,6,1.0),std::invalid_argument);
		CPPUNIT_ASSERT_THROW(printOlf("bla",zeros<fMat>(3,2),tCell.B().prod(tCell.Ap()),
					id,6,1.0),std::invalid_argument);
		CPPUNIT_ASSERT_THROW(printOlf("bla",zeros<fMat>(4),tCell.B().prod(tCell.Ap()),
					id,6,1.0),std::invalid_argument);
		CPPUNIT_ASSERT_THROW(printOlf("bla",tCell.B(),tCell.B().prod(tCell.Ap()),
					{"A1"},6,1.0),std::invalid_argument);
	}

	// write test file
	{
		const double f = 1.1;

		const auto b = tCell.getBonds(f,NN);
		const auto bl = b.radius();
		
		const auto Nb = b.Nindex(b.inds());
		const auto nn = *std::max_element(Nb.cbegin(),Nb.cend())/2;

		const std::string fileName = "outp/test_lattice_dat";

		printOlf(fileName,tCell.B(),tCell.B().prod(tCell.Ap()),
				id,nn,bl);
	}
}

void test_io::test_readOmf() {
	const std::string path = "data/omen/";

	CPPUNIT_ASSERT_THROW(readOmf("dsgjdkas"),std::invalid_argument);
	CPPUNIT_ASSERT_NO_THROW(readOmf(path+"mat_par"));

	const auto inp = readOmf(path+"mat_par");
	CPPUNIT_ASSERT_EQUAL(-.6236,inp.E.cb);
	CPPUNIT_ASSERT_EQUAL(-1.7451,inp.E.vb);
	CPPUNIT_ASSERT(std::vector<size_t>({3,5,3,3,5,3})==inp.Norb);
	CPPUNIT_ASSERT(std::vector<double>({1.1,2.1,3.1,3.4,5.4,1.0})==inp.mass);
}

void test_io::test_readOlf() {
	const std::string path = "data/omen/";

	CPPUNIT_ASSERT_THROW(readOlf("dsgjdkas"),std::invalid_argument);
	CPPUNIT_ASSERT_NO_THROW(readOlf(path+"lattice_dat"));
	
	const auto inp = readOlf(path+"lattice_dat");
	CPPUNIT_ASSERT(!inp.empty());
	CPPUNIT_ASSERT_EQUAL(size_t(3),inp.dim());
	CPPUNIT_ASSERT_EQUAL(size_t(144),inp.N());

	const fMat B(path+"ref/B.mat");
	CPPUNIT_ASSERT(T(B)==inp.B);

	const fMat Ap(path+"ref/Ap.mat");
	CPPUNIT_ASSERT(T(Ap).cAt(0)==inp.Ap.cAt(0));
	CPPUNIT_ASSERT(T(Ap).cAt(1)==inp.Ap.cAt(1));
	CPPUNIT_ASSERT(T(Ap).cAt(2)==inp.Ap.cAt(2));
	CPPUNIT_ASSERT(T(Ap).cAt(3)==inp.Ap.cAt(141));
	CPPUNIT_ASSERT(T(Ap).cAt(4)==inp.Ap.cAt(142));
	CPPUNIT_ASSERT(T(Ap).cAt(5)==inp.Ap.cAt(143));

	CPPUNIT_ASSERT_EQUAL(std::string("A6"),inp.id[0]);
	CPPUNIT_ASSERT_EQUAL(std::string("A6"),inp.id[1]);
	CPPUNIT_ASSERT_EQUAL(std::string("A6"),inp.id[2]);
	CPPUNIT_ASSERT_EQUAL(std::string("A5"),inp.id[141]);
	CPPUNIT_ASSERT_EQUAL(std::string("A5"),inp.id[142]);
	CPPUNIT_ASSERT_EQUAL(std::string("A5"),inp.id[143]);

	CPPUNIT_ASSERT_EQUAL(size_t(4),inp.nn);
	CPPUNIT_ASSERT_EQUAL(2.58634,inp.bl);

	CPPUNIT_ASSERT_EQUAL(inp.id.size(),inp.T.size());

	idv uid; uid.reserve(6);
	for (const auto& s: inp.id)
		if (std::find(uid.cbegin(),uid.cend(),s)==uid.cend())
			uid.push_back(s);
	for (size_t i=0; i!=inp.id.size(); ++i)
		CPPUNIT_ASSERT_EQUAL(inp.id[i],uid[inp.T[i]]);
}

void test_io::test_hrDim() {
	const std::string wpath = "data/w90/";

	// mos2
	{
		const auto r = hrDim(wpath+"mos2/wannier90_hr.dat");
		rv ref = {false,true,false};
		CPPUNIT_ASSERT(ref==r);
	}

	// inas mono
	{
		const auto r = hrDim(wpath+"inas_mono/wannier90_hr.dat");
		rv ref = {false,false,true};
		CPPUNIT_ASSERT(ref==r);
	}
}

void test_io::test_readB() {
	const std::string wpath = "data/w90/";
	const std::string opath = "data/omen/";

	// check exceptions
	{
		// check bad file
		CPPUNIT_ASSERT_THROW(readB(opath+"Layer_Matrix.dat"), std::runtime_error);
			
		// check non existant file
		CPPUNIT_ASSERT_THROW(readB("blub.dat"), std::invalid_argument);
			
		// check good file no throw
		CPPUNIT_ASSERT_NO_THROW(readB(wpath+"wannier90.wout"));
	}

	// check data is read correctly
	{
		// without P
		{
			const auto B = readB(wpath+"wannier90.wout");
			const fMat ck({3.4,	0.0,		0.0,
				       0.0,	29.893856,	0.0,
				      -1.7,	0.0,		2.944486},3,3);
			CPPUNIT_ASSERT(B==ck);
		}
	}

	// check minimal vs full
	{
		const auto Bfull = readB(wpath+"mos2/wannier90.wout");
		const auto Bmin = readB(wpath+"minimal/wannier90.wout");
		CPPUNIT_ASSERT(Bfull==Bmin);
	}
}

void test_io::test_readAp() {
	const std::string wpath = "data/w90/";
	const std::string opath = "data/omen/";

	// check exceptions
	{
		// check bad file
		CPPUNIT_ASSERT_THROW(readAp(opath+"Layer_Matrix.dat"), std::runtime_error);
			
		// check non existant file
		CPPUNIT_ASSERT_THROW(readAp("blub.dat"), std::invalid_argument);
	
		// check abl not sorted
		CPPUNIT_ASSERT_THROW(readAp(wpath+"wannier90.wout",true,
					{3,4,1,0}), std::invalid_argument);

		// check abl duplicate entries
		CPPUNIT_ASSERT_THROW(readAp(wpath+"wannier90.wout",true,
					{0,1,1,2}), std::invalid_argument);
			
		// check abl illegal entries
		CPPUNIT_ASSERT_THROW(readAp(wpath+"wannier90.wout",true,
					{0,1,2,172645}), std::invalid_argument);
	
		// check good file no throw
		CPPUNIT_ASSERT_NO_THROW(readAp(wpath+"wannier90.wout"));
	}

	// direct
	{
		const auto L = readAp(wpath+"wannier90.wout",true);
		const fMat ck({0.0,		0.05620,	0.0,
			       0.33333,		0.27477,	0.66667,
			       0.33333,		0.0,		0.66667,
			       0.0,		0.33097,	0.0,
			       0.0,		0.21857,	0.0,
			       0.33333,		0.11240,	0.66667},3,6);
		const idv ckid = {"Mo","Mo","S","S","S","S"};

		CPPUNIT_ASSERT(L.Ap()==ck);
		CPPUNIT_ASSERT(L.id()==ckid);
	}
	// direct, with abl
	{
		const auto L = readAp(wpath+"wannier90.wout",true,{1,3});
		const fMat ck({0.0,		0.05620,	0.0,
			       0.33333,		0.0,		0.66667,
			       0.0,		0.21857,	0.0,
			       0.33333,		0.11240,	0.66667},3,4);
		const idv ckid = {"Mo","S","S","S"};


		CPPUNIT_ASSERT(L.Ap()==ck);
		CPPUNIT_ASSERT(L.id()==ckid);
	}
	// direct, with abl, one type missing
	{
		const auto L = readAp(wpath+"wannier90.wout",true,{0,1});
		const fMat ck({0.33333,		0.0,		0.66667,
			       0.0,		0.33097,	0.0,
			       0.0,		0.21857,	0.0,
			       0.33333,		0.11240,	0.66667},3,4);
		const idv ckid = {"S","S","S","S"};

		CPPUNIT_ASSERT(L.Ap()==ck);
		CPPUNIT_ASSERT(L.id()==ckid);
	}

	// cartesian
	{
		const auto L = readAp(wpath+"wannier90.wout",false);
		const fMat ck({0.0,	1.68,		0.0,
			       0.0,	8.21385,	1.96299,
			       0.0,	0.0,		1.96299,
			       0.0,	9.89385,	0.0,
			       0.0,	6.53387,	0.0,
			       0.0,	3.36001,	1.96299},3,6);
		const idv ckid = {"Mo","Mo","S","S","S","S"};

		CPPUNIT_ASSERT(L.Ap()==ck);
		CPPUNIT_ASSERT(L.id()==ckid);
	}
	// cartesian, with abl;
	{
		const auto L = readAp(wpath+"wannier90.wout",false,{1,3});
		const fMat ck({0.0,	1.68,		0.0,
			       0.0,	0.0,		1.96299,
			       0.0,	6.53387,	0.0,
			       0.0,	3.36001,	1.96299},3,4);
		const idv ckid = {"Mo","S","S","S"};

		CPPUNIT_ASSERT(L.Ap()==ck);
		CPPUNIT_ASSERT(L.id()==ckid);
	}
	// cartesian, with abl one type missing
	{
		const auto L = readAp(wpath+"wannier90.wout",false,{0,1});
		const fMat ck({0.0,	0.0,		1.96299,
			       0.0,	9.89385,	0.0,
			       0.0,	6.53387,	0.0,
			       0.0,	3.36001,	1.96299},3,4);
		const idv ckid = {"S","S","S","S"};


		CPPUNIT_ASSERT(L.Ap()==ck);
		CPPUNIT_ASSERT(L.id()==ckid);
	}

	// compare direct and cartesian
	{
		const std::vector<size_t> abl = {0,3,5};

		const auto B = readB(wpath+"wannier90.wout");
		const auto Lc = readAp(wpath+"wannier90.wout",false,abl);
		const auto Lf = readAp(wpath+"wannier90.wout",true,abl);
		
		set_mtol(1e-4);	// shitty tolerance due to wannier90 output accuracy
		CPPUNIT_ASSERT(B.prod(Lf.Ap())==Lc.Ap());
		reset_mtol();
		CPPUNIT_ASSERT(Lf.id()==Lc.id());
	}
	
	// check minimal vs full
	{
		// cartesian
		{
			const auto Apfull = readAp(wpath+"mos2/wannier90.wout",false);
			const auto Apmin = readAp(wpath+"minimal/wannier90.wout",false);
			CPPUNIT_ASSERT(Apfull.Ap()==Apmin.Ap());
			CPPUNIT_ASSERT(Apfull.id()==Apmin.id());
		}
		
		// direct
		{
			const auto Apfull = readAp(wpath+"mos2/wannier90.wout",true);
			const auto Apmin = readAp(wpath+"minimal/wannier90.wout",true);
			CPPUNIT_ASSERT(Apfull.Ap()==Apmin.Ap());
			CPPUNIT_ASSERT(Apfull.id()==Apmin.id());
		}
	}
}

void test_io::test_readWp() {
	const std::string wpath = "data/w90/";
	const std::string opath = "data/omen/";
	
	// check exceptions
	{
		// check bad file
		CPPUNIT_ASSERT_THROW(readWp(opath+"Layer_Matrix.dat"), std::runtime_error);
			
		// check non existant file
		CPPUNIT_ASSERT_THROW(readWp("blub.dat"), std::invalid_argument);
	
		// check wbl not sorted
		CPPUNIT_ASSERT_THROW(readWp(wpath+"wannier90.wout",
					{3,4,1,0}), std::invalid_argument);

		// check wbl duplicate entries
		CPPUNIT_ASSERT_THROW(readWp(wpath+"wannier90.wout",
					{0,1,1,2}), std::invalid_argument);
			
		// check wbl illegal entries
		CPPUNIT_ASSERT_THROW(readWp(wpath+"wannier90.wout",
					{0,1,2,172645}), std::invalid_argument);

		// check good file no throw
		CPPUNIT_ASSERT_NO_THROW(readWp(wpath+"wannier90.wout"));
	}
	
	// no wbl, no P
	{
		const auto Wp = readWp(wpath+"wannier90.wout");
	
		const fMat ck({
		0.482391,  7.078464, -0.007794,
		0.000192, -0.423313,  1.965239,
		-0.007638, 10.308024, -0.023286,
		-0.553955,  2.917203,  2.264988,
		0.038522,  2.946627,  1.315318,
		0.003413,  9.467770,  0.644471,
		-1.139898,  0.450749, -0.680885,
		-0.016440,  0.432696,  1.340404,
		0.551736,  9.443958, -0.298955,
		-1.148999,  2.984042, -0.638365,
		0.113523,  6.162493, -0.105254,
		3.475106,  1.740480,  0.023430,
		-0.538264,  9.432773, -0.297512,
		-0.014971,  1.735985, -0.031433,
		-0.552809,  6.967889, -0.346652,
		-0.033350,  3.788433,  1.960502,
		-0.007120,  8.157331,  1.957693,
		-0.031401,  8.083511,  2.019429,
		1.156762,  0.446536, -0.657528,
		-0.057934,  1.731305,  0.006021,
		-1.700063,  8.081807, -0.995361,
		-0.061165,  6.746088,  0.445236},3,22);
		
		const fMat ckspread({
		1.98218953,1.36059869,1.41498919,1.92502688,1.98677354,
		1.89215823,1.86081759,1.94003709,1.86943488,1.98311332,
		1.59689773,2.37554007,1.87489639,2.14892524,1.91210902,
		1.51071199,2.14108509,3.73373851,1.86329425,2.18784727,
		3.56399432,2.31241476},1,22);
		
		CPPUNIT_ASSERT(Wp.Wp()==ck);
		CPPUNIT_ASSERT(Wp.s()==ckspread);
	}
		
	// with wbl
	{
		const auto Wp = readWp(wpath+"wannier90.wout",{0,5,13});
	
		const fMat ck({
		0.000192, -0.423313,  1.965239,
		-0.007638, 10.308024, -0.023286,
		-0.553955,  2.917203,  2.264988,
		0.038522,  2.946627,  1.315318,
		-1.139898,  0.450749, -0.680885,
		-0.016440,  0.432696,  1.340404,
		0.551736,  9.443958, -0.298955,
		-1.148999,  2.984042, -0.638365,
		0.113523,  6.162493, -0.105254,
		3.475106,  1.740480,  0.023430,
		-0.538264,  9.432773, -0.297512,
		-0.552809,  6.967889, -0.346652,
		-0.033350,  3.788433,  1.960502,
		-0.007120,  8.157331,  1.957693,
		-0.031401,  8.083511,  2.019429,
		1.156762,  0.446536, -0.657528,
		-0.057934,  1.731305,  0.006021,
		-1.700063,  8.081807, -0.995361,
		-0.061165,  6.746088,  0.445236},3,19);
		
		const fMat ckspread({
		1.36059869,1.41498919,1.92502688,1.98677354,
		1.86081759,1.94003709,1.86943488,1.98311332,
		1.59689773,2.37554007,1.87489639,1.91210902,
		1.51071199,2.14108509,3.73373851,1.86329425,2.18784727,
		3.56399432,2.31241476},1,19);
		
		CPPUNIT_ASSERT(Wp.Wp()==ck);
		CPPUNIT_ASSERT(Wp.s()==ckspread);
	}

	// cartesian
	{
		const auto Wpfull = readWp(wpath+"mos2/wannier90.wout");
		const auto Wpmin = readWp(wpath+"minimal/wannier90.wout");
		CPPUNIT_ASSERT(Wpfull.Wp()==Wpmin.Wp());
		CPPUNIT_ASSERT(Wpfull.s()==Wpmin.s());
	}
}

void test_io::test_readHr() {
	const std::string wpath = "data/w90/";
	
	// check exceptions
	{
		// check bad file
		CPPUNIT_ASSERT_THROW(readHr(wpath+"wannier90_hr.dat.cut"), std::runtime_error);
			
		// check non existant file
		CPPUNIT_ASSERT_THROW(readHr("blub.dat"), std::invalid_argument);
	
		// check wbl not sorted
		CPPUNIT_ASSERT_THROW(readHr(wpath+"wannier90_hr.dat",
					{3,4,1,0}), std::invalid_argument);

		// check wbl duplicate entries
		CPPUNIT_ASSERT_THROW(readHr(wpath+"wannier90_hr.dat",
					{0,1,1,2}), std::invalid_argument);
			
		// check wbl illegal entries
		CPPUNIT_ASSERT_THROW(readHr(wpath+"wannier90_hr.dat",
					{0,1,2,172645}), std::invalid_argument);

		// check good file no throw
		CPPUNIT_ASSERT_NO_THROW(readHr(wpath+"wannier90_hr.dat"));
	}
	
	// check reading with tol==0.0, no wbl
	{
		const auto H = readHr(wpath+"wannier90_hr.dat");
		const size_t NR = 39;
		const size_t Nw = 22;

		CPPUNIT_ASSERT_EQUAL(size_t(3),H.dim());
		CPPUNIT_ASSERT_EQUAL(NR,H.N());
		CPPUNIT_ASSERT(std::all_of(H.begin(),H.end(),[](const cMat& i){
			return i.M()==Nw && i.N()==Nw; }));
	
		CPPUNIT_ASSERT(H.cAt(NR/2)==zeros<fMat>(3,1));
		for (size_t i=0; i!=NR/2; ++i)
			CPPUNIT_ASSERT(H.cAt(i)==(-H.cAt(NR-1-i)));
	
		CPPUNIT_ASSERT(H[NR/2]==lm__::T(H[NR/2]));
		for (size_t i=0; i!=NR/2; ++i)
			CPPUNIT_ASSERT(H[i]==lm__::T(H[NR-1-i]));
		
		// compare to reading with wbl;
		{
			const std::vector<size_t> wbl = genRndIvec(genRndST(0,Nw),0,Nw-1);
				
			const auto Hwbl = readHr(wpath+"wannier90_hr.dat",wbl);
				
			CPPUNIT_ASSERT_EQUAL(H.dim(),Hwbl.dim());
			CPPUNIT_ASSERT_EQUAL(H.N(),Hwbl.N());
			CPPUNIT_ASSERT(H.R()==Hwbl.R());

			auto Htmp = H;
			auto ih = Htmp.begin();
			auto iwbl = Hwbl.begin();
			for (; ih!=Htmp.end(); ++ih, ++iwbl)
				CPPUNIT_ASSERT(*iwbl==(ih->cRm(wbl).rRm(wbl)));
		}

		// compare to reading with tol==0.1
		{
			const double tol=0.1;
			const auto Htol = readHr(wpath+"wannier90_hr.dat",{},
					[tol](const cMat& inp)->bool{
						for (const auto& i: inp)
							if (std::real(i)>=tol) return true;
						return false;
					});
				
			const size_t NR = 17;
			const fMat Rck({-2,0,-2,
					-2,0,-1,
					-2,0,0,
					-1,0,-2,
					-1,0,-1,
					-1,0,0,
					-1,0,1,
				   	 0,0,-1,
					 0,0,0,
					 0,0,1,
					 1,0,-1,
					 1,0,0,
					 1,0,1,
					 1,0,2,
					 2,0,0,
					 2,0,1,
					 2,0,2},3,NR);
			CPPUNIT_ASSERT(Rck==Htol.R());
			CPPUNIT_ASSERT_EQUAL(NR,Htol.size());
				
			for (size_t i=0; i<NR; ++i) {
				const auto j = std::find(H.ccBegin(),H.ccEnd(),Htol.cAt(i));
				CPPUNIT_ASSERT(H[size_t(j)]==Htol[i]);
			}
		}
	}

	// compare full vs minimal
	{
		const std::string path = "data/w90/";
		const auto Hfull = readHr(wpath+"mos2/wannier90_hr.dat");
		const auto Hmin = readHr(wpath+"minimal/wannier90_hr.dat");
		CPPUNIT_ASSERT(Hfull.R()==Hmin.R());
		CPPUNIT_ASSERT(Hfull.H()==Hmin.H());
	}
}

void test_io::test_readAMN() {
	const std::string wpath = "data/w90/mos2_dl/";

	const auto amn = readAMN(wpath+"wannier90.amn");
	CPPUNIT_ASSERT_EQUAL(size_t(49),amn.size());
	CPPUNIT_ASSERT(std::all_of(amn.cbegin(),amn.cend(),
		[](const auto& i)->bool{return i.M()==11 && i.N()==11;}));

	std::complex<double> ref1(-0.128058715695,0.149159827634);
	CPPUNIT_ASSERT_EQUAL(ref1,amn[21-1](8-1,10-1));

	std::complex<double> ref2(-0.311164198266,-0.042322857367);
	CPPUNIT_ASSERT_EQUAL(ref2,amn[36-1](2-1,9-1));
}

void test_io::test_readXYZ() {
	const std::string wpath = "data/w90/mos2_dl/";
	
	// check exceptions
	{
		// check bad file
		CPPUNIT_ASSERT_THROW(readXYZ(wpath+"wannier90_r.dat.cut"), std::runtime_error);
			
		// check non existant file
		CPPUNIT_ASSERT_THROW(readXYZ("blub.dat"), std::invalid_argument);
	
		// check wbl not sorted
		CPPUNIT_ASSERT_THROW(readXYZ(wpath+"wannier90_r.dat",0,
					{3,4,1,0}), std::invalid_argument);

		// check wbl duplicate entries
		CPPUNIT_ASSERT_THROW(readXYZ(wpath+"wannier90_r.dat",0,
					{0,1,1,2}), std::invalid_argument);
			
		// check wbl illegal entries
		CPPUNIT_ASSERT_THROW(readXYZ(wpath+"wannier90_r.dat",0,
					{0,1,2,172645}), std::invalid_argument);

		// check good file no throw
		CPPUNIT_ASSERT_NO_THROW(readXYZ(wpath+"wannier90_r.dat",0));
		CPPUNIT_ASSERT_NO_THROW(readXYZ(wpath+"wannier90_r.dat",1));
		CPPUNIT_ASSERT_NO_THROW(readXYZ(wpath+"wannier90_r.dat",2));
		CPPUNIT_ASSERT_NO_THROW(readXYZ(wpath+"wannier90_r.dat",3));
		CPPUNIT_ASSERT(readXYZ(wpath+"wannier90_r.dat",3).empty());
	}
	
	// check reading with tol==0.0, no wbl
	{
		const auto X = readX(wpath+"wannier90_r.dat");
		const size_t NR = 69;
		const size_t Nw = 22;

		CPPUNIT_ASSERT(X.dim()==3);
		CPPUNIT_ASSERT(X.N()==NR);
		CPPUNIT_ASSERT(std::all_of(X.begin(),X.end(),[](const cMat& i){
			return i.M()==Nw && i.N()==Nw; }));
		
		// R vectors symmetric
		CPPUNIT_ASSERT(X.cAt(NR/2)==zeros<fMat>(3,1));
		for (size_t i=0; i!=NR/2; ++i)
			CPPUNIT_ASSERT(X.cAt(i)==(-X.cAt(NR-1-i)));
	
		// compare to reading with wbl;
		{
			const std::vector<size_t> wbl = genRndIvec(genRndST(0,Nw),0,Nw-1);
				
			const auto Xwbl = readX(wpath+"wannier90_r.dat",wbl);
				
			CPPUNIT_ASSERT_EQUAL(X.dim(),Xwbl.dim());
			CPPUNIT_ASSERT_EQUAL(X.N(),Xwbl.N());
			CPPUNIT_ASSERT_EQUAL(X.R(),Xwbl.R());
		
			// R vectors symmetric
			CPPUNIT_ASSERT(Xwbl.cAt(NR/2)==zeros<fMat>(3,1));
			for (size_t i=0; i!=NR/2; ++i)
				CPPUNIT_ASSERT(Xwbl.cAt(i)==(-Xwbl.cAt(NR-1-i)));

			auto Xtmp = X;
			auto ih = Xtmp.begin();
			auto iwbl = Xwbl.begin();
			for (; ih!=Xtmp.end(); ++ih, ++iwbl)
				CPPUNIT_ASSERT(*iwbl==(ih->cRm(wbl).rRm(wbl)));
		}

		// compare to reading with tol
		{
			const double tol = genRndDouble(1e-3,1e-1);

			const auto Xtol = readX(wpath+"wannier90_r.dat",{},
					[tol](const cMat& inp)->bool{
						for (const auto& i: inp)
							if (std::real(i)>=tol) return true;
						return false;
					});
			const size_t NR = Xtol.N();
			
			// R vectors symmetric
			CPPUNIT_ASSERT(Xtol.cAt(NR/2)==zeros<fMat>(3,1));
			for (size_t i=0; i!=NR/2; ++i)
				CPPUNIT_ASSERT(Xtol.cAt(i)==(-Xtol.cAt(NR-1-i)));
			
			CPPUNIT_ASSERT_EQUAL(X.dim(),Xtol.dim());
			CPPUNIT_ASSERT(X.N()>=Xtol.N());

			// check what is there is the same
			for (auto i=Xtol.ccBegin(),ie=Xtol.ccEnd(); i!=ie; ++i) {
				const auto itr = std::lower_bound(X.ccBegin(),X.ccEnd(),*i,vcmp);
				CPPUNIT_ASSERT(X.ccEnd()!=itr);
				CPPUNIT_ASSERT_EQUAL(X[size_t(itr)],Xtol[size_t(i)]);
			}
		}
	}
	{
		const auto Y = readY(wpath+"wannier90_r.dat");
		const size_t NR = 69;
		const size_t Nw = 22;

		CPPUNIT_ASSERT(Y.dim()==3);
		CPPUNIT_ASSERT(Y.N()==NR);
		CPPUNIT_ASSERT(std::all_of(Y.begin(),Y.end(),[](const cMat& i){
			return i.M()==Nw && i.N()==Nw; }));
	
		// R vectors symmetric
		CPPUNIT_ASSERT(Y.cAt(NR/2)==zeros<fMat>(3,1));
		for (size_t i=0; i!=NR/2; ++i)
			CPPUNIT_ASSERT(Y.cAt(i)==(-Y.cAt(NR-1-i)));
		
		// compare to reading with wbl;
		{
			const std::vector<size_t> wbl = genRndIvec(genRndST(0,Nw),0,Nw-1);
				
			const auto Ywbl = readY(wpath+"wannier90_r.dat",wbl);
			
			// R vectors symmetric
			CPPUNIT_ASSERT(Ywbl.cAt(NR/2)==zeros<fMat>(3,1));
			for (size_t i=0; i!=NR/2; ++i)
				CPPUNIT_ASSERT(Ywbl.cAt(i)==(-Ywbl.cAt(NR-1-i)));
				
			CPPUNIT_ASSERT_EQUAL(Y.dim(),Ywbl.dim());
			CPPUNIT_ASSERT_EQUAL(Y.N(),Ywbl.N());
			CPPUNIT_ASSERT_EQUAL(Y.R(),Ywbl.R());

			auto Ytmp = Y;
			auto ih = Ytmp.begin();
			auto iwbl = Ywbl.begin();
			for (; ih!=Ytmp.end(); ++ih, ++iwbl)
				CPPUNIT_ASSERT(*iwbl==(ih->cRm(wbl).rRm(wbl)));
		}
		
		// compare to reading with tol
		{
			const double tol = genRndDouble(1e-3,1e-1);

			const auto Ytol = readY(wpath+"wannier90_r.dat",{},
					[tol](const cMat& inp)->bool{
						for (const auto& i: inp)
							if (std::real(i)>=tol) return true;
						return false;
					});
			const size_t NR = Ytol.N();
			
			// R vectors symmetric
			CPPUNIT_ASSERT(Ytol.cAt(NR/2)==zeros<fMat>(3,1));
			for (size_t i=0; i!=NR/2; ++i)
				CPPUNIT_ASSERT(Ytol.cAt(i)==(-Ytol.cAt(NR-1-i)));
			
			CPPUNIT_ASSERT_EQUAL(Y.dim(),Ytol.dim());
			CPPUNIT_ASSERT(Y.N()>=Ytol.N());

			// check what is there is the same
			for (auto i=Ytol.ccBegin(),ie=Ytol.ccEnd(); i!=ie; ++i) {
				const auto itr = std::lower_bound(Y.ccBegin(),Y.ccEnd(),*i,vcmp);
				CPPUNIT_ASSERT(Y.ccEnd()!=itr);
				CPPUNIT_ASSERT_EQUAL(Y[size_t(itr)],Ytol[size_t(i)]);
			}
		}
	}
	{
		const auto Z = readZ(wpath+"wannier90_r.dat");
		const size_t NR = 69;
		const size_t Nw = 22;

		CPPUNIT_ASSERT(Z.dim()==3);
		CPPUNIT_ASSERT(Z.N()==NR);
		CPPUNIT_ASSERT(std::all_of(Z.begin(),Z.end(),[](const cMat& i){
			return i.M()==Nw && i.N()==Nw; }));
	
		// R vectors symmetric
		CPPUNIT_ASSERT(Z.cAt(NR/2)==zeros<fMat>(3,1));
		for (size_t i=0; i!=NR/2; ++i)
			CPPUNIT_ASSERT(Z.cAt(i)==(-Z.cAt(NR-1-i)));
		
		// compare to reading with wbl;
		{
			const std::vector<size_t> wbl = genRndIvec(genRndST(0,Nw),0,Nw-1);
				
			const auto Zwbl = readZ(wpath+"wannier90_r.dat",wbl);
			
			// R vectors symmetric
			CPPUNIT_ASSERT(Zwbl.cAt(NR/2)==zeros<fMat>(3,1));
			for (size_t i=0; i!=NR/2; ++i)
				CPPUNIT_ASSERT(Zwbl.cAt(i)==(-Zwbl.cAt(NR-1-i)));
				
			CPPUNIT_ASSERT_EQUAL(Z.dim(),Zwbl.dim());
			CPPUNIT_ASSERT_EQUAL(Z.N(),Zwbl.N());
			CPPUNIT_ASSERT_EQUAL(Z.R(),Zwbl.R());

			auto Ztmp = Z;
			auto ih = Ztmp.begin();
			auto iwbl = Zwbl.begin();
			for (; ih!=Ztmp.end(); ++ih, ++iwbl)
				CPPUNIT_ASSERT(*iwbl==(ih->cRm(wbl).rRm(wbl)));
		}
		// compare to reading with tol
		{
			const double tol = genRndDouble(1e-3,1e-1);

			const auto Ztol = readZ(wpath+"wannier90_r.dat",{},
					[tol](const cMat& inp)->bool{
						for (const auto& i: inp)
							if (std::real(i)>=tol) return true;
						return false;
					});
			const size_t NR = Ztol.N();
			
			// R vectors symmetric
			CPPUNIT_ASSERT(Ztol.cAt(NR/2)==zeros<fMat>(3,1));
			for (size_t i=0; i!=NR/2; ++i)
				CPPUNIT_ASSERT(Ztol.cAt(i)==(-Ztol.cAt(NR-1-i)));
			
			CPPUNIT_ASSERT_EQUAL(Z.dim(),Ztol.dim());
			CPPUNIT_ASSERT(Z.N()>=Ztol.N());

			// check what is there is the same
			for (auto i=Ztol.ccBegin(),ie=Ztol.ccEnd(); i!=ie; ++i) {
				const auto itr = std::lower_bound(Z.ccBegin(),Z.ccEnd(),*i,vcmp);
				CPPUNIT_ASSERT(Z.ccEnd()!=itr);
				CPPUNIT_ASSERT_EQUAL(Z[size_t(itr)],Ztol[size_t(i)]);
			}
		}
	}

	// center matrix should hold wannier centers
	{
		const auto Wp = readWp(wpath+"wannier90.wout").Wp();
		const auto X = readX(wpath+"wannier90_r.dat");
		const auto Y = readY(wpath+"wannier90_r.dat");
		const auto Z = readZ(wpath+"wannier90_r.dat");
		
		fMat Wp_op(3,Wp.N());
		Wp_op.rAt(0) = lm__::diag(X.center());
		Wp_op.rAt(1) = lm__::diag(Y.center());
		Wp_op.rAt(2) = lm__::diag(Z.center());

		CPPUNIT_ASSERT(Wp==Wp_op);
	}
}

void test_io::test_readEig() {
	const std::string wpath = "data/w90/";
	
	// bad file exception
	CPPUNIT_ASSERT_THROW(readEig("blabla"),std::invalid_argument);

	// InAs_mono
	{
		CPPUNIT_ASSERT_NO_THROW(readEig(wpath+"inas_mono/wannier90.eig"));
		const auto E = readEig(wpath+"inas_mono/wannier90.eig");
			
		const double checksum = sum(E);
		const double refsum = 20.784455054499993;
			
		CPPUNIT_ASSERT_DELTA(refsum,checksum,mtol());
	}
		
	// mos2_dl
	{
		CPPUNIT_ASSERT_NO_THROW(readEig(wpath+"mos2/wannier90.eig"));
		const auto E = readEig(wpath+"mos2/wannier90.eig");
			
		const double checksum = sum(E);
		const double refsum = -3.776307192368124e+03;
			
		CPPUNIT_ASSERT_DELTA(refsum,checksum,mtol());
	}
}

void test_io::test_readWannierTransf() {
	const std::string wpath = "data/w90/";

	// bad file exception
	CPPUNIT_ASSERT_THROW(readWannierTransf("blabla"),std::invalid_argument);
	
	// inas_mono
	{
		CPPUNIT_ASSERT_NO_THROW(readWannierTransf(wpath+"inas_mono/wannier90_u.mat"));
		
		const auto WT = readWannierTransf(wpath+"inas_mono/wannier90_u.mat",
						  wpath+"inas_mono/wannier90_u_dis.mat");
		CPPUNIT_ASSERT(!WT.empty());
		CPPUNIT_ASSERT_EQUAL(size_t(64),WT.N());
		CPPUNIT_ASSERT_EQUAL(size_t(DIM__),WT.dim());
		CPPUNIT_ASSERT_EQUAL(size_t(16),WT.Nb());
		CPPUNIT_ASSERT_EQUAL(size_t(10),WT.Nw());
		CPPUNIT_ASSERT(std::all_of(WT.begin(),WT.end(),[&WT](const cMat& i)->bool
			{ return size(WT.front())==size(i); }));
		CPPUNIT_ASSERT(cunique(WT.k()));
	}
	
	// mos2_dl
	{
		CPPUNIT_ASSERT_NO_THROW(readWannierTransf(wpath+"mos2/wannier90_u.mat"));
		
		const auto WT = readWannierTransf(wpath+"mos2/wannier90_u.mat");
		CPPUNIT_ASSERT(!WT.empty());
		CPPUNIT_ASSERT_EQUAL(size_t(64),WT.N());
		CPPUNIT_ASSERT_EQUAL(size_t(DIM__),WT.dim());
		CPPUNIT_ASSERT_EQUAL(size_t(22),WT.Nb());
		CPPUNIT_ASSERT_EQUAL(size_t(22),WT.Nw());
		CPPUNIT_ASSERT(std::all_of(WT.begin(),WT.end(),[&WT](const cMat& i)->bool
			{ return size(WT.front())==size(i); }));
		CPPUNIT_ASSERT(cunique(WT.k()));
	}
}

void test_io::test_readChk() {
	const std::string wpath = "data/w90/";
	
	// bad file exception
	CPPUNIT_ASSERT_THROW(readWannierTransf("blabla"),std::invalid_argument);

	// inas_mono
	{
		CPPUNIT_ASSERT_NO_THROW(readChk(wpath+"chk/inas_mono/wannier90.chk.fmt"));
		
		const auto WT = readChk(wpath+"chk/inas_mono/wannier90.chk.fmt");
		const auto WTref = readWannierTransf(wpath+"inas_mono/wannier90_u.mat",
						  wpath+"inas_mono/wannier90_u_dis.mat");
		
		CPPUNIT_ASSERT(WTref.k()==WT.k());
		CPPUNIT_ASSERT(WTref.U()==WT.U());
	}
	
	// mos2_dl
	{
		CPPUNIT_ASSERT_NO_THROW(readChk(wpath+"chk/mos2_dl/wannier90.chk.fmt"));
		
		const auto WT = readChk(wpath+"chk/mos2_dl/wannier90.chk.fmt");
		const auto WTref = readWannierTransf(wpath+"mos2/wannier90_u.mat");
		
		CPPUNIT_ASSERT(WTref.k()==WT.k());
		CPPUNIT_ASSERT(WTref.U()==WT.U());
	}
}


const char* test_io::test_id() noexcept {
	return "test_io";
}

CppUnit::Test* test_io::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_writePOSCAR", &test_io::test_writePOSCAR));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_readEf", &test_io::test_readEf));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_genr", &test_io::test_genr));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_readLayerMatrix", &test_io::test_readLayerMatrix));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_printOmf", &test_io::test_printOmf));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_printOlf", &test_io::test_printOlf));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_readOmf", &test_io::test_readOmf));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_readOlf", &test_io::test_readOlf));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_hrDim", &test_io::test_hrDim));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_readB", &test_io::test_readB));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_readAp", &test_io::test_readAp));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_readWp", &test_io::test_readWp));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_readHr", &test_io::test_readHr));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_readAMN", &test_io::test_readAMN));
//	suite->addTest(new CppUnit::TestCaller<test_io>(
//		"test_readXYZ", &test_io::test_readXYZ));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_readEig", &test_io::test_readEig));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_readWannierTransf", &test_io::test_readWannierTransf));
	suite->addTest(new CppUnit::TestCaller<test_io>(
		"test_readChk", &test_io::test_readChk));
	

	return suite;
}
