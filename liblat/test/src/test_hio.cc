// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "test_hio.h"
#include "ll_io.h"
#include "ll_hio.h"
#include "ll_omen.h"
#include "ll_fn.h"
#include "ll_cell.h"
#include "ll_sparse.h"
#include "ll_testTools.h"
#include "lm_testTools.h"
#include <random>
#include <chrono>


using namespace lm__;
using namespace lm__::test;
using namespace ll__;
using namespace ll__::test;


void test_hio::test_writer() {

	// base functions and OMEN version
	{
		const fMat R({.0,.0,-1.0,
			      .0,.0,  .0,
			      .0,.0, 1.0},3);
		const size_t Nw = genRndST(1,100);
		const std::string prefix = "outp/";
		
		ll__::omen::writer writer(R,Nw,prefix);

		// base general information
		CPPUNIT_ASSERT_EQUAL(R.M(),writer.dim());
		CPPUNIT_ASSERT_EQUAL(R.N(),writer.NR());
		CPPUNIT_ASSERT_EQUAL(R,writer.R());
		CPPUNIT_ASSERT(std::vector<uint64_t>(3,0)==writer.nnz());
		CPPUNIT_ASSERT(writer.nnzTest());
		CPPUNIT_ASSERT(!writer.eof());

		// base current block information
		CPPUNIT_ASSERT_EQUAL(size_t(0),writer.c_ind());
		CPPUNIT_ASSERT_EQUAL(R.cFront(),writer.c_R());
		CPPUNIT_ASSERT_EQUAL(uint64_t(0),writer.c_nnz());
		CPPUNIT_ASSERT_EQUAL(size_t(0),writer.c_bytes());

		// OMEN general information
		CPPUNIT_ASSERT_EQUAL(Nw,writer.Nw());
		CPPUNIT_ASSERT_EQUAL(std::string("outp/H_3.bin"),writer.fileName(R.cAt(0)));
		CPPUNIT_ASSERT_EQUAL(std::string("outp/H_4.bin"),writer.fileName(R.cAt(1)));
		CPPUNIT_ASSERT_EQUAL(std::string("outp/H_5.bin"),writer.fileName(R.cAt(2)));
		CPPUNIT_ASSERT_EQUAL(prefix,writer.prefix());

		const std::vector<size_t> N = {genRndST(0,10),genRndST(0,10),genRndST(0,10)};
		const size_t BS = 2*sizeof(size_t)+sizeof(hel);
		
		// try going through blocks inserting random stuff
		for (size_t ind=0; ind!=R.N(); ++ind) {

			writer.newBlock();

			CPPUNIT_ASSERT_EQUAL(ind,writer.c_ind());
			CPPUNIT_ASSERT_EQUAL(R.cAt(ind),writer.c_R());
			CPPUNIT_ASSERT_EQUAL(writer.fileName(R.cAt(ind)),writer.c_fileName());
			CPPUNIT_ASSERT_EQUAL(writer.c_fileName(),writer.c_descr());

			std::vector<sphel> hdat; hdat.reserve(N[ind]);
			while (hdat.size()<hdat.capacity())
				hdat.push_back({genRndST(0,Nw-1),genRndST(0,Nw-1),
					{genRndDouble(-3.0,3.0),genRndDouble(-3.0,3.0)}});
			std::sort(hdat.begin(),hdat.end());

			writer.insert(hdat);
			CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
			CPPUNIT_ASSERT_EQUAL(size_t(BS*N[ind]),writer.c_bytes());
			CPPUNIT_ASSERT(!writer.eof());
			
			writer.flush();
		
			// read file into sparse and compare
			const ll_sparse H(writer.fileName(R.cAt(ind)),"OMEN");
			CPPUNIT_ASSERT_EQUAL(N[ind],H.nnz());

			auto i=hdat.cbegin();
			for (auto j=H.cbegin(),e=H.cend(); j!=e; ++j,++i) {
				CPPUNIT_ASSERT_EQUAL(size_t(i->m),j->m);
				CPPUNIT_ASSERT_EQUAL(size_t(i->n),j->n);
				CPPUNIT_ASSERT_EQUAL(i->h,j->h);
			}
		}

		CPPUNIT_ASSERT(writer.eof());
	}

	// to memory, fMat
	{
		const size_t dim = genRndST(1,10);
		const size_t NR  = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		fMat R = randi<fMat>(dim,NR,-3,3); std::sort(R.cBegin(),R.cEnd());
		std::vector<fMat> H(NR,zeros<fMat>(Nw,Nw));

		ll_writerMEM<fMat> writer(R,H);
		
		std::vector<size_t> N(NR);
		for (auto& i: N) i = genRndST(0,10);
		const size_t BS = 2*sizeof(size_t)+sizeof(double);

		// try going through blocks inserting random stuff
		for (size_t ind=0; ind!=R.N(); ++ind) {

			writer.newBlock();

			CPPUNIT_ASSERT_EQUAL(T(R.cAt(ind)).print(0),writer.c_descr());
			
			std::vector<sphel> hdat; hdat.reserve(N[ind]);
			while (hdat.size()<hdat.capacity()) {
				const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
						{genRndDouble(1.0,3.0),genRndDouble(1.0,3.0)}};
				const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
				if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
			}
			
			writer.insert(hdat);
			CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
			CPPUNIT_ASSERT_EQUAL(size_t(BS*N[ind]),writer.c_bytes());
			CPPUNIT_ASSERT(!writer.eof());

			writer.flush();
			
			// check nnz in matrix manually
			CPPUNIT_ASSERT_EQUAL(N[ind],uint64_t(sum(H[ind].neq(.0))));

			// check data is inserted correctly
			for (const auto& i: hdat)
				CPPUNIT_ASSERT_EQUAL(std::real(i.h),H[ind](i.m,i.n));
		}
		
		CPPUNIT_ASSERT(writer.eof());
	}
	// to memory, cMat
	{
		const size_t dim = genRndST(1,10);
		const size_t NR  = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		fMat R = randi<fMat>(dim,NR,-3,3); std::sort(R.cBegin(),R.cEnd());
		std::vector<cMat> H(NR,zeros<cMat>(Nw,Nw));

		ll_writerMEM<cMat> writer(R,H);
		
		std::vector<size_t> N(NR);
		for (auto& i: N) i = genRndST(0,10);
		const size_t BS = 2*sizeof(size_t)+sizeof(hel);

		// try going through blocks inserting random stuff
		for (size_t ind=0; ind!=R.N(); ++ind) {

			writer.newBlock();

			CPPUNIT_ASSERT_EQUAL(T(R.cAt(ind)).print(0),writer.c_descr());
			
			std::vector<sphel> hdat; hdat.reserve(N[ind]);
			while (hdat.size()<hdat.capacity()) {
				const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
						{genRndDouble(1.0,3.0),genRndDouble(1.0,3.0)}};
				const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
				if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
			}
			
			writer.insert(hdat);
			CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
			CPPUNIT_ASSERT_EQUAL(size_t(BS*N[ind]),writer.c_bytes());
			CPPUNIT_ASSERT(!writer.eof());

			writer.flush();
			
			// check nnz in matrix manually
			CPPUNIT_ASSERT_EQUAL(N[ind],uint64_t(sum(H[ind].neq(.0))));

			// check data is inserted correctly
			for (const auto& i: hdat)
				CPPUNIT_ASSERT_EQUAL(i.h,H[ind](i.m,i.n));
		}
		
		CPPUNIT_ASSERT(writer.eof());
	}
	// to memory, ll_sparse
	{
		const size_t dim = genRndST(1,10);
		const size_t NR  = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		fMat R = randi<fMat>(dim,NR,-3,3); std::sort(R.cBegin(),R.cEnd());
		std::vector<ll_sparse> H(NR,ll_sparse(Nw,Nw));

		ll_writerMEM<ll_sparse> writer(R,H);
		
		std::vector<size_t> N(NR);
		for (auto& i: N) i = genRndST(0,10);
		const size_t BS = sizeof(sphel);

		// try going through blocks inserting random stuff
		for (size_t ind=0; ind!=R.N(); ++ind) {

			writer.newBlock();

			CPPUNIT_ASSERT_EQUAL(T(R.cAt(ind)).print(0),writer.c_descr());
			
			std::vector<sphel> hdat; hdat.reserve(N[ind]);
			while (hdat.size()<hdat.capacity()) {
				const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
						{genRndDouble(1.0,3.0),genRndDouble(1.0,3.0)}};
				const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
				if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
			}
			
			writer.insert(hdat);
			CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
			CPPUNIT_ASSERT_EQUAL(size_t(BS*N[ind]),writer.c_bytes());
			CPPUNIT_ASSERT(!writer.eof());

			writer.flush();

			// check nnz in matrix manually
			CPPUNIT_ASSERT_EQUAL(N[ind],H[ind].nnz());

			// check data is inserted correctly
			for (const auto& i: hdat)
				CPPUNIT_ASSERT_EQUAL(i.h,H[ind](i.m,i.n));
		}
		
		CPPUNIT_ASSERT(writer.eof());
	}

	// R0, fMat
	{
		const size_t dim = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		fMat H = zeros<fMat>(Nw,Nw);

		ll_writerR0<fMat> writer(H,dim);
		
		const size_t N = genRndST(0,10);
		const size_t BS = 2*sizeof(size_t)+sizeof(double);

		// try inserting random stuff
		writer.newBlock();
		CPPUNIT_ASSERT_EQUAL(zeros<fMat>(1,dim).print(0),writer.c_descr());
			
		std::vector<sphel> hdat; hdat.reserve(N);
		while (hdat.size()<hdat.capacity()) {
			const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
					{genRndDouble(1.0,3.0),genRndDouble(1.0,3.0)}};
			const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
			if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
		}
			
		writer.insert(hdat);
		CPPUNIT_ASSERT_EQUAL(uint64_t(N),writer.c_nnz());
		CPPUNIT_ASSERT_EQUAL(size_t(BS*N),writer.c_bytes());
		CPPUNIT_ASSERT(!writer.eof());

		writer.flush();
			
		// check nnz in matrix manually
		CPPUNIT_ASSERT_EQUAL(N,uint64_t(sum(H.neq(.0))));

		// check data is inserted correctly
		for (const auto& i: hdat)
			CPPUNIT_ASSERT_EQUAL(std::real(i.h),H(i.m,i.n));
		
		CPPUNIT_ASSERT(writer.eof());
	}
	// R0, cMat
	{
		const size_t dim = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		cMat H = zeros<cMat>(Nw,Nw);

		ll_writerR0<cMat> writer(H,dim);
		
		const size_t N = genRndST(0,10);
		const size_t BS = 2*sizeof(size_t)+sizeof(hel);

		// try inserting random stuff
		writer.newBlock();
		CPPUNIT_ASSERT_EQUAL(zeros<fMat>(1,dim).print(0),writer.c_descr());
			
		std::vector<sphel> hdat; hdat.reserve(N);
		while (hdat.size()<hdat.capacity()) {
			const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
					{genRndDouble(1.0,3.0),genRndDouble(1.0,3.0)}};
			const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
			if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
		}
			
		writer.insert(hdat);
		CPPUNIT_ASSERT_EQUAL(uint64_t(N),writer.c_nnz());
		CPPUNIT_ASSERT_EQUAL(size_t(BS*N),writer.c_bytes());
		CPPUNIT_ASSERT(!writer.eof());

		writer.flush();
			
		// check nnz in matrix manually
		CPPUNIT_ASSERT_EQUAL(N,uint64_t(sum(H.neq(.0))));

		// check data is inserted correctly
		for (const auto& i: hdat)
			CPPUNIT_ASSERT_EQUAL(i.h,H(i.m,i.n));
		
		CPPUNIT_ASSERT(writer.eof());
	}
	// R0, ll_sparse
	{
		const size_t dim = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		ll_sparse H(Nw,Nw);

		ll_writerR0<ll_sparse> writer(H,dim);
		
		const size_t N = genRndST(0,10);
		const size_t BS = sizeof(sphel);

		// try inserting random stuff
		writer.newBlock();
		CPPUNIT_ASSERT_EQUAL(zeros<fMat>(1,dim).print(0),writer.c_descr());
			
		std::vector<sphel> hdat; hdat.reserve(N);
		while (hdat.size()<hdat.capacity()) {
			const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
					{genRndDouble(1.0,3.0),genRndDouble(1.0,3.0)}};
			const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
			if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
		}
			
		writer.insert(hdat);
		CPPUNIT_ASSERT_EQUAL(uint64_t(N),writer.c_nnz());
		CPPUNIT_ASSERT_EQUAL(size_t(BS*N),writer.c_bytes());
		CPPUNIT_ASSERT(!writer.eof());

		writer.flush();
			
		// check nnz in matrix manually
		CPPUNIT_ASSERT_EQUAL(N,H.nnz());

		// check data is inserted correctly
		for (const auto& i: hdat)
			CPPUNIT_ASSERT_EQUAL(i.h,H(i.m,i.n));
		
		CPPUNIT_ASSERT(writer.eof());
	}

	// w90
	{
		const std::string fileName = "outp/writerW90_hr.dat";
		
		const size_t dim =  genRndST(1,10);
		const size_t NR  =  genRndST(1,10);
		const size_t Nw  =  genRndST(50,100);
		const size_t prec = genRndST(3,12);

		fMat R = randi<fMat>(dim,NR,-3,3); std::sort(R.cBegin(),R.cEnd());
		std::vector<cMat> H(NR,zeros<cMat>(Nw,Nw));

		ll_writerW90 writer(fileName,R,Nw,prec);
		
		CPPUNIT_ASSERT_EQUAL(prec,writer.precision());
		CPPUNIT_ASSERT_EQUAL(fileName,writer.fileName());
		
		std::vector<size_t> N(NR);
		for (auto& i: N) i = genRndST(0,10);

		// try going through blocks inserting random stuff
		for (size_t ind=0; ind!=R.N(); ++ind) {

			writer.newBlock();

			CPPUNIT_ASSERT_EQUAL(T(R.cAt(ind)).print(0),writer.c_descr());
			
			std::vector<sphel> hdat; hdat.reserve(N[ind]);
			while (hdat.size()<hdat.capacity()) {
				const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
						{genRndDouble(1.0,3.0),genRndDouble(1.0,3.0)}};
				const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
				if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
			}
			
			writer.insert(hdat);
			CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
			CPPUNIT_ASSERT(!writer.eof());

			writer.flush();
		}
		
		CPPUNIT_ASSERT(writer.eof());

	}


	// binary all in one file, uint16, float
	{
		typedef uint16_t IT;
		typedef float FT;

		const std::string fileName = "outp/hrssr.bin";
		
		const size_t dim = genRndST(1,10);
		const size_t NR  = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		fMat R = randi<fMat>(dim,NR,-3,3); std::sort(R.cBegin(),R.cEnd());

		std::vector<size_t> N(NR);
		for (auto& i: N) i = genRndST(0,10);
		const size_t BS = 2*sizeof(IT)+sizeof(FT);

		std::vector<std::vector<sphel>> hdat_; hdat_.reserve(NR);
		
		{
			ll_writerBIN<IT,FT> writer(fileName,R,Nw);
			CPPUNIT_ASSERT_EQUAL(fileName,writer.fileName());

			// try going through blocks inserting random stuff
			for (size_t ind=0; ind!=R.N(); ++ind) {

				writer.newBlock();

				CPPUNIT_ASSERT_EQUAL(T(R.cAt(ind)).print(0),writer.c_descr());

				hdat_.push_back(std::vector<sphel>());
				std::vector<sphel>& hdat = hdat_.back();
				
				hdat.reserve(N[ind]);
				while (hdat.size()<hdat.capacity()) {
					const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
							{genRndDouble(-3.0,3.0),genRndDouble(-3.0,3.0)}};
					const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
					if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
				}
				
				writer.insert(hdat);
				CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
				CPPUNIT_ASSERT_EQUAL(size_t(BS*N[ind]),writer.c_bytes());
				CPPUNIT_ASSERT(!writer.eof());

				writer.flush();
			}
		
			CPPUNIT_ASSERT(writer.eof());
		}

		// read file into ll_sparse and compare
		const auto hr = readHrSparse(fileName);
		CPPUNIT_ASSERT_EQUAL(R,hr.R());
		for (size_t ind=0; ind!=R.N(); ++ind) {
			CPPUNIT_ASSERT_EQUAL(N[ind],hr[ind].nnz());
			
			auto i=hr[ind].cbegin();
			for (auto j=hdat_[ind].cbegin(),e=hdat_[ind].cend(); j!=e; ++j,++i) {
				CPPUNIT_ASSERT_EQUAL(size_t(j->m),size_t(i->m));
				CPPUNIT_ASSERT_EQUAL(size_t(j->n),size_t(i->n));
				CPPUNIT_ASSERT_EQUAL(FT(std::real(j->h)),FT(std::real(i->h)));
			}
		}
	}
	// binary all in one file, uint32, float
	{
		typedef uint32_t IT;
		typedef float FT;

		const std::string fileName = "outp/hrlsr.bin";
		
		const size_t dim = genRndST(1,10);
		const size_t NR  = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		fMat R = randi<fMat>(dim,NR,-3,3); std::sort(R.cBegin(),R.cEnd());

		std::vector<size_t> N(NR);
		for (auto& i: N) i = genRndST(0,10);
		const size_t BS = 2*sizeof(IT)+sizeof(FT);

		std::vector<std::vector<sphel>> hdat_; hdat_.reserve(NR);
		
		{
			ll_writerBIN<IT,FT> writer(fileName,R,Nw);
			CPPUNIT_ASSERT_EQUAL(fileName,writer.fileName());

			// try going through blocks inserting random stuff
			for (size_t ind=0; ind!=R.N(); ++ind) {

				writer.newBlock();

				CPPUNIT_ASSERT_EQUAL(T(R.cAt(ind)).print(0),writer.c_descr());

				hdat_.push_back(std::vector<sphel>());
				std::vector<sphel>& hdat = hdat_.back();
				
				hdat.reserve(N[ind]);
				while (hdat.size()<hdat.capacity()) {
					const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
							{genRndDouble(-3.0,3.0),genRndDouble(-3.0,3.0)}};
					const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
					if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
				}
				
				writer.insert(hdat);
				CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
				CPPUNIT_ASSERT_EQUAL(size_t(BS*N[ind]),writer.c_bytes());
				CPPUNIT_ASSERT(!writer.eof());

				writer.flush();
			}
		
			CPPUNIT_ASSERT(writer.eof());
		}

		// read file into ll_sparse and compare
		const auto hr = readHrSparse(fileName);
		CPPUNIT_ASSERT_EQUAL(R,hr.R());
		for (size_t ind=0; ind!=R.N(); ++ind) {
			CPPUNIT_ASSERT_EQUAL(N[ind],hr[ind].nnz());
			
			auto i=hr[ind].cbegin();
			for (auto j=hdat_[ind].cbegin(),e=hdat_[ind].cend(); j!=e; ++j,++i) {
				CPPUNIT_ASSERT_EQUAL(size_t(j->m),size_t(i->m));
				CPPUNIT_ASSERT_EQUAL(size_t(j->n),size_t(i->n));
				CPPUNIT_ASSERT_EQUAL(FT(std::real(j->h)),FT(std::real(i->h)));
			}
		}
	}
	// binary all in one file, uint16, std::complex<float>
	{
		typedef uint16_t IT;
		typedef std::complex<float> FT;

		const std::string fileName = "outp/hrssc.bin";
		
		const size_t dim = genRndST(1,10);
		const size_t NR  = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		fMat R = randi<fMat>(dim,NR,-3,3); std::sort(R.cBegin(),R.cEnd());

		std::vector<size_t> N(NR);
		for (auto& i: N) i = genRndST(0,10);
		const size_t BS = 2*sizeof(IT)+sizeof(FT);

		std::vector<std::vector<sphel>> hdat_; hdat_.reserve(NR);
		
		{
			ll_writerBIN<IT,FT> writer(fileName,R,Nw);
			CPPUNIT_ASSERT_EQUAL(fileName,writer.fileName());

			// try going through blocks inserting random stuff
			for (size_t ind=0; ind!=R.N(); ++ind) {

				writer.newBlock();

				CPPUNIT_ASSERT_EQUAL(T(R.cAt(ind)).print(0),writer.c_descr());

				hdat_.push_back(std::vector<sphel>());
				std::vector<sphel>& hdat = hdat_.back();
				
				hdat.reserve(N[ind]);
				while (hdat.size()<hdat.capacity()) {
					const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
							{genRndDouble(-3.0,3.0),genRndDouble(-3.0,3.0)}};
					const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
					if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
				}
				
				writer.insert(hdat);
				CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
				CPPUNIT_ASSERT_EQUAL(size_t(BS*N[ind]),writer.c_bytes());
				CPPUNIT_ASSERT(!writer.eof());

				writer.flush();
			}
		
			CPPUNIT_ASSERT(writer.eof());
		}

		// read file into ll_sparse and compare
		const auto hr = readHrSparse(fileName);
		CPPUNIT_ASSERT_EQUAL(R,hr.R());
		for (size_t ind=0; ind!=R.N(); ++ind) {
			CPPUNIT_ASSERT_EQUAL(N[ind],hr[ind].nnz());
			
			auto i=hr[ind].cbegin();
			for (auto j=hdat_[ind].cbegin(),e=hdat_[ind].cend(); j!=e; ++j,++i) {
				CPPUNIT_ASSERT_EQUAL(size_t(j->m),size_t(i->m));
				CPPUNIT_ASSERT_EQUAL(size_t(j->n),size_t(i->n));
				CPPUNIT_ASSERT_EQUAL(FT(j->h),FT(i->h));
			}
		}
	}
	// binary all in one file, uint32, std::complex<float>
	{
		typedef uint32_t IT;
		typedef std::complex<float> FT;

		const std::string fileName = "outp/hrssc.bin";
		
		const size_t dim = genRndST(1,10);
		const size_t NR  = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		fMat R = randi<fMat>(dim,NR,-3,3); std::sort(R.cBegin(),R.cEnd());

		std::vector<size_t> N(NR);
		for (auto& i: N) i = genRndST(0,10);
		const size_t BS = 2*sizeof(IT)+sizeof(FT);

		std::vector<std::vector<sphel>> hdat_; hdat_.reserve(NR);
		
		{
			ll_writerBIN<IT,FT> writer(fileName,R,Nw);
			CPPUNIT_ASSERT_EQUAL(fileName,writer.fileName());

			// try going through blocks inserting random stuff
			for (size_t ind=0; ind!=R.N(); ++ind) {

				writer.newBlock();

				CPPUNIT_ASSERT_EQUAL(T(R.cAt(ind)).print(0),writer.c_descr());

				hdat_.push_back(std::vector<sphel>());
				std::vector<sphel>& hdat = hdat_.back();
				
				hdat.reserve(N[ind]);
				while (hdat.size()<hdat.capacity()) {
					const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
							{genRndDouble(-3.0,3.0),genRndDouble(-3.0,3.0)}};
					const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
					if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
				}
				
				writer.insert(hdat);
				CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
				CPPUNIT_ASSERT_EQUAL(size_t(BS*N[ind]),writer.c_bytes());
				CPPUNIT_ASSERT(!writer.eof());

				writer.flush();
			}
		
			CPPUNIT_ASSERT(writer.eof());
		}

		// read file into ll_sparse and compare
		const auto hr = readHrSparse(fileName);
		CPPUNIT_ASSERT_EQUAL(R,hr.R());
		for (size_t ind=0; ind!=R.N(); ++ind) {
			CPPUNIT_ASSERT_EQUAL(N[ind],hr[ind].nnz());
			
			auto i=hr[ind].cbegin();
			for (auto j=hdat_[ind].cbegin(),e=hdat_[ind].cend(); j!=e; ++j,++i) {
				CPPUNIT_ASSERT_EQUAL(size_t(j->m),size_t(i->m));
				CPPUNIT_ASSERT_EQUAL(size_t(j->n),size_t(i->n));
				CPPUNIT_ASSERT_EQUAL(FT(j->h),FT(i->h));
			}
		}
	}
	
	
	// binary all in one file, uint16, double
	{
		typedef uint16_t IT;
		typedef double FT;

		const std::string fileName = "outp/hrsdr.bin";
		
		const size_t dim = genRndST(1,10);
		const size_t NR  = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		fMat R = randi<fMat>(dim,NR,-3,3); std::sort(R.cBegin(),R.cEnd());

		std::vector<size_t> N(NR);
		for (auto& i: N) i = genRndST(0,10);
		const size_t BS = 2*sizeof(IT)+sizeof(FT);

		std::vector<std::vector<sphel>> hdat_; hdat_.reserve(NR);
		
		{
			ll_writerBIN<IT,FT> writer(fileName,R,Nw);
			CPPUNIT_ASSERT_EQUAL(fileName,writer.fileName());

			// try going through blocks inserting random stuff
			for (size_t ind=0; ind!=R.N(); ++ind) {

				writer.newBlock();

				CPPUNIT_ASSERT_EQUAL(T(R.cAt(ind)).print(0),writer.c_descr());

				hdat_.push_back(std::vector<sphel>());
				std::vector<sphel>& hdat = hdat_.back();
				
				hdat.reserve(N[ind]);
				while (hdat.size()<hdat.capacity()) {
					const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
							{genRndDouble(-3.0,3.0),genRndDouble(-3.0,3.0)}};
					const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
					if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
				}
				
				writer.insert(hdat);
				CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
				CPPUNIT_ASSERT_EQUAL(size_t(BS*N[ind]),writer.c_bytes());
				CPPUNIT_ASSERT(!writer.eof());

				writer.flush();
			}
		
			CPPUNIT_ASSERT(writer.eof());
		}

		// read file into ll_sparse and compare
		const auto hr = readHrSparse(fileName);
		CPPUNIT_ASSERT_EQUAL(R,hr.R());
		for (size_t ind=0; ind!=R.N(); ++ind) {
			CPPUNIT_ASSERT_EQUAL(N[ind],hr[ind].nnz());
			
			auto i=hr[ind].cbegin();
			for (auto j=hdat_[ind].cbegin(),e=hdat_[ind].cend(); j!=e; ++j,++i) {
				CPPUNIT_ASSERT_EQUAL(size_t(j->m),size_t(i->m));
				CPPUNIT_ASSERT_EQUAL(size_t(j->n),size_t(i->n));
				CPPUNIT_ASSERT_EQUAL(FT(std::real(j->h)),FT(std::real(i->h)));
			}
		}
	}
	// binary all in one file, uint32, double
	{
		typedef uint32_t IT;
		typedef double FT;

		const std::string fileName = "outp/hrldr.bin";
		
		const size_t dim = genRndST(1,10);
		const size_t NR  = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		fMat R = randi<fMat>(dim,NR,-3,3); std::sort(R.cBegin(),R.cEnd());

		std::vector<size_t> N(NR);
		for (auto& i: N) i = genRndST(0,10);
		const size_t BS = 2*sizeof(IT)+sizeof(FT);

		std::vector<std::vector<sphel>> hdat_; hdat_.reserve(NR);
		
		{
			ll_writerBIN<IT,FT> writer(fileName,R,Nw);
			CPPUNIT_ASSERT_EQUAL(fileName,writer.fileName());

			// try going through blocks inserting random stuff
			for (size_t ind=0; ind!=R.N(); ++ind) {

				writer.newBlock();

				CPPUNIT_ASSERT_EQUAL(T(R.cAt(ind)).print(0),writer.c_descr());

				hdat_.push_back(std::vector<sphel>());
				std::vector<sphel>& hdat = hdat_.back();
				
				hdat.reserve(N[ind]);
				while (hdat.size()<hdat.capacity()) {
					const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
							{genRndDouble(-3.0,3.0),genRndDouble(-3.0,3.0)}};
					const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
					if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
				}
				
				writer.insert(hdat);
				CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
				CPPUNIT_ASSERT_EQUAL(size_t(BS*N[ind]),writer.c_bytes());
				CPPUNIT_ASSERT(!writer.eof());

				writer.flush();
			}
		
			CPPUNIT_ASSERT(writer.eof());
		}

		// read file into ll_sparse and compare
		const auto hr = readHrSparse(fileName);
		CPPUNIT_ASSERT_EQUAL(R,hr.R());
		for (size_t ind=0; ind!=R.N(); ++ind) {
			CPPUNIT_ASSERT_EQUAL(N[ind],hr[ind].nnz());
			
			auto i=hr[ind].cbegin();
			for (auto j=hdat_[ind].cbegin(),e=hdat_[ind].cend(); j!=e; ++j,++i) {
				CPPUNIT_ASSERT_EQUAL(size_t(j->m),size_t(i->m));
				CPPUNIT_ASSERT_EQUAL(size_t(j->n),size_t(i->n));
				CPPUNIT_ASSERT_EQUAL(FT(std::real(j->h)),FT(std::real(i->h)));
			}
		}
	}
	// binary all in one file, uint16, std::complex<double>
	{
		typedef uint16_t IT;
		typedef std::complex<double> FT;

		const std::string fileName = "outp/hrsdc.bin";
		
		const size_t dim = genRndST(1,10);
		const size_t NR  = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		fMat R = randi<fMat>(dim,NR,-3,3); std::sort(R.cBegin(),R.cEnd());

		std::vector<size_t> N(NR);
		for (auto& i: N) i = genRndST(0,10);
		const size_t BS = 2*sizeof(IT)+sizeof(FT);

		std::vector<std::vector<sphel>> hdat_; hdat_.reserve(NR);
		
		{
			ll_writerBIN<IT,FT> writer(fileName,R,Nw);
			CPPUNIT_ASSERT_EQUAL(fileName,writer.fileName());

			// try going through blocks inserting random stuff
			for (size_t ind=0; ind!=R.N(); ++ind) {

				writer.newBlock();

				CPPUNIT_ASSERT_EQUAL(T(R.cAt(ind)).print(0),writer.c_descr());

				hdat_.push_back(std::vector<sphel>());
				std::vector<sphel>& hdat = hdat_.back();
				
				hdat.reserve(N[ind]);
				while (hdat.size()<hdat.capacity()) {
					const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
							{genRndDouble(-3.0,3.0),genRndDouble(-3.0,3.0)}};
					const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
					if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
				}
				
				writer.insert(hdat);
				CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
				CPPUNIT_ASSERT_EQUAL(size_t(BS*N[ind]),writer.c_bytes());
				CPPUNIT_ASSERT(!writer.eof());

				writer.flush();
			}
		
			CPPUNIT_ASSERT(writer.eof());
		}

		// read file into ll_sparse and compare
		const auto hr = readHrSparse(fileName);
		CPPUNIT_ASSERT_EQUAL(R,hr.R());
		for (size_t ind=0; ind!=R.N(); ++ind) {
			CPPUNIT_ASSERT_EQUAL(N[ind],hr[ind].nnz());
			
			auto i=hr[ind].cbegin();
			for (auto j=hdat_[ind].cbegin(),e=hdat_[ind].cend(); j!=e; ++j,++i) {
				CPPUNIT_ASSERT_EQUAL(size_t(j->m),size_t(i->m));
				CPPUNIT_ASSERT_EQUAL(size_t(j->n),size_t(i->n));
				CPPUNIT_ASSERT_EQUAL(FT(j->h),FT(i->h));
			}
		}
	}
	// binary all in one file, uint32, std::complex<double>
	{
		typedef uint32_t IT;
		typedef std::complex<double> FT;

		const std::string fileName = "outp/hrsdc.bin";
		
		const size_t dim = genRndST(1,10);
		const size_t NR  = genRndST(1,10);
		const size_t Nw  = genRndST(50,100);

		fMat R = randi<fMat>(dim,NR,-3,3); std::sort(R.cBegin(),R.cEnd());

		std::vector<size_t> N(NR);
		for (auto& i: N) i = genRndST(0,10);
		const size_t BS = 2*sizeof(IT)+sizeof(FT);

		std::vector<std::vector<sphel>> hdat_; hdat_.reserve(NR);
		
		{
			ll_writerBIN<IT,FT> writer(fileName,R,Nw);
			CPPUNIT_ASSERT_EQUAL(fileName,writer.fileName());

			// try going through blocks inserting random stuff
			for (size_t ind=0; ind!=R.N(); ++ind) {

				writer.newBlock();

				CPPUNIT_ASSERT_EQUAL(T(R.cAt(ind)).print(0),writer.c_descr());

				hdat_.push_back(std::vector<sphel>());
				std::vector<sphel>& hdat = hdat_.back();
				
				hdat.reserve(N[ind]);
				while (hdat.size()<hdat.capacity()) {
					const sphel ne = {genRndST(0,Nw-1),genRndST(0,Nw-1),
							{genRndDouble(-3.0,3.0),genRndDouble(-3.0,3.0)}};
					const auto itr = std::lower_bound(hdat.begin(),hdat.end(),ne);
					if (itr==hdat.end() || *itr!=ne) hdat.insert(itr,ne);
				}
				
				writer.insert(hdat);
				CPPUNIT_ASSERT_EQUAL(uint64_t(N[ind]),writer.c_nnz());
				CPPUNIT_ASSERT_EQUAL(size_t(BS*N[ind]),writer.c_bytes());
				CPPUNIT_ASSERT(!writer.eof());

				writer.flush();
			}
		
			CPPUNIT_ASSERT(writer.eof());
		}

		// read file into ll_sparse and compare
		const auto hr = readHrSparse(fileName);
		CPPUNIT_ASSERT_EQUAL(R,hr.R());
		for (size_t ind=0; ind!=R.N(); ++ind) {
			CPPUNIT_ASSERT_EQUAL(N[ind],hr[ind].nnz());
			
			auto i=hr[ind].cbegin();
			for (auto j=hdat_[ind].cbegin(),e=hdat_[ind].cend(); j!=e; ++j,++i) {
				CPPUNIT_ASSERT_EQUAL(size_t(j->m),size_t(i->m));
				CPPUNIT_ASSERT_EQUAL(size_t(j->n),size_t(i->n));
				CPPUNIT_ASSERT_EQUAL(FT(j->h),FT(i->h));
			}
		}
	}
}

void test_hio::test_getConnectedGrid() {

	const std::string path = "data/w90/mos2/";
	const std::string prefix = "outp/";

	// input file for hbonds
	ll_wmatching_input inp;
	inp.wout = path+"wannier90.wout";
	inp.hrdat = path+"wannier90_hr.dat";
	inp.mode = {"atomic"};
	inp.tol = 0.0;
	inp.rm_unmatched = 1;
	inp.verbosity = 0;

	const ll_hbonds W(inp,std::cout); W.setQueryTolDirect(WTOL__);
	const auto hr = readHr(inp.hrdat);

	// find R grid and compare to reference
	const auto R = getConnectedGrid(W,rv(W.dim(),false),DBL_MAX,false);

	CPPUNIT_ASSERT_EQUAL(hr.R(),R);

}

void test_hio::test_hctor() {
	
	const std::string path = "data/w90/mos2/";
	const std::string prefix = "outp/";

	// input file for hbonds
	ll_wmatching_input inp;
	inp.wout = path+"wannier90.wout";
	inp.hrdat = path+"wannier90_hr.dat";
	inp.mode = {"atomic"};
	inp.tol = 0.0;
	inp.rm_unmatched = 1;
	inp.verbosity = 0;

	const ll_hbonds W(inp,std::cout); W.setQueryTolDirect(WTOL__);
	const auto hr_ = readHr(inp.hrdat);

	// generate hr like that from wannier90 and compare
	{
		R_H<cMat> hr(hr_.R(), std::vector<cMat>(hr_.size(),zeros<cMat>(hr_.Nw(),hr_.Nw())));
		hctor(W,ll_writerMEM<cMat>(hr),DBL_MAX,false);

		// check self adjointedness, i.e. H(R) = H(-R)^T
		const double htol = 1e-6;
		for (size_t i=0; i!=hr.size(); ++i) {
			CPPUNIT_ASSERT(max(abs(hr.H()[i] - T(hr.H()[hr.size()-1-i])))<htol);
		}

		// compare eigenvalues at gamma point
		const auto E_ = eigh(std::accumulate(hr_.cbegin(),hr_.cend(),zeros<cMat>(hr_.Nw(),hr_.Nw())));
		const auto E = eigh(std::accumulate(hr.cbegin(),hr.cend(),zeros<cMat>(hr.Nw(),hr.Nw())));
		CPPUNIT_ASSERT(max(abs(E_-E))<1e-6);
	}
	
	// use hctor implicitly via genHam to compute bandstructure on a super cell
	{
		fMat c;
		do {
			c = randi<fMat>(W.dim()-1,W.dim()-1,-2,2);
		} while (std::abs(det(c))<2.0);
		fMat C({c[0], 0.0, c[1],
			 0.0, 1.0,  0.0,
		        c[2], 0.0, c[3]},3);

		const auto scell = W.cell().copy().expand(C);
		const auto shr = genHam<cMat>(scell,W.r(),W,true);
		
		// compare bandstructures along trace
		fMat kpts({-.5,.0,.0,
			    .0,.0,.0,
			    .5,.0,.0,
			    .5,.0,.5,
			    .0,.0,.5,
			    .0,.0,.0},3);
		const size_t Np = genRndST(30,50);
		const auto PP = genPath(kpts,Np);

		const auto E = calcBS(shr,PP.path(),0);
		const auto E_ = calcFoldedBS(hr_,PP.path(),scell.B(),W.cell().B(),0);

		CPPUNIT_ASSERT(max(abs(E_-E))<1e-3);
	}
}


const char* test_hio::test_id() noexcept {
	return "test_hio";
}

CppUnit::Test* test_hio::suite() {
	CppUnit::TestSuite* suite = new CppUnit::TestSuite(test_id());
	
	suite->addTest(new CppUnit::TestCaller<test_hio>(
		"test_writer", &test_hio::test_writer));
	suite->addTest(new CppUnit::TestCaller<test_hio>(
		"test_getConnectedGrid", &test_hio::test_getConnectedGrid));
	suite->addTest(new CppUnit::TestCaller<test_hio>(
		"test_hctor", &test_hio::test_hctor));

	return suite;
}
