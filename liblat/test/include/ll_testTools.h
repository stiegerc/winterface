// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _LL_TESTTOOLS_
#define _LL_TESTTOOLS_

#include "ll_cell.h"
#include "lm_testTools.h"
#include "lm_fn.h"
#include <random>
#include <vector>
#include <cmath>


namespace ll__ {
namespace test {
	// random
	inline double genRndDouble(const double l=2.0, const double u=3.0) noexcept {
		assert(l<=u);
	
		std::random_device rd;
		std::mt19937_64 gen(rd());
		std::uniform_real_distribution<double> dis(l,u);
	
		return dis(gen);
	}
	inline std::vector<size_t> genRndaCv(const size_t N) noexcept {
		std::vector<size_t> buff(N);
		for (auto& i: buff) i=lm__::test::genRndST(0,N-1);
		std::sort(buff.begin(),buff.end());
		
		auto unq = buff;
		const auto e = std::unique(unq.begin(),unq.end());

		aCv res; res.reserve(std::distance(unq.begin(),e));
		for (auto i=unq.begin(); i!=e; ++i)
			res.push_back(std::count(buff.begin(),buff.end(),*i));

		return res;
	}
	inline rv genRndrv(const size_t D=3, const size_t Nt=NPOS__) noexcept {
		if (Nt==D) return rv(D,true);
		if (Nt==0) return rv(D,false);
		
		std::vector<bool> res; res.reserve(D);
		if (Nt==NPOS__)
			for (size_t i=0; i!=D; ++i)
				res.push_back(genRndDouble(0.0,1.0)<.5);
		else
			do {
				res.clear();
				for (size_t i=0; i!=D; ++i)
					res.push_back(genRndDouble(0.0,1.0)<.5);
			} while (Ntrue(res)!=Nt);
		return res;
	}
	inline fMat rnd_i(const size_t D=3, const int l=-3.0, const int u=3.0) noexcept {
		assert(l<u);
		
		const double d = std::pow((u-l)/2.0,D);
		fMat res;
		do res = randi<fMat>(D,D,l,u);
		while (std::abs(lm__::det(res))<d);
		return res;
	}
	inline std::vector<size_t> genRndIvec(const size_t N, const size_t mi, const size_t ma) noexcept {
		assert(mi<ma);

		std::vector<size_t> res(N); res.reserve(N);
		while (res.size()!=res.capacity())
			res.push_back(lm__::test::genRndST(mi,ma));

		std::sort(res.begin(),res.end());
		res.resize(std::distance(res.begin(),std::unique(res.begin(),res.end())));
		res.shrink_to_fit();

		return res;
	}
	inline aTv genRndaTv(const ll_cell& inp) noexcept {
		aTv res(inp.Nspecies());
		std::iota(res.begin(),res.end(),aT(0));

		for (size_t i=0,e=lm__::test::genRndST(0,inp.Nspecies()-1); i!=e; ++i)
			res.erase(res.begin()+lm__::test::genRndST(0,res.size()-1));
		return res;
	}
	inline wi genRndWi(const size_t N, const size_t n) noexcept {
		assert(n<=N);

		ll__::wi res(n);
		for (size_t i=0; i!=N; ++i)
			res[lm__::test::genRndST(0,n-1)].push_back(i);
		for (auto& i: res)
			if (i.empty()) {
				for (auto& j: res)
					if (j.size()>=2) {
						i.push_back(j.back());
						j.pop_back();
						break;
					}
			}
		return res;
	}
	inline fMat genRndR(const size_t D) noexcept {
		if (D==1) return fMat({1.0});
		auto res = gsorth(rand<fMat>(D,D));
		if (det(res)<0.0)
			std::iter_swap(res.cBegin(),res.cBegin()+1);
		return res;
	}
	inline char genRndLetter() noexcept { return char(lm__::test::genRndST(97,122)); }
	inline char genRndCapLetter() noexcept { return char(lm__::test::genRndST(65,90)); }
	inline std::string genRndAtomicString(const size_t wl=2) noexcept {
		if (!wl) return "";
		std::string res; res.reserve(wl);
		res.push_back(genRndCapLetter());
		while (res.size()!=wl)
			res.push_back(genRndLetter());
		return res;
	}
	inline idv genRndIdv(const size_t n, const size_t wl=2) noexcept {
		idv res(n);
		for (auto& i: res) i=genRndAtomicString(wl);
		std::sort(res.begin(),res.end());
		return res;
	}

	// fundamental
	inline ll_cell genCubic(const size_t D=3, const double a=genRndDouble()) noexcept {
		assert(D>0);
		using namespace lm__;

		auto B = a*eye<fMat>(D);
		auto Ap = zeros<fMat>(D,1);
		const aCv N = {1};

		return ll_cell(std::move(B),std::move(Ap),N,genRndIdv(1,2));
	}
	inline ll_cell genOrthorhombic(const double a=genRndDouble(), const double b=genRndDouble(),
					const double c=genRndDouble()) noexcept {
		using namespace lm__;
		
		auto B = diag(fMat({a,b,c}));
		auto Ap = zeros<fMat>(3,1);
		const aCv N = {1};

		return ll_cell(std::move(B),std::move(Ap),N,genRndIdv(1,2));
	}
	inline ll_cell genFCC(const size_t D=3, const double a=genRndDouble()) noexcept {
		assert(D>0);
		using namespace lm__;
		if (D==1) return genCubic(D,a);

		auto B = a/double(D-1)*(ones<fMat>(D)-eye<fMat>(D));
		auto Ap = zeros<fMat>(D,1);
		aCv N = {1};

		return ll_cell(std::move(B),std::move(Ap),N,genRndIdv(1,2));
	}
	inline ll_cell genBCC(const size_t D=3, const double a=genRndDouble()) noexcept {
		assert(D>0);
		using namespace lm__;
		if (D<3) return genCubic(D,a);

		auto B = a/double(D-2)*(ones<fMat>(D)-2.0*eye<fMat>(D));
		auto Ap = zeros<fMat>(D,1);
		aCv N = {1};

		return ll_cell(std::move(B),std::move(Ap),N,genRndIdv(1,2));
	}
	inline ll_cell genHexagonal(const double a=genRndDouble(),
				const double c=genRndDouble()) noexcept {
		fMat B({.5*a,-.5*a*sqrt(3.0),0.0,
			.5*a, .5*a*sqrt(3.0),0.0,
			 0.0,            0.0,  c},3);
		fMat Ap({0.0,        0.0,0.0,
			 2.0/3.0,2.0/3.0,0.0},3);
		aCv N = {1,1};
		
		return ll_cell(std::move(B),std::move(Ap),N,genRndIdv(2,2));
	}

	// compound
	inline ll_cell genZincblende(const size_t D=3, const double a=genRndDouble()) noexcept {
		assert(D>0);
		using namespace lm__;

		auto B = (D==1) ? a*ones<fMat>(D): a/double(D-1)*(ones<fMat>(D)-eye<fMat>(D));
		fMat Ap(D,0); Ap.reserve(2);
		Ap.push_back(zeros<fMat>(D,1));
		Ap.push_back(.25*ones<fMat>(D,1));
		aCv N = {1,1};

		return ll_cell(std::move(B),std::move(Ap),N,genRndIdv(2,2));
	}
	inline ll_cell genDiChalcogenide(const double a=genRndDouble(), const double c=genRndDouble(),
					const double d=genRndDouble(.1,.2)) noexcept {
		fMat B({.5*a,-.5*a*sqrt(3.0),0.0,
			.5*a, .5*a*sqrt(3.0),0.0,
			 0.0,            0.0,  c},3);
		fMat Ap({0.0,        0.0,  0.0,
			 2.0/3.0,2.0/3.0,    d,
			 2.0/3.0,2.0/3.0,1.0-d},3);
		aCv N = {1,2};
		
		return ll_cell(std::move(B),std::move(Ap),N,genRndIdv(2,2));
	}

	// true random
	inline ll_cell genRndCell(const size_t D=3, const size_t Nl=1, const size_t Nu=10) noexcept {
		assert(Nl<=Nu);
		assert(D>0);
		using namespace lm__;

		const auto N = lm__::test::genRndST(Nl,Nu);
		fMat B;
		do B = rand<fMat>(D,D);
		while (std::abs(det(B))<.5);
		
		// extend Ap with minimal spacing lambda
		fMat Ap(D,0); Ap.reserve(N);
		const auto msl = [&Ap,D](const double sp) -> void {
			const auto cont = rand<fMat>(D,1,0.0,1.0-2.0*mtol());
			for (auto ii=Ap.cBegin(), ie=Ap.cEnd(); ii!=ie; ++ii) {
				if (norm(cont-*ii)<sp) return;
			}
			Ap << cont;
		};
		const double sp = 1.0/(N*5.0);
		do msl(sp); while(Ap.N()!=Ap.ccap());
	
		auto NN = genRndaCv(N);
		auto id = genRndIdv(NN.size(),2);
		return ll_cell(std::move(B),std::move(Ap),std::move(NN),std::move(id));
	}

	// random
	inline ll_cell genRandom() noexcept {
		switch (lm__::test::genRndST(0,7)) {
			case 0: return genCubic();
			case 1: return genOrthorhombic();
			case 2: return genFCC();
			case 3: return genBCC();
			case 4: return genHexagonal();
			case 5: return genZincblende();
			case 6: return genDiChalcogenide();
			case 7: return genRndCell();
			default: return genZincblende();
		}
		return genZincblende();
	}
	inline ll_cell genRandom(const size_t D) noexcept {
		if (D!=3) {
			switch (lm__::test::genRndST(0,4)) {
				case 0: return genCubic(D);
				case 1: return genFCC(D);
				case 2: return genBCC(D);
				case 3: return genZincblende(D);
				case 4: return genRndCell(D);
				default: return genZincblende(D);
			}
		} else {
			return genRandom();
		}
	}
}
}

#endif // _LL_TESTTOOLS_
