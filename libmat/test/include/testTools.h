// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _TESTTOOLS_
#define _TESTTOOLS_

#include "lm_defs.h"
#include "lm_tMat.h"
#include "lm_types.h"
#include "lm_fn.h"
#include <cmath>
#include <cassert>
#include <cstring>
#include <random>
#include <algorithm>
#include <iostream>


// to_string for complex
namespace std {
	inline std::string to_string(const std::complex<RE__>& inp) {
		return "("+std::to_string(inp.real())+","+std::to_string(inp.imag())+")";
	}
}

// compare macro with message
#ifndef CMP_DELTA
#define CPPUNIT_ASSERT_DELTA(ref,act,delta) \
	CPPUNIT_ASSERT_MESSAGE("ref: "+std::to_string(ref)+", act: "+std::to_string(act),std::abs(ref-act)<delta);
#endif


namespace lm__ {
namespace test {
	// random size_t
	struct dims {
		size_t M;
		size_t N;
	};
	inline size_t genRndST(const size_t l=10, const size_t u=20) noexcept {
		assert(l<=u);
	
		std::random_device rd;
		std::mt19937_64 gen(rd());
		std::uniform_int_distribution<size_t> dis(l,u);
	
		return dis(gen);
	}
	inline dims genRndMEST(const size_t l=10, const size_t u=20) noexcept {
		dims res;
		do {
			res.M = genRndST(l,u);
			res.N = genRndST(l,u);
		} while (res.M==res.N);
		return res;
	}

	// random range
	inline void rnd(RE__* ptr, const size_t L, const RE__ l=0.0, const RE__ u=1.0) noexcept {
		assert(l<=u);
		std::random_device rd;
		std::mt19937_64 gen(rd());
		std::uniform_real_distribution<RE__> dis(u,l);
	
		std::for_each(ptr,ptr+L,[&](RE__& i){i=dis(gen);});
	}
	inline void rnd(CPX__* ptr, const size_t L, const RE__ l=0.0, const RE__ u=1.0) noexcept {
		assert(l<=u);

		std::random_device rd;
		std::mt19937_64 gen(rd());
		std::uniform_real_distribution<RE__> dis(u,l);
	
		auto ptr_ = reinterpret_cast<RE__*>(ptr);
		std::for_each(ptr_,ptr_+2*L,[&](RE__& i){i=dis(gen);});
	}

	// random matrix
	template<class MT>
	MT rnd(const size_t M, const size_t N, const RE__ l=0.0, const RE__ u=1.0) noexcept {
		MT res(M,N);
		rnd(res.data(),res.L(),l,u);
		return res;
	}
	template<class MT>
	MT rnd_b(const size_t M, const size_t N, const RE__ l=0.0, const RE__ u=1.0) noexcept {
		MT res(M,N);
		
		const size_t R = std::min(M,N);
		do rnd(res.data(),res.L(),l,u);
		while (rank(res)<R);
		return res;
	}
	template<class MT>
	MT rnd_sqb(const size_t M, const RE__ l=-1.0, const RE__ u=1.0) noexcept {
		MT res(M,M);
		do rnd(res.data(),res.L(),l,u);
		while (std::abs(lm__::det(res))<std::abs(std::pow((u-l)/RE__(3.0),M))
			&& std::abs(lm__::det(res))<mtol());
		return res;
	}

	// generate random index vector
	inline std::vector<size_t> genRndI(const size_t sl, const size_t su,
					   const size_t il, const size_t iu) noexcept {
		assert(su>=sl);
		assert(iu>=il);
		assert((iu-il)>=(su-sl));

		std::vector<size_t> res; res.reserve(genRndST(sl,su));
		std::vector<bool> ck(iu-il+1,false);
		while (res.size()<res.capacity()) {
			const size_t t = genRndST(il,iu);
			if (!ck[t-il]) {
				res.push_back(t);
				ck[t-il] = true;
			}
		}
	
		return res;
	}

	// conversion to TT
	template<class TT>
	inline TT cc(const RE__ inp) noexcept {
		return TT(inp);
	}
	template<class TT>
	inline TT cc(const CPX__& inp) noexcept {
		return TT(inp);
	}
	template<>
	inline RE__ cc<RE__>(const CPX__& inp) noexcept {
		return std::real(inp);
	}

}
}

#endif // _TESTTOOLS_
