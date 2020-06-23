// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "lm_tArray.h"
#include "lm_tArray.cc"


// basic properties
template<>
bool lm_tArray<RE__,RE__,CPX__>::cpx() noexcept {
	return false;
}


// printing
template<>
std::string lm_tArray<RE__,RE__,CPX__>::print(const size_t precision, const size_t blanks) const noexcept {
	std::vector<size_t> wmax(N(),0);
	for (size_t n=0; n!=N(); ++n)
		for (size_t m=0; m!=M(); ++m) {
			std::stringstream tmp;
			tmp.precision(precision); tmp.setf(std::ios::right); tmp.setf(std::ios::fixed);

			tmp << (*this)[n*M()+m];
			if (tmp.str().size()>wmax[n]) wmax[n]=tmp.str().size();
		}
	
	std::stringstream outp;
	for (size_t m=0; m!=M(); ++m) {
		outp << std::string(blanks,' ');
		for (size_t n=0; n!=N(); ++n) {
			outp.precision(precision); outp.width(wmax[n]);
			outp.setf(std::ios::right); outp.setf(std::ios::fixed);

			outp << (*this)[n*M()+m];
			outp << (n==N()-1 ? m==M()-1 ? "": "\n": " ");
		}
	}
		
	return outp.str();
}


// instantiation
template class lm_tArray<RE__,RE__,CPX__>;
template lm_tMat<RE__,RE__,CPX__> operator+<RE__,RE__,CPX__>(const RE__ lhs, const lm_tArray<RE__,RE__,CPX__>& rhs) noexcept;
template lm_tMat<CPX__,RE__,CPX__> operator+<RE__,RE__,CPX__>(const CPX__& lhs, const lm_tArray<RE__,RE__,CPX__>& rhs) noexcept;
template lm_tMat<RE__,RE__,CPX__> operator-<RE__,RE__,CPX__>(const RE__ lhs, const lm_tArray<RE__,RE__,CPX__>& rhs) noexcept;
template lm_tMat<CPX__,RE__,CPX__> operator-<RE__,RE__,CPX__>(const CPX__& lhs, const lm_tArray<RE__,RE__,CPX__>& rhs) noexcept;
template lm_tMat<RE__,RE__,CPX__> operator*<RE__,RE__,CPX__>(const RE__ lhs, const lm_tArray<RE__,RE__,CPX__>& rhs) noexcept;
template lm_tMat<CPX__,RE__,CPX__> operator*<RE__,RE__,CPX__>(const CPX__& lhs, const lm_tArray<RE__,RE__,CPX__>& rhs) noexcept;
template lm_tMat<RE__,RE__,CPX__> operator/<RE__,RE__,CPX__>(const RE__ lhs, const lm_tArray<RE__,RE__,CPX__>& rhs) noexcept;
template lm_tMat<CPX__,RE__,CPX__> operator/<RE__,RE__,CPX__>(const CPX__& lhs, const lm_tArray<RE__,RE__,CPX__>& rhs) noexcept;
