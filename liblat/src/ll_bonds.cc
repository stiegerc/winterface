// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "ll_bonds.h"

using namespace lm__;
using namespace ll__;

template<class ITYPE>
ll_bonds<ITYPE>::ll_bonds(ll_cell cell, const fMat& R,
		const std::function<bool(const fMat&,const i_i&)>& keep) noexcept:
			cell_(std::move(cell)) {
	assert(R==round(R));

	this->vec_.reserve(cell.N()*cell.N());
		
	// find valid bonds
	for (auto p1=this->cell().ccBegin(),e=this->cell().ccEnd(); p1!=e; ++p1)
		for (auto p2=this->cell().ccBegin(); p2!=e; ++p2) {

			const auto b = *p2 - *p1;
			this->vec_.push_back(i_i_R(size_t(p1),size_t(p2),dim()));
			
			this->vec_.back().reserve(R.N());
			for (auto r=R.ccBegin(),re=R.ccEnd(); r!=re; ++r)
				if (keep(this->cell().B().prod(b+*r),this->vec_.back()))
					this->vec_.back().push_back(*r);

			if (!this->vec_.back().empty())
				this->vec_.back().shrink_to_fit();
			else this->vec_.pop_back();
		}

	this->vec_.shrink_to_fit();
}
template ll_bonds<i_i_R>::ll_bonds(ll_cell cell, const fMat& R,
		const std::function<bool(const fMat&,const i_i&)>& keep) noexcept;
