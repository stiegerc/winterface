// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_BONDS_
#define _LL_BONDS_

#include "ll_types.h"
#include "ll_compound.h"
#include "ll_cell.h"
#include <iomanip>


namespace ll__ {
	fMat genNNmat(const size_t d) noexcept;
}


namespace ll__ {
	
	//! bond index pair class
	class i_i {
	public:
		/** @name constructors
		 */
		/** constructor from indices
		 * @param i1 first index
		 * @param i2 second index
		 */
		i_i(const size_t i1=NPOS__, const size_t i2=NPOS__) noexcept:
			i1_(i1), i2_(i2) {}

		/** @name information
		 */
		//! first index
		inline size_t i1() const noexcept { return i1_; }
		//! second index
		inline size_t i2() const noexcept { return i2_; }
		
		/** @name comparison
		 */
		//! index equality
		inline bool operator==(const i_i& rhs) const noexcept {
			return i1()==rhs.i1() && i2()==rhs.i2();
		}
		//! index lexicographical order operator<
		inline bool operator<(const i_i& rhs) const noexcept {
			return (i1()==rhs.i1()) ? i2()<rhs.i2(): i1()<rhs.i1();
		}
		//! derived comparison
		inline bool operator>(const i_i& rhs) const noexcept {
			return rhs<*this;
		}
		//! derived comparison
		inline bool operator<=(const i_i& rhs) const noexcept {
			return (*this==rhs)||(*this<rhs);
		}
		//! derived comparison
		inline bool operator>=(const i_i& rhs) const noexcept {
			return (rhs==*this)||(rhs<*this);
		}
		//! derived comparison
		inline bool operator!=(const i_i& rhs) const noexcept {
			return !(*this==rhs);
		}
	
	private:
		/** @name member variables
		 */
		size_t i1_;	//!< first index
		size_t i2_;	//!< second index
	
	public:
		/** @name friends
		 */
		//! streaming operator
		friend inline std::ostream& operator<<(std::ostream& os, const i_i& inp) noexcept {
			return (os<<"("<<inp.i1()<<":"<<inp.i2()<<")");
		}
		//! conjugate, i.e. flip indices
		friend inline i_i conj(const i_i& inp) noexcept {
			return i_i(inp.i2(),inp.i1());
		}
	};
	

	/** index pair including R vectors class
	 */
	class i_i_R final: public i_i, public mat_b<fMat,fArray> {
	public:
		/** @name constructors
		 */
		/** default constructor
		 * @param d_ dimension of R, i.e. of space
		 */
		i_i_R(const size_t d_=DIM__) noexcept: i_i_R(NPOS__,NPOS__,d_) {}
		/** constructor from indices
		 * @param i1 first index
		 * @param i2 second index
		 * @param d dimension of R, i.e. of space
		 */
		i_i_R(const size_t i1, const size_t i2, const size_t d=DIM__) noexcept:
			i_i_R(i1,i2,fMat(d,0)) {}
		/** constructor from indices and R vectors
		 * @param i1 first index
		 * @param i2 second index
		 * @param R R vectors
		 */
		i_i_R(const size_t i1, const size_t i2, fMat R) noexcept:
			i_i(i1,i2), mat_b(std::move(R)) {}


		/** @name information
		 */
		//! R vectors
		inline const fMat& R() const noexcept { return this->mat_; }


		/** friends
		 */
		//! streaming operator
		friend inline std::ostream& operator<<(std::ostream& os, const i_i_R& inp) noexcept {
			return (os<<"("<<inp.i1()<<","<<inp.i2()<<")\nR:\n"<<lm__::T(inp.R()));
		}
		//! conjugate, i.e. flips indices, invert R
		friend inline i_i_R conj(const i_i_R& inp) noexcept {
			return i_i_R(inp.i2(),inp.i1(),-inp.R());
		}
	};


	/** Simple bonds and indices class.
	 * This class holds the bonds in cartesian.
	 * As such it does not contain an ll_cell.
	 */
	class b_I: public mat_vec_b<fMat,fArray,std::vector<i_i>,i_i> {
	public:
		/** @name contructors
		 */
		/** constructor from bonds and index pairs
		 * @param bonds matrix of size DIMxNb containing bonds
		 * @param I index pairs corresponding to bonds
		 */
		inline b_I(fMat bonds, std::vector<i_i> I) noexcept: 
			mat_vec_b<fMat,fArray,std::vector<i_i>,i_i>(std::move(bonds),std::move(I)) {

			if (empty()) return;

			// sort bonds
			const auto J = aux::sorted_order(ccBegin(),ccEnd(),lm__::vcmp);
			aux::reorder(mat_.cBegin(),J);
			aux::reorder(vec_.begin(),J);
				
			// sort indices in equal bonds
			auto mi=mat_.cBegin(),mj=mat_.cBegin()+1,me=mat_.cEnd();
			auto vi=vec_.begin(), vj=vec_.begin()+1;
			for (; mi!=me; mi=mj,++mj,vi=vj,++vj) {
					while (mj!=me && *mi==*mj) ++mj,++vj;
					const auto J = aux::sorted_order(vi,vj);
					aux::reorder(mi,J);
					aux::reorder(vi,J);
				}
		}


		/** @name data access
		 */
		//! bond vectors
		inline const fMat& bonds() const noexcept { return mat_; }
		//! bond indices
		inline const std::vector<i_i>& I() const noexcept { return vec_; }


		/** @name printing
		 */
		/** print bonds to stream
		 * @param os stream to print to
		 * @param prec floating point precision
		 */
		inline std::ostream& print(std::ostream& os, const size_t prec=PPREC__) const noexcept {
			const std::string s = lm__::T(mat_).print(prec);
			auto j=ccBegin(); size_t p1=0, p2=s.find('\n');
			for (auto i=cbegin(),ie=cend(); i!=ie; ++i,++j,p1=p2+1,p2=s.find('\n',p1))
				os << s.substr(p1,p2-p1) << " | "
				   << std::fixed << std::setprecision(prec) << std::setw(prec+4)
				   << lm__::norm(*j) << ", " << *i << "\n";
			return os;
		}
		//! streaming operator
		friend inline std::ostream& operator<<(std::ostream& os, const b_I& inp) noexcept {
			inp.print(os);
			return os;
		}
	};
}



/** Central bonds class. ll_bonds works in tandem with ll_cell, since the start and
 * ending indices correspond to the contents of a unit cell. A bond is composed of
 * - a starting index i1
 * - an ending index i2
 * - a bond vector, decomposed into the 'principal bond' component,
 *   		    i.e. a bond inside the unit cell, and an R vector
 *   		    component pointing to the target unit cell
 * The data is totally sorted. First by {i1,i2} lexicographically and for each
 * index pair, the R vectors are sorted lexicographically.
 * An important property is whether the bonds are (inversion) symmetric
 * or not, i.e. whether for each {i1,i2,b}, {i2,i1,-b} is included, especially
 * when interactions are included (in ll_hbonds).
 */
template <class ITYPE=ll__::i_i_R>
class ll_bonds: public ll__::vec_cb<std::vector<ITYPE>,ITYPE> {
public:
	/** @name types
	 */
	typedef lm__::fArray fArray;		//!< real array
	typedef lm__::fCol fCol;		//!< real column
	typedef lm__::fMat fMat;		//!< real matrix
	typedef ll__::idv idv;			//!< identity string vector
	typedef ll__::i_i i_i;			//!< index pair
	typedef ll__::i_i_R i_i_R;		//!< index pair with R vectors
	typedef ll__::b_I b_I;			//!< simple bonds
	typedef ll__::c_fColItr c_fColItr;	//!< column const_iterator
	typedef typename std::vector<ITYPE>::
		const_iterator vitr;		//!< iterator to bond entry

public:
	/** @name constructors
	 */
	//! default constructor
	inline ll_bonds<ITYPE>() noexcept: cell_() {}
	//! constructor from unit cell
	inline ll_bonds<ITYPE>(ll_cell cell) noexcept: cell_(std::move(cell)) {}
	//! constructor from unit cell and bond entries
	inline ll_bonds<ITYPE>(ll_cell cell, std::vector<ITYPE> dat) noexcept:
		ll__::vec_cb<std::vector<ITYPE>,ITYPE>(std::move(dat)), cell_(std::move(cell)) {}
	/** constructor from cell and R vectors
	 * @param cell the underlying unit cell
	 * @param R R vectors for bonds
	 * @param keep lambda to decide which bonds to keep
	 */
	ll_bonds<ITYPE>(ll_cell cell, const fMat& R,
		const std::function<bool(const fMat&, const i_i&)>& keep) noexcept;
	/** constructor from cell and R vectors
	 * @param cell the underlying unit cell
	 * @param R R vectors for bonds
	 * @param r cutoff for bond length above which bonds are discarded
	 */
	inline ll_bonds<ITYPE>(const ll_cell& cell, const fMat& R, const double r) noexcept:
		ll_bonds<ITYPE>(cell,R,[r](const fMat& b,const i_i& j)->bool{
			const double bn = lm__::norm(b);
			return lm__::ops::leq(bn,r) && lm__::ops::nz(bn);
		}) {}
	/** constructor from two sets of index pairs.
	 * Includes only bonds of length <= r and going from the set of
	 * indices I1 -> I2 or I2 -> I1, i.e. not I1(+)I2 -> I1(+)I2
	 * @param cell the underlying unit cell
	 * @param R R vectors for bonds
	 * @param r cutoff for bond length above which bonds are discarded
	 * @param I1 first set of index pairs
	 * @param I2 second set of index pairs
	 */
	inline ll_bonds<ITYPE>(const ll_cell& cell, const fMat& R, const double r,
			const std::vector<size_t>& I1, const std::vector<size_t>& I2) noexcept:
		ll_bonds<ITYPE>(cell,R,[r,&I1,&I2](const fMat& b,const i_i& j)->bool{
			assert(std::is_sorted(I1.cbegin(),I1.cend()));
			assert(std::is_sorted(I2.cbegin(),I2.cend()));

			// bond length check
			const double bn = lm__::norm(b);
			if (lm__::ops::z(bn) || lm__::ops::gt(bn,r))
				return false;

			// type check bond from I1 to I2
			{
				const auto itr1 = std::lower_bound(I1.cbegin(),I1.cend(),j.i1());
				const auto itr2 = std::lower_bound(I2.cbegin(),I2.cend(),j.i2());
				if (itr1!=I1.cend() && *itr1==j.i1() && itr2!=I2.cend() && *itr2==j.i2())
					return true;
			}
			// type check bond from I2 to I1
			{
				const auto itr1 = std::lower_bound(I1.cbegin(),I1.cend(),j.i2());
				const auto itr2 = std::lower_bound(I2.cbegin(),I2.cend(),j.i1());
				if (itr1!=I1.cend() && *itr1==j.i2() && itr2!=I2.cend() && *itr2==j.i1())
					return true;
			}
			
			// bad index pair
			return false;
		}) {}


	/** @name iterators
	 */
	using ll__::vec_cb<std::vector<ITYPE>,ITYPE>::begin;
	using ll__::vec_cb<std::vector<ITYPE>,ITYPE>::end;
	using ll__::vec_cb<std::vector<ITYPE>,ITYPE>::cbegin;
	using ll__::vec_cb<std::vector<ITYPE>,ITYPE>::cend;
	//! iterator to beginning of the range of index i
	inline auto begin(const size_t i) noexcept -> decltype(this->cbegin()) {
		return std::lower_bound(this->begin(),this->end(),i,
			[](const auto& j, const auto i)->bool{return j.i1()<i;});
	}
	//! iterator past the end of the range of index i
	inline auto end(const size_t i) noexcept -> decltype(this->cend()) {
		return (i==NPOS__) ? this->end():
			std::lower_bound(this->begin(),this->end(),i+1,
			[](const auto& j, const auto i)->bool{return j.i1()<i;});
	}
	//! const_iterator to beginning of the range of index i
	inline auto begin(const size_t i) const noexcept -> decltype(this->cbegin()) {
		return std::lower_bound(this->cbegin(),this->cend(),i,
			[](const auto& j, const auto i)->bool{return j.i1()<i;});
	}
	//! const_iterator past the end of the range of index i
	inline auto end(const size_t i) const noexcept -> decltype(this->cend()) {
		return (i==NPOS__) ? this->cend():
			std::lower_bound(this->cbegin(),this->cend(),i+1,
			[](const auto& j, const auto i)->bool{return j.i1()<i;});
	}
	//! const_iterator to beginning of the range of index i
	inline auto cbegin(const size_t i) const noexcept -> decltype(this->cbegin()) {
		return this->begin(i);
	}
	//! const_iterator past the end of the range of index i
	inline auto cend(const size_t i) const noexcept -> decltype(this->cend()) {
		return this->end(i);
	}


	/** @name data access
	 */
	//! underlying unit cell
	inline const ll_cell& cell() const noexcept { return cell_; }
	//! dimension of space
	inline size_t dim() const noexcept { return cell().dim(); }
	//! iterator past the end of the unit cell
	const c_fColItr e() const noexcept { return cell_.ccEnd(); }
	

	/** @name information
	 */
	//! checks whether there are no bonds
	inline bool empty() const noexcept { return !this->cardinality(); }
	//! returns the maximum spacial tolerance for which bonds are distinguishable (in direct)
	inline static double maxtol() noexcept { return .499; }
	/** returns a restriction matrix for space. A spacial dimension is considered restricted if
	 * none of the R vectors have entries greater other than 0 in it.
	 */
	inline fMat rmat() const noexcept {
		return std::accumulate(cbegin(),cend(),lm__::ones<fMat>(dim(),1),
			[](const fMat& s, const auto& i) -> fMat { return s & ~lm__::nany(i.R()); });
	}
	/** returns a restriction vector for space. A spacial dimension is considered restricted if
	 * none of the R vectors have entries greater other than 0 in it.
	 */
	inline ll__::rv r() const noexcept {
		const auto r = this->rmat();
		ll__::rv res(r.size());
		std::transform(r.cbegin(),r.cend(),res.begin(),[](const double i)->bool{return i;});
		return res;
	}
	/** construct the principal bond for two indices
	 * @param i1 first index
	 * @param i2 second index
	 */
	inline fMat bond(const size_t i1, const size_t i2) const noexcept {
		return cell().Ap().cAt(i2) - cell().Ap().cAt(i1);
	}
	/** constuct the principal bond for an index pair
	 * @param I the index pair {i1,i2}
	 */
	inline fMat bond(const i_i& I) const noexcept {
		return cell().Ap().cAt(I.i2()) - cell().Ap().cAt(I.i1());
	}
	//! returns the tolerance level used for querying bonds
	inline double queryTol() const noexcept { return tol_; }
	//! total number bonds in the class
	inline size_t cardinality() const noexcept {
		return std::accumulate(this->cbegin(),this->cend(),size_t(0),
			[](const size_t s, const auto& i)->size_t{return s+i.N();});
	}
	//! returns the indices for which there exist bonds
	inline std::vector<size_t> inds() const noexcept {
		std::vector<size_t> res; res.reserve(2*this->size());
		for (const auto& i: *this)
			res.push_back(i.i1()), res.push_back(i.i2());
		std::sort(res.begin(),res.end());
		res.resize(std::distance(res.begin(),std::unique(res.begin(),res.end())));
		return res;
	}
	/** returns the index for an id string
	 * @param s 'atomic' id string, e.g. 'Mo'
	 */
	inline size_t ind(const std::string& s) const noexcept {
		const auto itr = std::lower_bound(cell().id().cbegin(),cell().id().cend(),s);
		return itr!=cell().id().cend() && s==*itr ?
			std::distance(cell().id().cbegin(),itr): NPOS__;
	}
	/** returns the indices for a range of id strings
	 * @param id vector of 'atomic' id strings, e.g. {'Mo','S'}
	 */
	inline std::vector<size_t> ind(const idv& id) const noexcept {
		std::vector<size_t> res; res.reserve(id.size());
		for (const auto& s: id)
			res.push_back(ind(s));
		return res;
	}
	//! returns the id string for an index
	inline const std::string& id(const size_t i) const noexcept {
		return cell().id(i);
	}
	//! returns the is strings for a range of indices
	inline idv id(const std::vector<size_t>& I) const noexcept {
		return cell().id(I);
	}
	/** returns the number of bonds starting or ending with an index.
	 * Each occurrence is counted, i.e. both directions inverse bond pairs
	 * and twice if the bond starts and ends at the same index. The number is
	 * therefore always an even number.
	 */
	inline size_t Nindex(const size_t ind) const noexcept {
		size_t res=0;
		std::for_each(this->cbegin(),this->cend(),[&res,ind](const auto& i)->void{
			if (i.i1()==ind) res+=i.N();
			if (i.i2()==ind) res+=i.N();
		});
		return res;
	}
	/** returns the numbers of bonds starting or ending with an index of a range of indices.
	 * Each occurrence is counted, i.e. both directions inverse bond pairs
	 * and twice if the bond starts and ends at the same index. The number is
	 * therefore always an even number.
	 */
	inline std::vector<size_t> Nindex(std::vector<size_t> inds) const noexcept {
		std::transform(inds.cbegin(),inds.cend(),inds.begin(),[this](const size_t i)->size_t{
			return Nindex(i);
		});
		return inds;
	}
	/** returns an upper bound for the bond length.
	 * e.g. all bonds fit within a sphere of this radius.
	 * Does not considered starting an ending points of points, just the bond vectors.
	 */
	inline double radius() const noexcept {
		double res = 0.0;
		if (this->empty()) return res;

		// strict radius from bond vectors
		for (auto i=this->cbegin(),ie=this->cend(); i!=ie; ++i) {
			const auto b = this->bond(*i);
			for (auto r=i->ccBegin(),re=i->ccEnd(); r!=re; ++r) {
				const double d = lm__::norm(cell().B().prod(*r+b));
				if (res<d) res=d;
			}
		}

		// add due to tolerance level
		res += std::abs(lm__::mtol())*lm__::norm(lm__::nsum(this->cell().B()));

		return res;
	}
	/** returns the sides of a cuboid encompassing all bonds.
	 * Does not considered starting an ending points of points, just the bond vectors.
	 */
	inline fMat boundary() const noexcept {
		fMat res = lm__::zeros<fMat>(this->dim(),1);
		for (auto i=this->cbegin(),ie=this->cend(); i!=ie; ++i) {
			const auto bnd = this->bond(*i);
			for (auto b=i->ccBegin(),be=i->ccEnd(); b!=be; ++b) {
				const auto m = lm__::nmax(lm__::abs(cell().B().prod(*b+bnd))).mat;
				auto mi = m.cbegin();
				for (auto ri=res.begin(),re=res.end(); ri!=re; ++ri,++mi)
					if (*ri<*mi) *ri=*mi;
			}
		}

		// add due to tolerance level
		res += std::abs(lm__::mtol())*lm__::nmax(lm__::abs(this->cell().B().prod(
				ll__::genNNmat(this->dim())))).mat;
		
		return res;
	}
	//! returns unique R vectors, i.e. the range of unit cells where bonds exist
	inline fMat range() const noexcept {
		fMat res(this->dim(),0); res.reserve(this->cardinality());
		for (auto i=this->cbegin(),ie=this->cend(); i!=ie; ++i)
			for (auto r=i->ccBegin(),re=i->ccEnd(); r!=re; ++r) {
				const auto itr = std::lower_bound(
					res.ccBegin(),res.ccEnd(),*r,lm__::vcmp);
				if (itr==res.ccEnd() || *itr!=*r)
					res.cInsert(itr,*r);
			}
		res.shrink_to_fit();
		return res;
	}
	/** returns the unique R vectors, i.e. the range of unit cells where bonds exist
	 * but view in a different basis B
	 * @param B other basis for which the range is computed
	 */
	inline fMat range(const fMat& B) const noexcept {
		assert(B.M()==this->dim() && B.N()==this->dim());
		assert(lm__::rank(B)==this->dim());
		if (this->empty()) return fMat(this->dim(),0);

		const fMat TRM = B.leftDivide(this->cell().B());
		fMat res(this->dim(),0); res.reserve(this->cardinality());
		
		// zero R is always included
		res.push_back(lm__::zeros<fMat>(this->dim(),1));
		
		// check end points of bonds in other basis B
		for (auto i=this->cbegin(),ie=this->cend(); i!=ie; ++i) {
			const fCol& p2 = this->cell().Ap().cAt(i->i2());
			
			for (auto r=i->ccBegin(),re=i->ccEnd(); r!=re; ++r) {
				const fMat p2_ = lm__::floor(TRM.prod(p2 + *r));
				
				// insert +
				{
					const auto itr = std::lower_bound(
						res.ccBegin(),res.ccEnd(),p2_,lm__::vcmp);
					if (itr==res.ccEnd() || *itr!=p2_)
						res.cInsert(itr,p2_);
				}
				// insert -
				{
					const auto itr = std::lower_bound(
						res.ccBegin(),res.ccEnd(),-p2_,lm__::vcmp);
					if (itr==res.ccEnd() || *itr!=-p2_)
						res.cInsert(itr,-p2_);
				}
			}
		}
		res.shrink_to_fit();
		return res;
	}
	/** returns the neighborhood for a bond index. i.e. all the bonds starting from this index
	 * @param ind starting index
	 * @param Rcut cutoff radius for the bonds
	 * @param sptol spacial tolerance used when sorting the bonds
	 */
	inline b_I neighborhood(const size_t ind, const double Rcut, const double sptol) const noexcept {
		assert(ind<cell().N());
		
		std::vector<i_i> I; I.reserve(this->Nindex(ind));
		fMat b(this->dim(),0); b.reserve(I.capacity());
		
		for (auto i=this->cbegin(ind),ie=this->cend(ind); i!=ie; ++i) {
			
			auto pb = this->bond(i->i1(),i->i2());
			for (auto j=i->R().ccBegin(),je=i->R().ccEnd(); j!=je; ++j) {
				const auto cb = this->cell().B().prod(pb + *j);
				if (lm__::norm(cb)<Rcut)
					b.push_back(cb),
					I.push_back({i->i1(),i->i2()});
			}
		}

		lm__::set_mtol(sptol);
		const auto J = aux::sorted_order(b.ccBegin(),b.ccEnd(),lm__::vcmp);
		lm__::reset_mtol();
		aux::reorder(b.cBegin(),J);
		aux::reorder(I.begin(),J);

		return {std::move(b),std::move(I)};
	}
	/** checks whether the bonds are symmetric. i.e. whether for each {i1,i2,R},
	 * {i2,i1,-R} is also included.
	 */
	inline bool symmetric() const noexcept {
		for (const auto& j: *this) {
			// find reverse index pair
			const auto jtr = std::lower_bound(
				this->cbegin(),this->cend(),i_i(j.i2(),j.i1()));
			if (jtr==this->cend() || jtr->i2()!=j.i1() || jtr->i1()!=j.i2())
				return false;
			
			// check all inverted R are present
			for (auto i=j.ccBegin(),e=j.ccEnd(); i!=e; ++i)
				if (!std::binary_search(jtr->ccBegin(),jtr->ccEnd(),-*i))
					return false;
		}
		return true;
	}
	

	/** @name modification
	 */
	/** rotate the bonds
	 * @param R rotation matrix
	 */
	inline void rotate(const fMat& R) noexcept { cell_.rotate(R); }
	/** scale the bonds
	 * @param f scaling factor
	 */
	inline void scale(const double f) noexcept { cell_.scale(f); }
	/** scale the bonds
	 * @param f array of scaling factors, one for each dimension
	 */
	inline void scale(const fArray& f) noexcept { cell_.scale(f); }
	/** stress the bonds
	 * @param S stress tensor
	 */
	inline void stress(const fMat& S) noexcept { cell_.stress(S); }


	/** @name conversion
	 */
	//! return the bonds in simple format
	inline b_I simple() const noexcept {
		fMat bnd(this->dim(),0); bnd.reserve(cardinality());
		std::vector<i_i> I; I.reserve(bnd.ccap());
		
		for (auto i=this->cbegin(),ie=this->cend(); i!=ie; ++i) {
			const auto b = this->bond(*i);
			for (auto n=i->ccBegin(),ne=i->ccEnd(); n!=ne; ++n) {
				bnd.push_back(cell().B().prod(*n+b));
				I.push_back(i_i(i->i1(),i->i2()));
			}
		}

		return b_I(std::move(bnd),std::move(I));
	}
	//! return the bonds starting at index i in simple format
	inline b_I simple(const size_t i) const noexcept {
		fMat bnd(this->dim(),0); bnd.reserve(Nindex(i));
		std::vector<i_i> I; I.reserve(bnd.ccap());

		for (auto j=this->cbegin(i), je=this->cend(i); j!=je; ++j) {
			const auto b = this->bond(*j);
			for (auto n=j->ccBegin(),ne=j->ccEnd(); n!=ne; ++n) {
				bnd.push_back(cell().B().prod(*n+b));
				I.push_back(i_i(j->i1(),j->i2()));
			}
		}
		
		return b_I(std::move(bnd),std::move(I));
	}
	/** find the centers of bonds
	 * @return a unit cell holding the bond centers as 'atoms'
	 */
	inline ll_cell getBondCenters() const noexcept {
		if (this->empty()) return ll_cell(cell().B());

		fMat bc(dim(),0); bc.reserve(this->cardinality()/2);
		ll__::idv bid; bid.reserve(bc.ccap());
		const std::regex r_("[_]+");

		// find centers
		for (auto i=this->cbegin(),ie=this->cend(); i!=ie; ++i) {
			if (i->i1()>i->i2()) continue; // skip inverted partners
			const auto cpbc = cell().cAt(i->i1()) + this->bond(*i)*.5;

			// if i1==i2 inverted partners are among the 2nd half of R vectors
			for (auto r=i->ccBegin(),re=(i->i1()==i->i2() ? i->ccBegin()+i->N()/2: i->ccEnd());
					r!=re; ++r) {
				const auto cbc = (cpbc + *r*.5)%1.0;
				const auto itr = std::lower_bound(bc.ccBegin(),bc.ccEnd(),cbc);
				if (itr!=bc.ccEnd() && *itr==cbc) {
					
					// center already included, dump all to buff vector
					std::vector<std::string> buff;
					for (std::sregex_token_iterator
						i(bid[(size_t)itr].cbegin(),bid[(size_t)itr].cend(),r_,-1), e;
						i!=e; ++i) buff.push_back(*i);
					
					const std::string id1 = cell().id(cell().type(i->i1()));
					const std::string id2 = cell().id(cell().type(i->i2()));
					buff.push_back(id1<id2 ? '('+id1+':'+id2+')': '('+id2+':'+id1+')');

					// sort, remove duplicates and reconstruct
					std::sort(buff.begin(),buff.end());
					buff.resize(distance(buff.begin(),
						 std::unique(buff.begin(),buff.end())));
					bid[(size_t)itr] = "";
					for (const auto& s: buff)
						bid[(size_t)itr] += (s+"_");
					bid[(size_t)itr].pop_back();

				} else {
					// new center
					bc.cInsert(itr,cbc);
					
					const std::string id1 = cell().id(cell().type(i->i1()));
					const std::string id2 = cell().id(cell().type(i->i2()));
					bid.insert(bid.cbegin()+(size_t)itr,
						id1<id2 ? '('+id1+':'+id2+')': '('+id2+':'+id1+')');
				}
			}
		}

		// add parenthesis around ids with '_' in them
		for (auto& s: bid)
			if (s.find('_')!=std::string::npos) s = '('+s+')';

		// return
		return ll_cell(cell().B(),std::move(bc),std::move(bid));
	}
	

	/** @name searching
	 */
	//! set bond query tolerance in direct coordinates
	inline void setQueryTolDirect(const double tol) const noexcept {
		assert(tol>=.0 && tol<.5); tol_ = tol;
	}
	//! set bond query tolerance in cartesian coordinates
	inline void setQueryTolCartesian(const double tol) const noexcept {
		setQueryTolDirect(cell().directTol(tol));
	}
	/** search for a bond index pair in the data
	 * @return iterator to the matching index pair
	 */
	inline vitr search(const i_i& j) const noexcept {
		const auto itr = std::lower_bound(this->cbegin(),this->cend(),j);
		return itr!=this->cend() && *itr==j ? itr: this->cend();
	}
	/** search for an R vector belonging to an index pair and a bond vector
	 * @param itr iterator to the matching index pair
	 * @param b bond vector in cartesian coordinates
	 * @return a column iterator pointing to the matching R vector
	 */
	inline c_fColItr search(const vitr& itr, const fArray& b) const noexcept {
		assert(b.M()==this->dim());
		assert(b.N()==1);
		if (itr==this->cend()) return e();

		const auto b_ = this->cell().B().leftDivide(b) - this->bond(*itr);
		lm__::set_mtol(tol_);
		const auto jtr = std::lower_bound(itr->ccBegin(),itr->ccEnd(),b_);
		if (jtr!=itr->ccEnd() && *jtr==b_) {
			lm__::reset_mtol();
			return jtr;
		}  else {
			lm__::reset_mtol();
			return e();
		}
	}
	

	/** @name printing
	 */
	/** print bonds to stream
	 * @param os the stream to print into
	 * @param mode printing mode
	 */
	inline std::ostream& print(std::ostream& os, const size_t mode=0) const noexcept {
		switch (mode) {
		case 0:
		{
			for (const auto& i: *this)
				os << i << "\nb:" << lm__::T(this->bond(i)) << "\n";
		}
		break;
		case 1:
		{
			os << simple();
		}
		break;
		}
		return os;
	}
	/** print bonds to file
	 * @param fileName name of the file to print into
	 * @param mode printing mode
	 * @param prec floating point precision
	 */
	inline void printToFile(const std::string& fileName, const size_t mode=0, const size_t prec=PPREC__) const {

		std::ofstream file;
		file.open(fileName);
		if (!file.good()) throw(std::invalid_argument("failed to open '"+fileName+"'"));
		
		print(file,mode,prec);
		file.close();
	}
	//! streaming operator
	friend inline std::ostream& operator<<(std::ostream& os, const ll_bonds& inp) noexcept {
		return inp.print(os);
	}


protected:
	/** @name member variables
	 */
	ll_cell cell_;			//!< underlying unit cell
	mutable double tol_ = WTOL__;	//!< querying tolerance
};

#endif // _LL_BONDS_

/** @}
 */
