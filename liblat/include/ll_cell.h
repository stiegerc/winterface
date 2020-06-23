// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_CELL_
#define _LL_CELL_

#include "ll_defs.h"
#include "ll_types.h"
#include "ll_compound.h"
#include "lm_fn.h"
#include "ll_lambda.h"
#include <utility>
#include <cmath>
#include <cassert>
#include <vector>
#include <string>
#include <regex>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <set>



#ifndef NLZ__
#define NLZ__ 2		//!< number of leading zeros when indexing
#endif


namespace ll__ {
	class i_i;
	class i_i_R;
}
template<class ITYPE> class ll_bonds;


/** This class defines a unit cell of an (atomic) lattice.
 * By convention, all coordinates of atomic positions are expressed in
 * a basis B, called direct coordinates, and inside the interval [0,1)^N.
 * Additionally, the data is completely ordered thus guaranteeing equal
 * state for equal input (upto a permutation). The ordering is first done
 * along id strings (e.g. 'Mo') and then along coordinates, in both cases
 * lexicographically. If two distinguished atomic types are given the same
 * id string, they are indexed, e.g. 'Mo-001', 'Mo-002', ... \n
 * However, id strings are recommended but optional. Furthermore, id strings
 * can be composited from other id strings such as Mo-001,S-002 -> (Mo-001:S-002)
 * which is the id strings given to the bond center between these two atoms.
 * In principle, this can go on indefinitely, e.g. ((Mo-001:S-002):S-001) etc.
 * Of particular importance is the changeBasis or expandBasis function, which
 * expresses the same lattice in terms of a different periodicity as expressed
 * in the basis B. The new and old periodicity do not need to be of the same
 * symmetry group.
 */
class ll_cell: public ll__::mat_cb<ll__::fMat,ll__::fArray> {
public:
	/** @name types
	 */
	typedef ll__::fArray fArray;			//!< real Array, fCol, fRow or fMat
	typedef ll__::fCol fCol;			//!< real column
	typedef ll__::fMat fMat;			//!< real matrix
	typedef ll__::aT aT;				//!< atomic type
	typedef ll__::aTv aTv;				//!< atomic type vector
	typedef ll__::aC aC;				//!< atomic count
	typedef ll__::aCv aCv;				//!< atomic count vector
	typedef ll__::rv rv;				//!< restriction vector
	typedef ll__::idv idv;				//!< identity string vector
	typedef ll__::i_i i_i;				//!< index pair
	typedef ll__::i_i_R i_i_R;			//!< index pair and R vectors


public:
	/** @name contructors
	 */
	//! default constructor
	inline ll_cell() noexcept {}
	/** constructor from basis only
	 * @param B basis
	 */
	explicit ll_cell(fMat B) noexcept: B_(std::move(B)) {}
	/** constructor from atomic counts.
	 * This constructor will index id strings if duplicates are found,
	 * since the type is determined by the atomic count range.
	 * @param B basis
	 * @param Ap atomic positions, must be in basis B and in [0,1)
	 * @param N atomic count vector, need sum(N) = size(Ap,2)
	 * @param id optional id strings
	 */
	explicit ll_cell(fMat B, fMat Ap, aCv N, idv id={}) noexcept;
	/** constructor from id string vector.
	 * This constructor will consider duplicate id strings as of
	 * the same atomic type and thus group positions into types
	 * rather than index id strings.
	 * @param B basis
	 * @param Ap atomic positions, must be in basis B and in [0,1)
	 * @param id range of strings
	 */
	explicit ll_cell(fMat B, fMat Ap, idv id) noexcept;
	/** constructor from file.
	 * This constructor can red POSCAR, wannier90 wout files
	 * and OMEN lattice files.
	 * @param fileName name of the file
	 * @param d dimension of space
	 */
	explicit ll_cell(const std::string& fileName, size_t d=DIM__);
	//! copy constructor
	inline ll_cell(const ll_cell& rhs) noexcept:
		mat_cb<fMat,fArray>(rhs.mat_), B_(rhs.B_), T_(rhs.T_), id_(rhs.id_) {}
	//! move constructor
	inline ll_cell(ll_cell&& rhs) noexcept: ll_cell() { swap(*this,rhs); }
	/** constructor from atomic count vector and id string.
	 * This constructor will generate indexed id strings from id for all atomic
	 * types determined by the atomic count vector.
	 * @param B basis
	 * @param Ap atomic positions
	 * @param N atomic count vector
	 * @param id id string
	 */
	inline explicit ll_cell(fMat B, fMat Ap, aCv N, const std::string& id) noexcept {
		idv tmp(N.size()); size_t ind=1;
		std::generate(tmp.begin(),tmp.end(),[sid=softstripId(id),&ind]()
						  {return appendIndex_(sid,ind++);});
		*this = ll_cell(std::move(B),std::move(Ap),std::move(N),std::move(tmp));
	}


	/** @name assignment
	 */
	//! swap function
	inline friend void swap(ll_cell& lhs, ll_cell& rhs) noexcept {
		swap(lhs.B_,rhs.B_);
		swap(lhs.mat_,rhs.mat_);
		swap(lhs.T_,rhs.T_);
		swap(lhs.id_,rhs.id_);
	}
	//! assignment operator
	inline ll_cell& operator=(const ll_cell& rhs) noexcept {
		B_ = rhs.B_;
		mat_ = rhs.mat_;
		T_ = rhs.T_;
		id_ = rhs.id_;
		return *this;
	}
	//! move assignment operator
	inline ll_cell& operator=(ll_cell&& rhs) noexcept {
		swap(*this,rhs); return *this;
	}
	

	/** @name iterators
	 */
	using ll__::mat_cb<fMat,fArray>::cBegin;
	using ll__::mat_cb<fMat,fArray>::ccBegin;
	using ll__::mat_cb<fMat,fArray>::cEnd;
	using ll__::mat_cb<fMat,fArray>::ccEnd;
	/** const_iterator to the beginning of the range of a type
	 * @param t atomic type
	 */
	inline auto cBegin(const aT t) const noexcept {
		return validType(t) ? Ap().cBegin()+T_[t]: Ap().cEnd();
	}
	/** const_iterator past the end of the range of a type
	 * @param t atomic type
	 */
	inline auto cEnd(const aT t) const noexcept {
		return validType(t) ? Ap().cBegin()+T_[t+1]: Ap().cEnd();
	}
	/** const_iterator to the beginning of the range of a type
	 * @param t atomic type
	 */
	inline auto ccBegin(const aT t) const noexcept { return cBegin(t); }
	/** const_iterator past the end of the range of a type
	 * @param t atomic type
	 */
	inline auto ccEnd(const aT t) const noexcept { return cEnd(t); }


	/** @name data access
	 */
	//! basis const reference
	inline const fMat& B() const noexcept { return B_; }
	//! atomic positions const reference
	inline const fMat& Ap() const noexcept { return mat_; }


	/** @name memory management
	 */
	//! basis rvalue reference
	inline fMat&& moveB() noexcept { return std::move(B_); }
	//! atomic positions rvalue reference
	inline fMat&& moveAp() noexcept { return std::move(mat_); }
	//! id range rvalue reference
	inline idv&& moveId() noexcept { return std::move(id_); }


	/** @name information
	 */
	//! returns the number of distinguished atomic types 
	inline size_t Nspecies() const noexcept { return T_.size()-1; }
	//! returns the volume of the cell
	inline double vol() const noexcept { return ll__::vol(B()); }
	//! returns the normalized volume of the cell
	inline double volnorm() const noexcept { return ll__::volnorm(B()); }
	//! returns the sign of the cell, e.g. sign(det(B))
	inline double sign() const noexcept { return empty() ? 1.0: lm__::det(B())<0.0 ? -1.0: 1.0; }
	/** returns an iterator to a position
	 * @param inp position to look for
	 * @param t atomic type
	 * @return column const_iterator
	 */
	lm__::c_fColItr find(const fArray& inp, const aT t) const noexcept;
	/** returns whether a column is a 'lattice' vector.
	 * A lattice vector is a vector by which all positions in the cell
	 * may be shifted and resulting in the same unit cell.
	 * @param inp column to check
	 */
	bool lVec(const fCol& inp) const noexcept;
	/** returns whether the columns of a matrix are 'lattice' vectors.
	 * A lattice vector is a vector by which all positions in the cell
	 * may be shifted and resulting in the same unit cell.
	 * @param inp matrix to check
	 */
	bool lVec(const fMat& inp) const noexcept;
	//! returns whether the unit cell is primitive
	bool primitive() const noexcept;
	/** checks whether a basis holds a valid periodicity
	 * @param rhs basis to check
	 */
	bool validBasis(fMat rhs) const noexcept;
	/** attempts to determine spacial restrictions by looking
	 * for unusally long columns in B and large gaps in positions
	 * by comparing to the bond length between nn atoms.
	 * @param f factor for the bond length of nn positions
	 */
	rv r(const double f) const noexcept;
	/** returns a tolerance level specified for cartesian coordinates
	 * in direct coordinates.
	 * @param tol tolerance level
	 */
	double directTol(const double tol) const noexcept;


	/** @name valid types and indices
	 */
	//! returns the valid atomic types
	inline aTv types() const noexcept { return ll__::rg(0,1,Nspecies()); }
	/** returns a range of indices where atoms of a specific type are located.
	 * @param t atomic type
	 */
	inline std::vector<size_t> inds(const aT t) const noexcept {
		if (!validType(t)) return {};
		std::vector<size_t> res(Ntype(t));
		std::iota(res.begin(),res.end(),T_[t]);
		return res;
	}
	/** returns a range of indices where atoms of specific types are located.
	 * @param T range of atomic types
	 */
	inline std::vector<size_t> inds(const aTv& T) const noexcept {
		std::vector<size_t> res; res.reserve(
			 std::accumulate(T.cbegin(),T.cend(),size_t(0),
			[this](const size_t s, const aT t)->size_t{return s+Ntype(t);}));
		for (const auto t: T) for (const auto i: inds(t)) res.push_back(i);
		return res;
	}
	//! returns the valid indices for atomic positions
	inline std::vector<size_t> inds() const noexcept { return ll__::rg(0,1,N()); }
	/** checks whether an atomic type is valid.
	 * @param t atomic type
	 */
	inline bool validType(const aT t) const noexcept { return t<Nspecies(); }
	/** check whether an index is valid.
	 * @param i index
	 */
	inline bool validInd(const size_t i) const noexcept { return i<N(); }
	//! returns the fundamental types, i.e. types that are not composite
	inline aTv fundamentalTypes() const noexcept {
		if (!this->id().empty()) return types();
		aTv res = types();
		res.resize(std::distance(res.begin(), std::remove_if(res.begin(),res.end(),
			[this](const aT t)->bool{return isComposite(id(t));})));
		return res;
	}
	//! returns the composite types
	inline aTv compositeTypes() const noexcept {
		if (!this->id().empty()) return {};
		aTv res = types();
		res.resize(std::distance(res.begin(), std::remove_if(res.begin(),res.end(),
			[this](const aT t)->bool{return !isComposite(id(t));})));
		return res;
	}
	/** returns the 'equal' types, i.e. types equal upto an index on the id string
	 * @return a vector of a vector of atomic types. Within each vector of atomic
	 * types, the types are 'equal'
	 */
	std::vector<aTv> equalTypes() const noexcept;


	/** @name id information
	 */
	//! returns id strings const reference
	inline const idv& id() const noexcept { return id_; }
	/** returns id string for an atomic type
	 * @param t atomic type
	 */
	inline const std::string& id(const aT t) const noexcept {
		return t<id().size() ? id()[t]: e_;
	}
	/** returns a range of id strings for a range of atomic types
	 * @param T range of atomic types
	 */
	inline idv id(const aTv& T) const noexcept {
		idv res(T.size());
		std::transform(T.cbegin(),T.cend(),res.begin(),
			[this](const auto t)->std::string{return id(t);});
		return res;
	}
	/** returns a 'soft' stripped id string.
	 * Stripping mean removing indices from id strings. Soft stripping
	 * will remove only the last index. For fundamental ids the difference
	 * is null, but for composite ids it matters.
	 * E.g. (Mo-001:S-002)-003 -> (Mo-001:S-002)
	 * @param id to soft strip
	 */
	inline static std::string softstripId(const std::string& id) noexcept {
		const size_t i = id.find_last_of(")-");
		return i!=std::string::npos && id[i]=='-' ? id.substr(0,i): id;
	}
	/** returns the index of an id string, i.e. what is removed when
	 * 'soft' stripping an id string.
	 * @param id id string to extract index from.
	 */
	inline static size_t getIndex(const std::string& id) noexcept {
		const size_t i = id.find_last_of(")-");
		return i!=std::string::npos && id[i]=='-' ? std::stod(id.substr(i+1)): 0;
	}
	/** returns fully stripped id stings.
	 * This will strip all indices, even from composite id strings.
	 * E.g. (Mo-001:S-002)-003 -> (Mo:S)
	 * @param id id string to strip
	 */
	inline static std::string stripId(const std::string& id) noexcept {
		const std::regex rgx("-[0-9]+");
		return std::accumulate(
			std::sregex_token_iterator(id.cbegin(),id.cend(),rgx,-1),
			std::sregex_token_iterator(),std::string(),
			[](const std::string& a, const std::string& b){return a+b;});
	}
	/** returns whether an id string is composite.
	 * @param id id string to check
	 */
	inline static bool isComposite(const std::string& id) noexcept {
		return id.find_first_of("():_")!=std::string::npos;
	}
	//! returns the fundamental ids
	idv fundamentalIds() const noexcept;
	//! returns the composite ids
	idv compositeIds() const;
	

	/** @name mass queries
	 */
	/** returns the mass of the fundamental type of an id string.
	 * Composite ids will return 0 mass.
	 * @param id id string to check
	 */
	static inline double mass(const std::string& id) noexcept {
		const auto itr = ll__::atom::mass.find(stripId(id));
		return (itr==ll__::atom::mass.cend()) ? 0.0: itr->second;
	}
	/** returns the masses of the fundamental types of a range of id strings
	 * Composite ids will return 0 mass.
	 * @param id range of id strings to check
	 */
	static inline std::vector<double> mass(const idv& id) noexcept {
		std::vector<double> res(id.size());
		std::transform(id.cbegin(),id.cend(),res.begin(),
			[](const std::string& s)->double{return mass(s);});
		return res;
	}
	//! returns the masses for all id strings
	inline std::vector<double> mass() const noexcept { return mass(id()); }


	/** @name type information
	 */
	/** returns the atomic type for an index
	 * @param i index
	 */
	inline aT type(const size_t i) const noexcept {
		const auto j = std::lower_bound(T_.cbegin()+1,T_.cend(),i+1);
		return j==T_.cend() ? NPOS__: std::distance(T_.cbegin()+1,j);
	}
	/** returns a range of atomic types for a range of atomic indices
	 * @param I range of indices
	 */
	inline aTv type(std::vector<size_t> I) const noexcept {
		std::transform(I.cbegin(),I.cend(),I.begin(),
			[this](const size_t i) -> aT {return type(i);});
		return I;
	}
	//! returns the atomic types
	inline aTv type() const noexcept { return type(inds()); }
	/** returns the atomic type for an id string
	 * @param s id string
	 */
	inline aT type(const std::string& s) const noexcept {
		const auto itr = std::lower_bound(id().cbegin(),id().cend(),s);
		return (itr==id().cend() || *itr!=s) ? NPOS__: std::distance(id().cbegin(),itr);
	}
	/** returns a range of atomic types for a range of id strings
	 * @param id range of id strings
	 */
	inline aTv type(const idv& id) const noexcept {
		std::vector<size_t> res(id.size());
		std::transform(id.cbegin(),id.cend(),res.begin(),
			[this](const auto& i)->size_t{return type(i);});
		return res;
	}
	/** returns a range of atomic types for whom the stripped id string is equal to s
	 * @param s id string
	 */
	inline aTv strippedType(const std::string& s) const noexcept {
		aTv res;
		for (size_t i=0; i!=id().size(); ++i)
			if (stripId(id()[i])==s) res.push_back(i);
		return res;
	}
	/** returns the number of positions of an atomic type
	 * @param t atomic type
	 */
	inline aC Ntype(const aT t) const noexcept { return validType(t) ? T_[t+1]-T_[t]: 0; }
	/** returns a range of atomic counts for a range of atomic types
	 * @param ts range of atomic types
	 * @return range of atomic counts
	 */
	inline aCv Ntype(aTv ts) const noexcept {
		std::transform(ts.cbegin(),ts.cend(),ts.begin(),
			[this](const aT t) -> aC {return Ntype(t);});
		return ts;
	}
	/** returns the atomic counts for all atomic types
	 * @return range of atomic counts
	 */
	inline aCv Ntype() const noexcept { return Ntype(types()); }
	/** returns the least frequent type of a range of atomic types, i.e. of the lowest atomic count
	 * @param ts range of atomic types
	 * @return atomic type
	 */
	inline aT leastFreqType(const aTv& ts) const noexcept {
		const auto N = Ntype(ts);
		return std::distance(N.cbegin(),std::min_element(N.cbegin(),N.cend()));
	}
	//! returns the least frequent atomic type of all atomic types
	inline aT leastFreqType() const noexcept { return leastFreqType(types()); }
	/** returns the most frequent type of a range of atomic types, i.e. of the highest atomic count
	 * @param ts range of atomic types
	 * @return atomic type
	 */
	inline aT mostFreqType(const aTv& ts) const noexcept {
		const auto N = Ntype(ts);
		return std::distance(N.cbegin(),std::max_element(N.cbegin(),N.cend()));
	}
	//! returns the most frequent atomic type of all atomic types
	inline aT mostFreqType() const noexcept { return mostFreqType(types()); }
	/** returns a range indices for an atomic type
	 * @param t atomic type
	 * @return vector of indices
	 */
	inline std::vector<size_t> ind(const aT t) const noexcept {
		return validType(t) ? ll__::rg(T_[t],1,T_[t+1]-T_[t]): std::vector<size_t>();
	}


	/** @name modification
	 */
	/** sets the id string vector to a new one
	 * @param Nid new id vector
	 */
	inline ll_cell& setId(idv Nid) noexcept {
		assert(Nid.size()==Nspecies());
		id_ = std::move(Nid);
		indexId_();
		return *this;
	}
	/** stresses the unit cell, i.e. applies a stress tensor to B
	 * @param S stress tensor, size(S) = size(B)
	 */
	inline ll_cell& stress(const fMat& S) noexcept {
		assert(S.square() && S.M()==dim());
		assert(lm__::det(S)>lm__::mtol());
		B_ = S.prod(B_);
		return *this;
	}
	/** swap two dimensions (columns) of B.
	 * This will cause a resorting of the resulting atomic positions.
	 * @param n1 index for dimension 1
	 * @param n2 index for dimension 2
	 */
	ll_cell& swapDim(const size_t n1, const size_t n2) noexcept;
	/** invert a dimension (a column) of B.
	 * This will cause a resorting of the resulting atomic positions.
	 * @param n dimension index
	 */
	ll_cell& invDim(const size_t n) noexcept;
	/** orient the cell to specific sign.
	 * This is done by inverting a dimension (1D) or swapping dimensions (>=2D).
	 * This will cause a resorting of the resulting atomic positions.
	 * @param sgn sign to be enforced
	 * @param r restriction vector, to determine 'frozen' dimensions
	 */
	ll_cell& orient(const double sgn, const rv& r) noexcept;
	/** permute the dimensions (columns) of B
	 * This will cause a resorting of the resulting atomic positions.
	 * @param P permutation matrix, size(P) = size(B)
	 */
	ll_cell& permute(const fMat& P) noexcept;
	/** rotate unit cell.
	 * This will NOT cause a resorting of the resulting atomic positions.
	 * @param R rotation matrix, size(R) = size(B)
	 */
	ll_cell& rotate(const fMat& R) noexcept;
	/** scale the unit cell, i.e. B = f*B
	 * This will NOT cause a resorting of the resulting atomic positions.
	 * @param f scaling factor
	 */
	ll_cell& scale(const double f) noexcept;
	/** scale the unit cell, i.e. b_i = f_i*b_i
	 * This will NOT cause a resorting of the resulting atomic positions.
	 * @param f scaling array, one factor per dimension
	 */
	ll_cell& scale(const fArray& f) noexcept;
	/** change the underlying basis, i.e periodicity.
	 * The new basis must be valid for the same lattice but not of the same symmetry group.
	 * @param NB new basis
	 */
	ll_cell& changeBasis(const fMat& NB) noexcept;
	/** change the underlying basis, i.e periodicity by an expansion.
	 * A valid expansion is a matrix of integers and |det|>=1.
	 * This is equivalent to change basis with NB = B*EC
	 * @param EC expansion matrix
	 */
	inline ll_cell& expand(const fMat& EC) noexcept {
		assert(EC==lm__::round(EC));
		return changeBasis(B().prod(EC));
	}
	//! attempts to make the unit cell primitive, i.e. reduce the volume to the minimum.
	inline ll_cell& makePrimitive() noexcept {
		return makePrimitive_(lm__::zeros<fMat>(dim()),rv(dim(),false),
				[](const fMat& vec)->bool{return true;});
	}
	/** attempts to make the unit cell primitive, i.e. reduce the volume to the minimum.
	 * @param r restriction vector to indicate 'frozen' columns of B
	 */
	ll_cell& makePrimitive(const rv& r) noexcept;
	/** attempts to make the unit cell primitive, i.e. reduce the volume to the minimum in a subspace.
	 * This will completely freeze one or more dimenions and run the algorithm on the remaining
	 * unit cell of lower dimenions.
	 * @param r restriction vector to indicate 'frozen' dimenions
	 */
	ll_cell& makePrimitiveInSubspace(const rv& r) noexcept;
	/** shift the atomic positions by a constant vector.
	 * This will NOT cause a resorting of the resulting atomic positions.
	 */
	ll_cell& shift(const fCol& sh) noexcept;
	/** shift the atomic positions by a constant vector.
	 * This will NOT cause a resorting of the resulting atomic positions.
	 */
	inline ll_cell& shift(const fMat& sh) noexcept { assert(sh.L()==dim()); return shift(sh.cAt(0)); }
	/** attempts to automatically shift atomic positions.
	 * This will try to maximize the 'leftness' of the unit cell, i.e. it will attempt to place
	 * the largest gap between positions at the right of the interval [0,1)
	 * @param r restriction vector ro indicate 'frozen' dimensions
	 */
	ll_cell& autoShift(const rv& r) noexcept;
	/** attempts to automatically shift atomic positions.
	 * This will try to maximize the 'leftness' of the unit cell, i.e. it will attempt to place
	 * the largest gap between positions at the right of the interval [0,1)
	 */
	inline ll_cell& autoShift() noexcept { return autoShift(rv(dim(),false)); }
	/** this will 'diversify' a range of atomic types. Henceforth they will be considered
	 * distinguished and the corresponding id strings will be reindexed. Indexed ids always
	 * start from -001
	 * @param ts range of atomic types
	 */
	ll_cell& diversify(const aTv& ts) noexcept;
	/** this will 'diversify' all atomic types. Henceforth they will be considered
	 * distinguished and the corresponding id strings will be reindexed. Indexed ids always
	 * start from -001
	 */
	inline ll_cell& diversify() noexcept { return diversify(types()); }
	/** this will 'collectivize' a range of atomic types. Henceforth they will be considered
	 * of the same type. The new id string is that of the first in the range
	 * @param ts range of atomic types
	 */
	ll_cell& collectivize(const aTv& ts) noexcept;
	/** this will 'collectivize' all atomic types. Henceforth they will be considered
	 * of the same type. The new id string is that of the first in the range
	 */
	inline ll_cell& collectivize() noexcept { return collectivize(types()); }
	/** this will attempt to merge two unit cells into one. Same id strings (including index)
	 * will be considered the same type.
	 * @param inp the cell wot merge this one with
	 */
	ll_cell& merge(const ll_cell& inp) noexcept;
	/** this will simplify composite ids since they can get very long in some cases.
	 * They will be replaced with '(CT)' and indexed afterwards.
	 */
	inline ll_cell& simplifyCompositeIds() noexcept {
		std::replace_if(id_.begin(),id_.end(),
			[](const std::string& s)->bool{return isComposite(s);},"(CT)");
		indexId_(); return *this;
	}


	/** @name conversion
	 */
	//! returns a copy of this
	inline ll_cell copy() const noexcept { return *this; }
	//! returns the reciprocal basis, i.e. 2*pi*B^-T
	inline fMat getRB() const noexcept {
		assert(!empty());
		return 2.0*M_PI*lm__::inv(B()).T();
	}
	/** returns atomic positions of an atomic type
	 * @param t atomic type
	 */
	inline fMat getAp(const aT t) const noexcept {
		return validType(t) ? Ap().get(0,T_[t],dim(),Ntype(t)): fMat(dim(),0);
	}
	/** returns atomic positions for a range of atomic types
	 * @param ts range of atomic types
	 */
	fMat getAp(const aTv& ts) const noexcept;
	//! returns atomic positions in cartesian coordinates
	inline fMat getcAp() const { return B().prod(Ap()); }
	/** returns atomic positions of an atomic type in cartesian coordinates
	 * @param t atomic type
	 */
	inline fMat getcAp(const aT t) const { return B().prod(getAp(t)); }
	/** returns atomic positions for a range of atomic types in cartesian coordinates
	 * @param ts range of atomic types
	 */
	inline fMat getcAp(const aTv& ts) const noexcept { return B().prod(getAp(ts)); }
	/** returns a subcell with only a subset of atomic types included.
	 * @param ts range of atomic types
	 */
	ll_cell getSubCell(const aTv& ts) const noexcept;
	/** returns bonds between atomic positions.
	 * Bonds across boundaries to neighbor cells can be included or excluded with an
	 * appropriate next neighbor matrix.
	 * @param NN next neighbor matrix
	 * @param keep lambda to determine which bonds to include
	 */
	ll_bonds<i_i_R> getBonds(const fMat& NN,
		const std::function<bool(const fMat&,const i_i&)>& keep) const noexcept;
	/** returns bonds between atomic positions of a certain type with a certain type.
	 * Bonds across boundaries to neighbor cells can be included or excluded with an
	 * appropriate next neighbor matrix.
	 * @param T1 range of atomic types (start points)
	 * @param T2 range of atomic types (ending points)
	 * @param f tolerance factor, f<0: upper limit is |f| times shortest bond length
	 *                            f>0: upper limit is f
	 * @param NN next neighbor matrix
	 */
	ll_bonds<i_i_R> getBonds(const aTv& T1, const aTv& T2,
				 const double f, const fMat& NN) const noexcept;
	/** returns bonds between atomic positions.
	 * Bonds across boundaries to neighbor cells can be included or excluded with an
	 * appropriate next neighbor matrix.
	 * @param f tolerance factor, f<0: upper limit is |f| times shortest bond length
	 *                            f>0: upper limit is f
	 * @param NN next neighbor matrix
	 */
	ll_bonds<i_i_R> getBonds(const double f, const fMat& NN) const noexcept;


	/** @name comparison
	 */
	/** returns a permutation matrix if the other cell is the same upto a permutation
	 * of the unit cell. Returns empty otherwise.
	 * @param inp other cell to compare to
	 * @return permutation matrix of size(B) or empty matrix
	 */
	fMat getPmat(const ll_cell& inp) const noexcept;
	/** returns a shift vector if the other cell is the same upto a permutation or a shift.
	 * Returns empty otherwise.
	 * @param inp other cell to compare to
	 * @param P permutation matrix
	 * @return vector of size 1xDIM, or empty
	 */
	fMat getAvec(const ll_cell& inp, const fMat& P) const noexcept;
	/** returns a shift vector if the other cell is the same upto a shift.
	 * Returns empty otherwise.
	 * @param inp other cell to compare to
	 * @return vector of size 1xDIM, or empty
	 */
	inline fMat getAvec(const ll_cell& inp) const noexcept { 
		return getAvec(inp,lm__::eye<fMat>(dim()));
	}
	/** check whether another cell describes the same lattice as this one
	 * @param inp other cell to compare to
	 */
	bool sameLattice(const ll_cell& inp) const noexcept;
	//! compares all containers directly
	inline bool operator==(const ll_cell& rhs) const noexcept {
		return B()==rhs.B() && Ap()==rhs.Ap() && T_==rhs.T_ && id()==rhs.id();
	}
	//! compares all containers directly
	inline bool operator!=(const ll_cell& rhs) const noexcept { return !(*this==rhs); }


	/** @name printing
	 */
	/** print the unit cell to a string in VASP POSCAR format
	 * @param direct print in direct coordinates
	 * @param prec floating point precision
	 */
	std::string print(const bool direct=true, const size_t prec=PPREC__) const noexcept;
	//! print the unit cell directly to a file
	void printToFile(const std::string& fileName) const;
	//! print to file overload
	void writeToFile(const std::string& fileName) const { printToFile(fileName); }


	/** @name friends
	 */
	friend class test_cell_assign;			//!< friend test class
	friend class test_cell_information;		//!< friend test class
	friend class test_cell_type_information;	//!< friend test class
	friend class test_cell_modification;		//!< friend test class
	friend class test_cell_conversion;		//!< friend test class


private:
	/** @name helpers
	 */
	//! general make primitive function
	ll_cell& makePrimitive_(fMat C, const rv& r, const std::function<bool(const fMat&)>& p) noexcept;
	

protected:
	/** @name members variables
	 */
	fMat B_;			//!< basis 
	std::vector<size_t> T_;		//!< atomic type vector
	idv id_;			//!< id strings vector
	static const std::string e_;	//!< empty string


protected:
	/** @name non const iterators
	 */
	/** iterator to the beginning of the range of a type
	 * @param t atomic type
	 */
	inline auto cBegin(const aT t) noexcept {
		assert(validType(t)); return mat_.cBegin()+T_[t];
	}
	/** iterator past the end of the range of a type
	 * @param t atomic type
	 */
	inline auto cEnd(const aT t) noexcept {
		assert(validType(t)); return mat_.cBegin()+T_[t+1];
	}
	/** sort Ap_ inside each type and check
	 * @param ts range of atomic types
	 */
	inline void sortAp_(const aTv& ts) noexcept {
		for (const auto t: ts) std::sort(cBegin(t),cEnd(t));
	}
	/** check if Ap_ is sorted
	 * @param ts range of atomic types
	 */
	inline bool Apsorted_(const aTv& ts) const noexcept {
		for (const auto t: ts)
			if (!std::is_sorted(ccBegin(t),ccEnd(t)))
				return false;
		return true;
	}
	//! index non unique ids
	void indexId_() noexcept;
	/** appends indices to an id string
	 * @param id id string
	 * @param i index to append
	 */
	inline static std::string appendIndex_(const std::string& id, const size_t i) noexcept {
			std::stringstream sstr;
			sstr << id << '-' << std::setfill('0')
					  << std::setw(NLZ__+1) << i;
			return sstr.str();
	}
	//! detect file type
	size_t detectFileType(const std::string& fileName) const;
};


//! streaming operator
std::ostream& operator<<(std::ostream& os, const ll_cell& inp) noexcept;

#endif // _LL_CELL_

/** @}
 */
