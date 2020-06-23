// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_HBONDSS_
#define _LL_HBONDSS_

#include "ll_hbonds.h"


/** Input for multiple wbh adaptations.
 * As a child of aux_parser, it redefines its own versions of
 * parseKey_ and printHelp_.
 */
struct ll_hbondss_input: public virtual aux_parser {
public:
	/** @name filenames
	 */
	std::vector<std::string> wbh = {WBH__};		//!< filenames for each input ll_bonds
	std::vector<std::string> wbh_out = {};		//!< filenames for each output ll_bonds

	/** @name parameters
	 */
	double sptol = 1e-1;				//!< spacial tolerance when querying for bonds
	double perimeter_radius = .0;			//!< radius factor to use when detecting bonds across
							//!< the perimeter between wbhs
	double signs_tol = 1e-3;			//!< tolerance level for sign equalizing

public:
	/** parseKey_ redefinition
	 */
	struct parseKey_: public virtual aux_parser::parseKey_ {
	/** constructor defining keys
	 */
	inline explicit parseKey_(ll_hbondss_input& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator i, std::sregex_token_iterator e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			// filenames
			case "wbh"_h: PARSE__(p.wbh); return;
			case "wbh_out"_h: PARSE__(p.wbh_out); return;
			
			// parameter
			case "sptol"_h: PARSE__(p.sptol); return;
			case "perimeter_radius"_h: PARSE__(p.perimeter_radius); return;
			case "signs_tol"_h: PARSE__(p.signs_tol); return;
		}
	}
	};

public:
	/** printHelp_ redefinition
	 */
	struct printHelp_: public virtual aux_parser::printHelp_ {
	/** constructor defining help messages
	 */
	inline explicit printHelp_(const ll_hbondss_input& p, std::ostream& os):
			aux_parser::printHelp_(p,os) {
		printHelpTuple_(os,std::make_tuple(
		
		TOPIC__("HBONDSS: FILENAMES")
		,HELP__("wbh",p.wbh,
			"names of the input wannier bonds hamiltonians")
		,HELP__("wbh_out",p.wbh_out,
			"names of the output wannier bonds hamiltonians\n"
			"after equalizing shifts")

		,TOPIC__("HBONDSS: PARAMETERS")
		,HELP__("sptol",p.sptol,
			"spacial tolerance to use when looking up interactions along bonds\n"
			"in real space in Angstroms")
		,HELP__("perimeter_radius",p.perimeter_radius,
			"cutoff radius to use when searching for bonds across the perimeter\n"
			">0.0: factor for the bond length cutoff of the intrinsic radius\n"
			"<0.0: direct cutoff in Angstroms\n"
			"=0.0: turn off adaptation")
		,HELP__("signs_tol",p.signs_tol,
			"tolerance level for as to which elements in the test hamiltonians\n"
			"should be discared when equalizing wannier signs")
		));
	}
	};
};



/** multiple ll_hbonds.
 * This class is designed with a broadly equivalent interface as 'll_hbonds' and as such
 * either can be used in many contexts. The generalized bond indices used here hold the
 * wbh index in the first byte and the regular index in the bytes above. The naked numbers
 * therefore appear to not make much sense. In any case, this class can be used without
 * adding significant overhead compared to ll_hbonds or even as a complete replacement
 * for the case of just one wbh. The approach here is required when generating Hamiltonian
 * matrices from multiple Wannierizations.
 */
class ll_hbondss final: public ll__::vec_cb<std::vector<ll_hbonds>,ll_hbonds> {
public:
	/** @name types
	 */
	typedef lm__::fMat fMat;		//!< real matrix
	typedef lm__::cMat cMat;		//!< complex matrix
	typedef lm__::fArray fArray;		//!< real array
	typedef ll__::idv idv;			//!< range of id strings
	typedef ll__::i_i i_i;			//!< index pair {i1,i2}


protected:
	/** @name helpers for handling types
	 */
	//! extracts wbh index from generalized index
	static inline size_t j_(const size_t i) noexcept {
		return i&size_t(255);
	}
	//! extracts normal index from generalized index
	static inline size_t i_(const size_t i) noexcept {
		return i>>8;
	}
	//! composes generalized index from wbh index and normal index
	static inline size_t compInd_(const size_t j, const size_t i) noexcept {
		assert(j<256);
		return (i<<8) + j;
	}
	//! decomposes generalized index into wbh index and normal index
	static inline auto decompInd_(const size_t i) noexcept {
		struct res_st { size_t j,i; };
		return res_st {j_(i),i_(i)};
	}
	//! decomposes generalized index starting from an id string
	static inline auto decompId_(const std::string& s) {
		struct res_st { size_t j; std::string id; };
		
		const size_t pos = s.find(':');
		if (pos == std::string::npos)
			return res_st{0,s};
		
		const std::string js = s.substr(0,pos);
		const std::string id = s.substr(pos+1);

		try {
			return res_st{std::stoul(js),std::move(id)};
		} catch(const std::exception& e) {
			return res_st{0,s};
		}
	}


public:
	/** @name constructors
	 */
	using ll__::vec_cb<std::vector<ll_hbonds>,ll_hbonds>::vec_cb;
	/** constructor from user inpput
	 * @param inp user input struct
	 * @param os stream to print into
	 */
	ll_hbondss(const ll_hbondss_input& inp, std::ostream& os);


	/** @name information
	 */
	//! reference to empty interaction matrix
	inline const cMat& eH() const noexcept { return this->front().eH(); }
	//! dimension of space
	inline size_t dim() const noexcept { return this->front().dim(); }
	//! bond radius encompassing bonds in all wbh
	inline double radius() const noexcept {
		double res = 0.0;
		for (const auto& w: *this) {
			const double r = w.radius();
			if (r>res) res=r;
		}
		return res;
	}
	/** returns a restriction matrix for space. A spacial dimension is considered restricted if
	 * none of the R vectors have entries greater other than 0 in it.
	 */
	inline fMat rmat() const noexcept {

		auto res = lm__::ones<fMat>(dim(),1);
		for (const auto& r: *this)
			res &= r.rmat();
		return res;
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
	//! returns the number of orbitals for an index
	inline size_t Norb(const size_t i) const noexcept {
		return j_(i)<this->size() ? (*this)[j_(i)].Norb(i_(i)): 0;
	}
	//! returns the number of orbitals for a range of indices
	inline std::vector<size_t> Norb(std::vector<size_t> inds) const noexcept {
		std::transform(inds.cbegin(),inds.cend(),inds.begin(),[this](const size_t ind)->size_t{
			return this->Norb(ind);
		});
		return inds;
	}
	//! returns the indices for which there exist bonds
	inline std::vector<size_t> inds() const noexcept {

		std::vector<size_t> res;
		res.reserve(std::accumulate(this->cbegin(),this->cend(),size_t(0),
				[](const size_t s, const auto& w) -> size_t
				{ return s + w.cell().N(); }));
		
		for (size_t j=0; j!=this->size(); ++j)
			for (const size_t i: (*this)[j].inds())
				res.push_back(compInd_(j,i));

		return res;
	}
	/** returns the index for an id string
	 * @param s 'atomic' id string, e.g. 'Mo'
	 */
	inline size_t ind(const std::string& s) const {
		const auto I = decompId_(s);
		if (I.j>=this->size()) return NPOS__;
		const size_t res = (*this)[I.j].ind(I.id);
		return res==NPOS__ ? NPOS__: compInd_(I.j,res);
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
	/** returns the id string for a wbh index and a normal bond index
	 * @param j wbh index
	 * @param i normal bond index
	 */
	inline std::string id(const size_t j, const size_t i) const noexcept {
		return j<this->size() ? std::to_string(j)+':'+(*this)[j].id(i): "";
	}
	/** returns the id string for a generalized bond index
	 * @param i generalized bond index
	 */
	inline std::string id(const size_t i) const noexcept {
		return id(j_(i),i_(i));
	}
	/** returns a range of id strings for a range of generalized bond indices
	 * @param I range of generalized bond indices
	 */
	inline idv id(const std::vector<size_t>& I) const noexcept {
		idv res(I.size());
		std::transform(I.cbegin(),I.cend(),res.begin(),
			[this](const size_t i)->std::string{ return id(i); });
		return res;
	}
	//! returns the tolerance used for querying bonds
	inline std::vector<double> queryTol() const noexcept {
		std::vector<double> res; res.reserve(this->size());
		for (const auto& w: *this) res.push_back(w.queryTol());
		return res;
	}
	/** returns the neighborhood for a generalized bond index. i.e. all the bonds
	 * 			starting from this index
	 * @param ind starting index
	 * @param Rcut cutoff radius for the bonds
	 * @param sptol spacial tolerance used when sorting the bonds
	 */
	ll__::b_I neighborhood(const size_t ind, const double Rcut,
			const double sptol) const noexcept {
		assert(j_(ind)<this->size());
		
		const size_t j = j_(ind);
		auto res = (*this)[j_(ind)].neighborhood(i_(ind),Rcut,sptol);
		for (auto& ii: res)
			ii = {compInd_(j,ii.i1()),compInd_(j,ii.i2())};

		return res;
	}


	/** @name information from structure
	 */
	//! index struct for matches across domain borders
	union pbrt {
		struct {
			size_t j1;		//!< section1 index
			size_t i1;		//!< index inside section 1
			size_t j2;		//!< section2 index
			size_t i2;		//!< index inside section 2
			size_t a2;		//!< actual target in section2
			size_t a1;		//!< actual target in section1
		};
		std::array<size_t,6> dat;	//!< access the struct as one array
		
		//! linear easy access operator
		inline size_t operator[](const size_t i) const noexcept {
			assert(i<6); return dat[i];
		}
		//! easy access operator
		inline size_t operator()(const size_t i, const size_t j) const noexcept {
			assert(i<3); assert(j<2); return dat[i*2+j];
		}
	};
	//! range of index matches and corresponding bonds
	struct pb {
		std::vector<pbrt> RT;	//!< roots and targets
		fMat bnds12;		//!< corresponding bonds
	};
	/** function determining the connections at the perimeter between Wannier domains.
	 * @param Ap atomic positions for the structure under investigation
	 * @param I generalized bond indices corresponding to Ap
	 * @param inp user input struct
	 * @param os stream to print into
	 * @return bonds going across perimeters, the roots, fake targets and actual target indices
	 */
	pb getPerimeterConnections(const fMat& Ap, const std::vector<size_t>& I,
			const ll_hbondss_input& inp, std::ostream& os) const;
	/** returns a vector of index pairs indicating possible substitute atoms.
	 * This can be used for 'inexact matching' where some atoms are allowed to
	 * substitute others.
	 * @param Ap a list of atomic positions.
	 * @param I indices corresponding to Ap
	 * @param f factor for the cutoff, f>0: Rcut=f*radius(), f<0: Rcut=-f
	 * @return sorted vector of index pairs {i1,i2} where an atom at index i1 may
	 * 	   be substituted with one of index i2 and vice versa.
	 */
	std::vector<std::vector<i_i>> getSubstitutes(const fMat& Ap,
			const std::vector<size_t>& I,const double f=.5) const noexcept;
	
	
	/** @name modification
	 */
	/** function to equalize relative energy shifts between Wannier domains.
	 * @param B bonds across perimiters as produced by getPerimeterConnections
	 * @param inp user input script
	 * @param os stream to print into
	 */
	void equalizeShifts(const pb& B, const ll_hbondss_input& inp, std::ostream& os);
	/** function to equalize Wannier signs. Here signs of Wannier signs are adapted
	 * to match for bond indices matched by substitution across perimiters of Wannier
	 * domains.
	 * @param B bonds across perimiters as produced by getPerimeterConnections
	 * @param SI range of substitute index pairs, as produced by getSubstitutes
	 * @param inp user input script
	 * @param os stream to print into
	 */
	void equalizeSigns(const pb& B, const std::vector<std::vector<i_i>>& SI,
			const ll_hbondss_input& inp, std::ostream& os);
	/** function to flip Wannier signs in the interaction matrices. This must be done on the
	 * off-diagonal element of interaction matrices.
	 * @param ind bond index
	 * @param F vector of bools as to which signs to flip
	 */
	inline void flipWannierSigns(const size_t ind, const std::vector<bool>& F) noexcept {
		assert(j_(ind)<this->size());
		this->vec_[j_(ind)].flipWannierSigns(i_(ind),F);
	}
	/** function to 'fix' the diagonals of interactions matrices. This adds a small value
	 * to each entry on the diagonals (because OMEN likes this).
	 * @param tol small parameter
	 */
	inline void fixDiagonals(const double tol) noexcept {
		std::for_each(this->begin(),this->end(),[tol](auto& W)
			{ W.fixDiagonals(tol); });
	}


	/** @name searching
	 */
	//! set bond query tolerance in direct coordinates
	inline void setQueryTolDirect(const double tol) const noexcept {
		assert(tol>=.0);
		for (auto& w: *this) w.setQueryTolDirect(tol);
	}
	//! set bond query tolerance in cartesian coordinates
	inline void setQueryTolCartesian(const double tol) const noexcept {
		for (auto& w: *this) w.setQueryTolCartesian(tol);
	}
	/** function to query for an interaction. This is 'exact' matching of bonds.
	 * @param j an index pair {i1,i2}
	 * @param b a bond vector
	 * @return the corresponding interaction matrix
	 */
	inline const cMat& getInteraction(const i_i& j, const fArray& b) const noexcept {
		return j.i1()==NPOS__ || j.i2()==NPOS__ || j_(j.i1())!=j_(j.i2()) ? eH(): 
			(*this)[j_(j.i1())].getInteraction({i_(j.i1()),i_(j.i2())},b);
	}
	/** function to query for an interaction. This is 'inexact' matching of bonds
	 * @param j in index pair {i1,i2}
	 * @param b a bond vector
	 */
	ll_hbonds::am_ getApproximateInteraction(const i_i& j, const fArray& b) const noexcept;


	/** @name printing
	 */
	//! streaming operator
	friend inline std::ostream& operator<<(std::ostream& os, const ll_hbondss::pb& inp) noexcept {
		auto b = inp.bnds12.ccBegin();
		for (auto i=inp.RT.cbegin(),e=inp.RT.cend(); i!=e; ++i,++b)
			os << "\n{" << i->j1 << "," << i->i1 << "," << i->a2 << "} "
			   <<   "{" << i->j2 << "," << i->i2 << "," << i->a1 << "} "
			   << lm__::T(*b);
		return os;
	}


	/** @name friend test class
	 */
	//! friend test class
	friend class test_hbondss_all;
};

#endif // _LL_HBONDSS_

/** @}
 */
