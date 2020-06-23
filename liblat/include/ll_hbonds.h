// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_HBONDS_
#define _LL_HBONDS_

#include "ll_types.h"
#include "aux_parser.h"
#include "ll_bonds.h"


/** Input for the wannier matching parameters.
 * As a child of aux_parser, it redefines its own versions of
 * parseKey_ and printHelp_.
 */
struct ll_wmatching_input: public virtual aux_parser {
public:
	/** @name filenames
	 */
	std::string wout = WOUT__;			//!< wannier90.wout file
	std::string hrdat = HR__;			//!< wannier90_hr.dat file
	
	/** @name material generation
	 */
	std::vector<std::string> mode = {"atomic"};	//!< mode string
	double tol = 1e-4;				//!< tolerance as to when hamiltonian data should be discarded
	std::vector<size_t> wbl = {};			//!< wannier black list
	std::vector<size_t> abl = {};			//!< atomic black list
	double bond_factor = -1.2;			//!< bond factor
	size_t rm_unmatched = 1;			//!< remove centers with no wannier centers matched to it
	bool simplify_composite = false;		//!< switch to simplify composite ids
	bool balance_matches = false;			//!< try to balance matches among same type
	size_t minpts = 2;				//!< minimum number of points per cluster
	double wradius = 1.5;				//!< search radius for clustering


public:
	/** parseKey_ redefinition
	 */
	struct parseKey_: public virtual aux_parser::parseKey_ {
	/** constructor defining keys
	 */
	inline explicit parseKey_(ll_wmatching_input& p, const uint32_t key,
		std::ifstream& file, size_t& lcnt,
		std::sregex_token_iterator i, std::sregex_token_iterator e):
			aux_parser::parseKey_(p,key,file,lcnt,i,e) {
		using namespace aux;
		
		switch(key) {
			// filenames
			case "wout"_h: PARSE__(p.wout); return;
			case "hrdat"_h: PARSE__(p.hrdat); return;
			
			// material generation	
			case "mode"_h: PARSE__(p.mode); return;
			case "tol"_h: PARSE__(p.tol); return;
			case "wbl"_h: PARSE__(p.wbl); return;
			case "abl"_h: PARSE__(p.abl); return;
			case "bond_factor"_h: PARSE__(p.bond_factor); return;
			case "rm_unmatched"_h: PARSE__(p.rm_unmatched); return;
			case "simplify_composite"_h: PARSE__(p.simplify_composite); return;
			case "balance_matches"_h: PARSE__(p.balance_matches); return;
			case "minpts"_h: PARSE__(p.minpts); return;
			case "wradius"_h: PARSE__(p.wradius); return;
		}
	}
	};

public:
	/** printHelp_ redefinition
	 */
	struct printHelp_: public virtual aux_parser::printHelp_ {
	/** constructor defining help messages
	 */
	inline explicit printHelp_(const ll_wmatching_input& p, std::ostream& os):
			aux_parser::printHelp_(p,os) {
		printHelpTuple_(os,std::make_tuple(
		
		TOPIC__("WMATCHING: FILENAMES")
		,HELP__("wout",p.wout,
			"name of the wannier90 *.wout file")
		,HELP__("hrdat",p.hrdat,
			"name of the wannier90 *.hrdat file")

		,TOPIC__("WMATCHING: MATERIAL GENERATION")
		,HELP__("mode",p.mode,
			"mode to use for wannier center matching to positions\n"
			" - atomic: use atomic positions as matching centers\n"
			" - atomic_bc: like atomic, but add bonds centers as well\n"
			" - cluster: use clustering algorithm to identify clusters of\n"
			"            wannier centers, then use centers of mass as positions\n"
			" - wannier: use the wannier centers themselves, i.e. no matching\n"
			" - 'filename': use positions as defined by a POSCAR file 'filename'")
		,HELP__("tol",p.tol,
			"cutoff tolerance for wannier90 hamiltonian data")
		,HELP__("wbl",p.wbl,NPOS__,
			"wannier functions blacklist, as to which hamiltonian entries and centers\n"
			"should be discarded")
		,HELP__("abl",p.abl,NPOS__,
			"atomic positions blacklist, as to which atomic centers should be discarded")
		,HELP__("bond_factor",p.bond_factor,
			"bonds length to use while searching for nearest neighbors\n"
			" - >0.0: direct cutoff in Angstroms\n"
			" - <1.0: adaptive cutoff, i.e. with bf=-1.1 it will keep\n"
			"           looking for bonds until it can't find anything\n"
			"           that is not larger than 10% more")
		,HELP__("rm_unmatched",p.rm_unmatched,
			"bitmask to keep or discard positions without wannier centers matched to it:\n"
			" - 1 (01): remove unmatched composite positions\n"
			" - 2 (10): remove unmatched fundamental positions\n"
			" - 3 (11): remove all unmatched positions")
		,HELP__("simplify_composite",p.simplify_composite,
			"switch to force simplify composite ids,\n"
			"set this to true if you get very long ids")
		,HELP__("balance_matches",p.balance_matches,
			"try to balance out the number of centers matched to atoms of the same type")
		,HELP__("minpts",p.minpts,
			"minimum number of points per cluster, cluster mode only")
		,HELP__("wradius",p.wradius,
			"search radius for clustering, cluster mode only")
		));
	}
	};
};



namespace ll__ {
	
	/** index pair including R vectors and interaction matrices
	 */
	class i_i_R_H final: public i_i, public mat_vec_b<fMat,fArray,std::vector<cMat>,cMat> {
	public:
		/** @name constructors
		 */
		//! default constructor
		i_i_R_H(const size_t d=DIM__) noexcept: i_i(), mat_vec_b(fMat(d,0),{}) {}
		/** constructor from indices
		 * @param i1 first index
		 * @param i2 second index
		 * @param d dimension of R, i.e. of space
		 */
		i_i_R_H(const size_t i1, const size_t i2, const size_t d) noexcept:
			i_i(i1,i2), mat_vec_b(fMat(d,0),{}) {}
		/** constructor from indices, R vectors and interaction matrices
		 * @param i1 first index
		 * @param i2 second index
		 * @param R R vectors
		 * @param H vector of interaction matrices
		 */
		i_i_R_H(const size_t i1, const size_t i2, fMat R, std::vector<cMat> H) noexcept:
			i_i(i1,i2), mat_vec_b(std::move(R),std::move(H)) {}
		

		/** @name information
		 */
		//! R vectors reference
		inline fMat& R() noexcept { return this->mat_; }
		//! R vectrors const reference
		inline const fMat& R() const noexcept { return this->mat_; }
		//! interaction matrrices reference
		inline std::vector<cMat>& H() noexcept { return this->vec_; }
		//! interaction matrices const reference
		inline const std::vector<cMat>& H() const noexcept { return this->vec_; }
		//! iterator to center interaction matrix (at R=0)
		inline auto center() noexcept {
			const auto itr = std::lower_bound(this->ccBegin(),this->ccEnd(),0.0);
			return (itr==this->ccEnd() || *itr!=0.0) ? this->end(): this->begin() + (size_t)itr;
		}
		//! const_iterator to center interaction matrix (at R=0)
		inline auto center() const noexcept {
			const auto itr = std::lower_bound(this->ccBegin(),this->ccEnd(),0.0);
			return (itr==this->ccEnd() || *itr!=0.0) ? this->cend(): this->cbegin() + (size_t)itr;
		}
	

		/** @name friends
		 */
		//! streaming operator
		friend inline std::ostream& operator<<(std::ostream& os, const i_i_R_H& inp) noexcept {
			return (os<<inp.i1()<<":"<<inp.i2()<<"\n"<<lm__::T(inp.R()));
		}
		//! conjugate, i.e. flip indices, invert R and transpose H
		friend inline i_i_R_H conj(const i_i_R_H& inp) noexcept {
			std::vector<cMat> cH;
			std::transform(inp.cbegin(),inp.cend(),cH.begin(),[](const cMat& inp)->cMat{
				return lm__::T(inp);
			});
			return i_i_R_H(inp.i2(),inp.i1(),-inp.R(),std::move(cH));
		}
	};
}



/** bonds with interactions class.
 * Each bond has an interaction matrix attached, which is directly constructed
 * from Wannier90 output. The contents of the underlying unit cell are generated
 * through a process we call Wannier matching. Wannier centers are grouped and
 * matched to atomic or possibly fictitious positions where they play the role
 * of 'atomic' orbitals. The number of orbitals for each 'atom' in the unit cell is 
 * therefore the number of Wannier functions matched to it. The physics of the
 * interactions between these atoms is described in terms of the Wannier functions
 * matched to them.
 */
class ll_hbonds final: public ll_bonds<ll__::i_i_R_H> {
public:
	// types
	typedef ll__::fMat fMat;	//!< real matrix
	typedef ll__::cMat cMat;	//!< complex matrix
	typedef ll__::fCol fCol;	//!< real column
	typedef ll__::wi wi;		//!< wannier matching
	typedef ll__::R_H<> R_H;	//!< R_H Hamiltonian container
	typedef ll__::idv idv;		//!< range of id strings

public:
	/** @name constructors
	 */
	//! default constructor
	inline ll_hbonds() noexcept: ll_bonds<ll__::i_i_R_H>(), Norb_() {}
	/** main constructor from Wannier90 input.
	 * @param inp the input script holding parameters for the matching process
	 * @param os stream to print into
	 */
	ll_hbonds(const ll_wmatching_input& inp, std::ostream& os);
	/** constructor from generated material and matched Wannier centers.
	 * @param cell unit cell holding 'atomic' positions corresponding to the Wannier center distribution
	 * @param Wp Wannier centers
	 * @param I Wannier matching indices, one vector per 'atomic' position in cell, one index per center in Wp
	 * @param hr Wannier90 Hamiltonian as generated by readHr
	 * @param transf generic lambda for manipulating interaction matrices
	 */
	ll_hbonds(ll_cell cell, const fMat& Wp, const wi& I, const R_H& hr,
	const std::function<cMat(cMat&&)>& transf=[](cMat&& inp)->cMat{return std::move(inp);}) noexcept;
	/** constructor from generated material and matched Wannier centers.
	 * @param cell unit cell holding 'atomic' positions corresponding to the Wannier center distribution
	 * @param Wp Wannier centers
	 * @param I Wannier matching indices, one vector per 'atomic' position in cell, one index per center in Wp
	 * @param hr Wannier90 Hamiltonian as generated by readHr
	 * @param keep lambda defining which interactions to keep
	 */
	inline ll_hbonds(ll_cell cell, const fMat& Wp, const wi& I, const R_H& hr,
			const std::function<bool(const cMat&)>& keep) noexcept:
		ll_hbonds(std::move(cell),Wp,I,hr,
			[keep](cMat&& inp)->cMat{return keep(inp) ? inp: cMat(); }) {}
	/** constructor from generated material and matched Wannier centers.
	 * @param cell unit cell holding 'atomic' positions corresponding to the Wannier center distribution
	 * @param Wp Wannier centers
	 * @param I Wannier matching indices, one vector per 'atomic' position in cell, one index per center in Wp
	 * @param hr Wannier90 Hamiltonian as generated by readHr
	 * @param tol tolerance level, interaction matrices with no entry above this will be discarded
	 */
	inline ll_hbonds(ll_cell cell, const fMat& Wp, const wi& I, const R_H& hr, const double tol) noexcept:
		ll_hbonds(std::move(cell),Wp,I,hr,[tol](cMat&& inp)->cMat{
			if (!std::any_of(inp.cbegin(),inp.cend(),[tol](const auto& inp)->bool
				{return ll__::cmph(inp,tol);})) return {};
			
			auto d=inp.cdbegin(), de=inp.cdend();
			for (auto i=inp.begin(),ie=inp.end(); i!=ie; ++i)
				if (d!=de && d==i) ++d;
				else { if (!ll__::cmph(*i,tol)) *i=0.0; }
			return std::move(inp);
		}) {}
	/** constructor from pregenerated data.
	 * @param cell unit cell holding 'atomic' positions
	 * @param Norb vector defining the number of Wannier centers or 'orbitals' matched to each position
	 * @param dat vector of index pairs, R vectors and interactions {i1,i2,R,H}
	 */
	inline ll_hbonds(ll_cell cell, std::vector<size_t> Norb, std::vector<ll__::i_i_R_H> dat) noexcept:
		ll_bonds<ll__::i_i_R_H>(std::move(cell),std::move(dat)), Norb_(std::move(Norb)) {}
	/** constructor from file
	 * @param fileName name of the file to read the wbh from
	 */
	ll_hbonds(const std::string& fileName);


public:
	/** @name data access
	 */
	//! reference to empty interaction matrix
	const cMat& eH() const noexcept { return eH_; }
	

	/** @name information
	 */
	//! total number of Hamiltonian entries in the class
	inline size_t NHentries() const noexcept {
		return std::accumulate(this->cbegin(),this->cend(),size_t(0),
		[](const size_t s, const auto& i)->size_t{
			return s+std::accumulate(i.cbegin(),i.cend(),size_t(0),
			[](const size_t s, const auto& i)->size_t{
				return s+i.size();
			});
		});
	}
	//! number of orbitals vector
	inline const std::vector<size_t>& Norb() const noexcept { return Norb_; }
	//! returns the number of orbitals for an index
	inline size_t Norb(const size_t i) const noexcept {
		return i<Norb_.size() ? Norb_[i]: 0;
	}
	//! returns the numbers of orbitals for a range of indices
	inline std::vector<size_t> Norb(std::vector<size_t> inds) const noexcept {
		std::transform(inds.cbegin(),inds.cend(),inds.begin(),[this](const size_t ind)->size_t{
			return this->Norb(ind);
		});
		return inds;
	}
	//! returns the number of orbitals for an id strings
	inline size_t Norb(const std::string& s) const noexcept { return Norb(this->ind(s)); }
	//! returns the number of orbitals for a range of id strings
	inline std::vector<size_t> Norb(const idv& id) const noexcept { return Norb(this->ind(id)); }
	/** returns a 'similarized cell' for which the id strings are replaced by strings of the
	 * form A{Norb}, e.g. 'A3' or 'A5'. All atoms with the same number of orbitals are thus
	 * treated as the same type.
	 */
	inline ll_cell symCell() const noexcept {
		if (empty()) return ll_cell(cell().B());

		idv id; id.reserve(cell().N());
		for (size_t i=0; i!=cell().N(); ++i)
			id.push_back("A"+std::to_string(Norb(i)));

		return ll_cell(cell().B(),cell().Ap(),std::move(id));
	}
	/** returns a vector of index pairs indicating possible substitute atoms.
	 * This can be used for 'inexact matching' where some atoms are allowed to
	 * substitute others.
	 * @param Ap a list of atomic positions.
	 * @param I indices corresponding to Ap
	 * @param f factor for the cutoff, f>0: Rcut=f*radius(), f<0: Rcut=-f
	 * @return sorted vector of index pairs {i1,i2} where an atom at index i1 may
	 * 	   be substituted with one of index i2 and vice versa.
	 */
	std::vector<i_i> getSubstitutes(const fMat& Ap, const std::vector<size_t>& I,
			const double f=.5) const noexcept;


	/** @name searching
	 */
	/** function to query for an interaction. This is 'exact' matching of bonds.
	 * @param j an index pair {i1,i2}
	 * @param b a bond vector
	 * @return the corresponding interaction matrix
	 */
	inline const cMat& getInteraction(const i_i& j, const fArray& b) const noexcept {
		const auto jtr = this->search(j);
		const auto itr = this->search(jtr,b);
		return itr!=this->e() ? jtr->H()[(size_t)itr]: eH();
	}
	//! helper struct for approximate interaction queries
	struct am_ {
		cMat H;		//!< The interaction matrix
		size_t pi;	//!< The substitute index for the bond in positive direction: 'b'
		size_t mi;	//!< The substitute index for the bond in negative direction: '-b'
		//! cast to cMat
		inline operator cMat() const noexcept { return H; }
	};
	/** function to query for an interaction. This is 'inexact' matching of bonds
	 * @param j in index pair {i1,i2}
	 * @param b a bond vector
	 */
	am_ getApproximateInteraction(const i_i& j, const fArray& b) const noexcept;


	/** @name modification
	 */
	//! function to filter bonds using a generic lambda on indices, start and end points, R and H
	size_t filter(const std::function<bool(
				const cMat& H, const fMat& B,
				const size_t i1, const size_t i2,
				const fMat& p1, const fMat& p2,
				const fCol& R)>& eval) noexcept;
	//! function to filter bonds using a lambda on just {i1,i2} and start and end points {p1,p2}
	inline size_t filter(const std::function<bool(
				const size_t i1, const size_t i2,
				const fMat& p1, const fMat& p2)>& eval) noexcept {
		return filter([&eval](const cMat& H, const fMat& B,
					const size_t i1, const size_t i2,
					const fMat& p1, const fMat& p2,
					const fCol& R) -> bool
				{ return eval(i1,i2,B.prod(p1),B.prod(p2+R)); });
	}
	/** function to shift the spectrum, i.e. shift bandstructures produces from bond interactions
	 * in energy. This shift must be applied to the self-interactions, i.e. bonds of length 0.
	 * @param sh energy shift including a possible imaginary component.
	 */
	inline void shiftSpectrum(const std::complex<double>& sh) noexcept {
		for (const auto i: this->inds()) {
			const auto itr = std::lower_bound(this->vec_.begin(),
							  this->vec_.end(),ll__::i_i{i,i});
			const auto jtr = itr->center();
			if (jtr!=itr->end())
				for (auto d=jtr->dbegin(),e=jtr->dend(); d!=e; ++d)
					*d += sh;
		}
	}
	/** function to flip Wannier signs in the interaction matrices. This must be done on the
	 * off-diagonal element of interaction matrices.
	 * @param ind bond index
	 * @param F vector of bools as to which signs to flip
	 */
	void flipWannierSigns(const size_t ind, const std::vector<bool>& F) noexcept;
	/** function to 'fix' the diagonals of interactions matrices. This adds a small value
	 * to each entry on the diagonals (because OMEN likes this).
	 * @param tol small parameter
	 */
	inline void fixDiagonals(const double tol) noexcept {
		for (auto& b: *this)
		for (auto& h: b)
			for (auto i=h.dbegin(),e=h.dend(); i!=e; ++i)
				if (std::abs(*i)<tol)
					*i = tol;
	}


	/** printing
	 */
	/** function to write the wbh to a file. The Fermi energy can be included since this way
	 * all the information to generate Hamiltonians/bandstructures is included in the same file.
	 * @param fileName name of the file
	 * @param Ef Fermi energy
	 */
	void writeToFile(const std::string& fileName, const double Ef) const;
	//! function to write the wbh to a file.
	inline void writeToFile(const std::string& fileName) const {
		writeToFile(fileName,std::nan(""));
	}
	/** print bonds to a stream.
	 * @param os stream to print into
	 * @param mode printing mode
	 */
	std::ostream& print(std::ostream& os, const size_t mode=0) const noexcept;
	//! stream operator
	friend inline std::ostream& operator<<(std::ostream& os, const ll_hbonds& inp) noexcept {
		return inp.print(os);
	}


protected:
	/** @name member variables
	 */
	static const cMat eH_;		//!< empty interaction matrix
	std::vector<size_t> Norb_;	//!< Norb vector
};

#endif // _LL_HBONDS_

/** @}
 */
