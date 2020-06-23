// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_TYPES_
#define _LL_TYPES_

#include "lm_defs.h"
#include "ll_defs.h"
#include "lm_tMat.h"
#include "aux_io.h"
#include "ll_compound.h"
#include <vector>
#include <map>


class test_io;
template<class MT> class ll_writerMEM;

namespace ll__ {
	
	using namespace lm__;
	using namespace aux;

	typedef size_t aT;					//!< atomic type
	typedef size_t aC;					//!< atomic count
	typedef size_t wT;					//!< wannier type
	typedef size_t wC;					//!< wannier count
	typedef std::complex<double> hel;			//!< hamiltonian element

	typedef std::vector<bool> rv;				//!< restriction vector
	typedef std::vector<aT> aTv;				//!< atomic type vector
	typedef std::vector<aC> aCv;				//!< atomic count vector
	typedef std::vector<std::string> idv;			//!< atomic id strings
	typedef std::vector<wT> wTv;				//!< wannier type vector
	typedef std::vector<wTv> wi;				//!< wannier indexing
	typedef std::vector<cMat> cMatv;			//!< complex matrix vector
	typedef std::vector<fMat> fMatv;			//!< real matrix vector


	//! operator! for rv
	inline std::vector<bool> operator!(rv inp) noexcept {
		inp.flip(); return inp;
	}
	//! atomic properties
	namespace atom {
		//! atomic mass
		extern const std::map<std::string, double> mass;
	}


	//! raw poscar info
	struct psc {
		fMat B;				//!< basis
		fMat Ap;			//!< atomic positions
		aCv N;				//!< atomic count for each atomic type
		idv id;				//!< id strings for each atomic type
	};
	//! valence, conduction band edges
	struct vb_cb {
		double vb;			//!< valence band edge
		double cb;			//!< conduction band edge
	};
	//! distance and bonds
	struct d_b {
		double val;			//!< distance
		fMat bonds;			//!< bonds
	};
	//! omen material file
	struct omf {
		vb_cb E;			//!< band edges
		std::vector<size_t> Norb;	//!< number of orbitals
		std::vector<double> mass;	//!< atomic masses
		//! number of entries
		inline size_t N() const noexcept { return Norb.size(); }
		//! check whether struct is empty
		inline bool empty() const noexcept { return Norb.empty(); }
	};
	//! omen lattice file
	struct olf {
		fMat B;				//!< basis
		fMat Ap;			//!< atomic positions
		idv id;				//!< id strings for each atomic type
		aTv T;				//!< atomic types
		size_t nn;			//!< number of next neighbors
		double bl;			//!< bond length encompassing nn bonds
		//! dimension of space
		inline size_t dim() const noexcept { return B.M(); }
		//! number of atomic positions
		inline size_t N() const noexcept { return Ap.N(); }
		//! check whether struct is empty
		inline size_t empty() const noexcept { return Ap.empty(); }
	};
	//! sparse hamiltonian element
	struct sphel {
		size_t m;			//!< row index
		size_t n;			//!< column index
		hel h;				//!< hamiltonian element
		
		//! comparison
		inline bool operator<(const sphel& rhs) const noexcept {
			return m<rhs.m ? true: m>rhs.m ? false: n<rhs.n;
		}
		//! equality
		inline bool operator==(const sphel& rhs) const noexcept {
			return m==rhs.m && n==rhs.n;
		}
		//! not equality
		inline bool operator!=(const sphel& rhs) const noexcept {
			return !(*this==rhs);
		}
		//! streaming operator
		friend inline std::ostream& operator<<(std::ostream& os, const sphel& inp) noexcept {
			return (os<<"("<<inp.m<<","<<inp.n<<"): "
				  << std::real(inp.h) << (std::imag(inp.h)<.0 ? "-": "+")
				  << std::abs(std::imag(inp.h)) << "i");
		}
	};
	//! atomic positions and types
	class Ap_T final: public mat_vec_b<fMat,fArray,aTv,aT> {
	public:
		/** @name constructors
		 */
		using mat_vec_b<fMat,fArray,aTv,aT>::mat_vec_b;


		/** @name data access
		 */
		//! atomic positions reference
		inline fMat& Ap() noexcept { return mat_; }
		//! atomic positions const reference
		inline const fMat& Ap() const noexcept { return mat_; }
		//! atomic types reference
		inline aTv& T() noexcept { return vec_; }
		//! atomic types const reference
		inline const aTv& T() const noexcept { return vec_; }
	};
	//! atomic positions and id strings
	class Ap_id final: public mat_vec_b<fMat,fArray,idv,std::string> {
	public:
		/** @name constructors
		 */
		using mat_vec_b<fMat,fArray,idv,std::string>::mat_vec_b;


		/** @name data access
		 */
		//! atomic positions reference
		inline fMat& Ap() noexcept { return mat_; }
		//! atomic positions const reference
		inline const fMat& Ap() const noexcept { return mat_; }
		//! id strings reference
		inline idv& id() noexcept { return vec_; }
		//! id strings const reference
		inline const idv& id() const noexcept { return vec_; }
	};
	//! wannier positions and spread
	class Wp_s final: public mat_vec_b<fMat,fArray,fMat,double> {
	public:
		/** @name constructors
		 */
		using mat_vec_b<fMat,fArray,fMat,double>::mat_vec_b;

		/** @name data access
		 */
		//! Wannier centers reference
		inline fMat& Wp() noexcept { return mat_; }
		//! Wannier centers const reference
		inline const fMat& Wp() const noexcept { return mat_; }
		//! Wannier spread reference
		inline fMat& s() noexcept { return vec_; }
		//! Wannier spread const reference
		inline const fMat& s() const noexcept { return vec_; }
	};
	//! H blocks and R vectors
	template<class MT=cMat>
	class R_H final: public mat_vec_cb<fMat,fArray,std::vector<MT>,MT> {
	public:
		/** @name constructors
		 */
		using mat_vec_cb<fMat,fArray,std::vector<MT>,MT>::mat_vec_cb;
		//! constructor from R vectors and Hamiltonian matrices
		inline R_H(fMat R_, std::vector<MT> H_) noexcept:
			mat_vec_cb<fMat,fArray,std::vector<MT>,MT>(std::move(R_),std::move(H_)) {
			assert(this->empty() ||
			     (!this->empty() && this->front().square() &&
			       std::all_of(this->cbegin(),this->cend(),
				[this](const auto& h)->bool{return this->front().S()==h.S();})));
		}


		/** @name data access
		 */
		//! R vectors const reference
		inline const fMat& R() const noexcept { return this->mat_; }
		//! Hamiltonians const reference
		inline const std::vector<MT>& H() const noexcept { return this->vec_; }
		//! center Hamiltonian const reference (at R=0)
		inline const MT& center() const noexcept { return this->operator[](this->N()/2); }


		/** @name information
		 */
		//! number of Wannier functions
		inline wC Nw() const noexcept { return this->empty() ? 0: this->front().M(); }
		//! number of R vectors
		inline size_t NR() const noexcept { return this->empty() ? 0: this->R().N(); }


		/** @name conversion
		 */
		//! restricted space matrix
		rv r() const noexcept {
			rv res; res.reserve(R().M());
			for (auto i=R().crBegin(),e=R().crEnd(); i!=e; ++i)
				res.push_back(!any(*i));
			return res;
		}
	

		/** @name printing
		 */
		/** write R and H to files.
		 * Filenames are seed_R.bin and seed_H%04d.bin
		 * @param seed the seen for the filenames
		 */
		inline void writeToFile(const std::string& seed) const {
			constexpr size_t NTZ = 4;
			for (size_t i=0; i!=this->size(); ++i) {
				std::stringstream sstr;
				sstr << seed << "_H" << std::setw(NTZ) << std::setfill('0')
				     << (i+1) << ".bin";
				this->vec_[i].writeToFile(sstr.str());
			}
			this->R().writeToFile(seed+"_R.bin");
		}


		/** @name friends
		 */
		//! test friend class
		friend class ::test_io;
		//! test friend class
		friend class ::ll_writerMEM<MT>;
	};
	//! wannier 90 transformation matrices
	class k_U final: public mat_vec_cb<fMat,fArray,std::vector<cMat>,cMat> {
	public:
		/** @name constructors
		 */
		using mat_vec_cb<fMat,fArray,std::vector<cMat>,cMat>::mat_vec_cb;
		//! constructor from k points and U matrices
		inline k_U(fMat k_, std::vector<cMat> U_) noexcept:
			mat_vec_cb<fMat,fArray,std::vector<cMat>,cMat>(std::move(k_),std::move(U_)) {
			assert(empty() ||
			(!empty() && std::all_of(cbegin(),cend(),
			[this](const auto& u)->bool{return front().S()==u.S();})));
		}


		/** @name data access
		 */
		//! k points const reference
		inline const fMat& k() const noexcept { return mat_; }
		//! U matrices const reference
		inline const std::vector<cMat>& U() const noexcept { return vec_; }


		/** @name information
		 */
		//! number of bands included in the Wannierization
		inline size_t Nb() const noexcept { return empty() ? 0: front().M(); }
		//! number of Wannier functions
		inline wC Nw() const noexcept { return empty() ? 0: front().N(); }
	};
	//! path and position
	class p_p final: public mat_vec_cb<fMat,fArray,fMat,double> {
	public:
		/** @name constructors
		 */
		using mat_vec_cb<fMat,fArray,fMat,double>::mat_vec_cb;
		/** constructor from path, position and number of points
		 * @param path_ the actual points along a path
		 * @param pos_ the distance from the first point for each point on the path
		 * @param Nps_ vector holding the number points between each fixpoint
		 */
		inline p_p(fMat path_, fMat pos_, std::vector<size_t> Nps_) noexcept:
			mat_vec_cb<fMat,fArray,fMat,double>(std::move(path_),std::move(pos_)),
			Nps_(std::move(Nps_)) {
			assert(N()==std::accumulate(Nps().cbegin(),Nps().cend(),size_t(0)));
		}


		/** @name data access
		 */
		//! path const reference
		inline const fMat& path() const noexcept { return mat_; }
		//! position const reference
		inline const fMat& pos() const noexcept { return vec_; }
		//! number of points const reference
		inline const std::vector<size_t>& Nps() const noexcept { return Nps_; }


	protected:
		/** @name member variables
		 */
		std::vector<size_t> Nps_;	//!< Number of points
	};
}

#endif // _LL_TYPES_

/** @}
 */
