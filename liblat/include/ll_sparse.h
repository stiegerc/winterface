// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_SPARSE_
#define _LL_SPARSE_

#include "ll_types.h"
#include "ll_compound.h"
#include "aux_io.h"
#include "aux_parser.h"
#include <string>
#include <fstream>


class ll_sparse;
namespace ll__ {
	R_H<ll_sparse> readHrSparse(const std::string& fileName);
}

template <class MT> class ll_writerMEM;
template <class MT> class ll_writerR0;


/** simple sparse matrix class. This is essentially a wrapper for a sorted std::vector
 * holding bond indices and interaction data {i1,i2,h}. It is mostly used to read
 * from files and to do comparisons
 */
class ll_sparse final: public ll__::vec_cb<std::vector<ll__::sphel>,ll__::sphel> {
public:
	/** @name types
	 */
	typedef ll__::hel hel;		//!< Hamiltonian element
	typedef ll__::sphel el;		//!< sparse matrix element


	/** @name constructors
	 */
	//! default constructor
	inline explicit ll_sparse(const size_t M=0): ll_sparse(M,M) {}
	/** constructor from number of rows and columns
	 * @param M number of rows
	 * @param N number of columns
	 */
	inline explicit ll_sparse(const size_t M, const size_t N): M_(M), N_(N) {}
	/** constructor from a file.
	 * @param fileName name of the file
	 * @param mode opening mode
	 */
	inline explicit ll_sparse(const std::string& fileName, const std::string& mode="OMEN") {
		using namespace aux;

		switch (fnvHash(mode.c_str())) {
			case "OMEN"_h:
			{
				auto file = aux::openFile<std::ifstream>(fileName, std::ios::binary);

				// header
				double head[3];
				file.read((char*) &head, 3*sizeof(double));
				M_ = size_t(head[0]); N_ = M_;
				
				this->vec_.reserve(size_t(head[1]));
				const size_t ip=head[2];

				// data
				double el_[4];
				while (this->vec_.size()!=this->vec_.capacity()) {
					file.read((char*) &el_, 4*sizeof(double));
					this->vec_.push_back({size_t(el_[0])-ip,size_t(el_[1])-ip,{el_[2],el_[3]}});
				}

				file.close();
				assert(std::is_sorted(begin(),end(),[](const el& i, const el& j)->bool
					{return i.m==j.m ? i.n<j.n: i.m<j.m;}));
			}
			break;
			case "text"_h:
			{
				auto file = aux::openFile<std::ifstream>(fileName);
				std::string line;

				// M,N,size
				{
					uint32_t M,N,S;
					std::getline(file,line);			
					sscanf(line.c_str(),"%u %u %u",&M,&N,&S);
					M_ = M, N_ = N; this->vec_.reserve(S);
				}
				
				// read sparse data
				uint32_t m,n; double re,im;
				while (this->vec_.size()<this->vec_.capacity()) {
					std::getline(file,line);
					sscanf(line.c_str(),"%u %u %lf %lf",&m,&n,&re,&im);
					this->vec_.push_back({m,n,{re,im}});
				}
			}
			break;
			default:
				throw(std::invalid_argument("sparse file mode \'"+mode+"\' invalid"));
		}
	}


	/** @name information
	 */
	//! number of rows
	inline size_t M() const noexcept { return N_; }
	//! number of columns
	inline size_t N() const noexcept { return N_; }
	//! number of total elements
	inline size_t L() const noexcept { return M()*N(); }
	//! number of rows and columns
	inline lm__::lm_size S() const noexcept { return {M(),N()}; }
	//! number of total elements
	inline size_t size() const noexcept { return M()*N(); }
	//! number of non-zero entries
	inline size_t nnz() const noexcept { return this->vec_.size(); }
	//! returns whether this is a square matrix
	inline bool square() const noexcept { return M()==N(); }


	/** data access
	 */
	/** indexed data access
	 * @param m row index
	 * @param n column index
	 */
	inline hel operator()(const size_t m, const size_t n) const noexcept {
		assert(m<M());
		assert(n<N());
	
		struct ipair { size_t m,n; };
		const auto itr = std::lower_bound(cbegin(),cend(),ipair{m,n},
			[](const ll__::sphel& h, const ipair& j)->bool{return h.m==j.m ? h.n<j.n: h.m<j.m;});
		return (itr!=cend() && (itr->m==m && itr->n==n)) ? itr->h: 0.0;
	}
	/** linear indexed data access
	 * @param i linear index
	 */
	inline hel operator[](const size_t i) const noexcept {
		assert(i<L());
		return this->operator()(i/M(),i%M());
	}


	/** comparison
	 */
	//! equality operator
	inline bool operator==(const ll_sparse& rhs) const noexcept {
		if (rhs.M()!=M() || rhs.nnz()!= nnz()) return false;
		for (auto i=cbegin(),e=cend(),j=rhs.cbegin(); i!=e; ++i,++j)
			if (*i!=*j) return false;
		return true;
	}
	//! inequality operator
	inline bool operator!=(const ll_sparse& rhs) const noexcept {
		return !(*this==rhs);
	}


	/** @name modification
	 */
	//! transpose the matrix
	inline const ll_sparse& T() noexcept {
		for (auto& i: *this) {
			std::swap(i.m,i.n);
			i.h = std::conj(i.h);
		}
		std::sort(begin(),end());
		return *this;
	}


	/** @name conversion
	 */
	/** function to extract a submatrix
	 * @param m row starting index
	 * @param n column starting index
	 * @param lm number of rows to include
	 * @param ln number of columns to include
	 */
	inline ll_sparse get(const size_t m, const size_t n,
			     const size_t lm=NPOS__, const size_t ln=NPOS__) const noexcept {
		assert(m<M());
		assert(n<N());

		const size_t ubm = (m+lm>M()) ? M(): m+lm;
		const size_t ubn = (n+ln>N()) ? N(): n+ln;
		ll_sparse res(ubn-n);
		
		for (auto i = std::lower_bound(cbegin(),cend(),m,
			[](const el& i, const size_t m)->bool{return i.m<m;}),
			  e = std::lower_bound(cbegin(),cend(),ubm,
			[](const el& i, const size_t m)->bool{return i.m<m;});
			  i!=e; ++i) {
		
			if (res.vec_.capacity()==res.vec_.size())
				res.vec_.reserve(res.vec_.size()+res.N());
			if (i->n>=n && i->n<ubn) res.vec_.push_back({i->m-m,i->n-n,i->h});
		}

		return res;
	}
	//! conversion to full matrix
	inline ll__::cMat full() const noexcept {
		ll__::cMat res(M(),N());
		for (const auto& i: vec_)
			res(i.m,i.n) = i.h;
		return res;
	}


	/** @name printing
	 */
	/** function to print the matrix to a textfile
	 * @param fileName name of the file
	 * @param prec floating point precision
	 */
	inline void printToFile(const std::string& fileName, const size_t prec) const {
		auto file = aux::openFile<std::ofstream>(fileName);
		file << M() << " " << N() << " " << this->vec_.size() << "\n";
		for (const auto& i: this->vec_)
			file << std::setw(8) << i.m
			     << std::setw(8) << i.n
			     << std::setw(prec+8) << std::fixed
			     << std::setprecision(prec) << std::real(i.h)
			     << std::setw(prec+8) << std::fixed
			     << std::setprecision(prec) << std::imag(i.h);
	}
	//! streaming operator
	friend inline std::ostream& operator<<(std::ostream& os, const ll_sparse& inp) noexcept {
		for (const auto& i: inp)
			os << i << "\n";
		return os;
	}


protected:
	/** @name member variables
	 */
	size_t M_;	//!< number of rows
	size_t N_;	//!< number of columns


	/** @name friends
	 */
	//! read sparse from file into R_H container function
	friend ll__::R_H<ll_sparse> ll__::readHrSparse(const std::string& fileName);
	//! writer to memory
	friend class ll_writerMEM<ll_sparse>;
	//! writer to memory for only R=0
	friend class ll_writerR0<ll_sparse>;
};

#endif // _LL_SPARSE_

/** @}
 */
