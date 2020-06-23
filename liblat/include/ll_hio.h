// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_HIO_
#define _LL_HIO_

#include "lm_fn.h"
#include "ll_fn.h"
#include "ll_types.h"
#include "ll_compound.h"
#include "ll_sparse.h"
#include "ll_hbondss.h"
#include <string>
#include <fstream>
#include <cassert>
#include <cfloat>
#include <omp.h>



//! base Hamiltonian writer virtual class
class ll_writer {
public:
	/** @name types
	 */
	typedef std::vector<ll__::sphel> hdat;	//!< Hamiltonian data
	
public:
	/** @name contructors
	 */
	/** constructor from R vectors
	 * @param R R vectors
	 */
	inline explicit ll_writer(lm__::fMat R) noexcept:
			R_(std::move(R)), r_(R_.ccBegin()),
			nnz_(R_.N(),0), n_(nnz_.begin()) {
		assert(lm__::all(R_.eq(lm__::round(R_))));
		assert(std::is_sorted(R_.ccBegin(),R_.ccEnd()));
	}
	//! virtual destructor
	virtual ~ll_writer() noexcept {};


	/** @name general information
	 */
	//! dimension of space
	inline size_t dim() const noexcept { return R_.M(); }
	//! number of R vectors
	inline size_t NR() const noexcept { return R_.N(); }
	//! R vectors const reference
	inline const lm__::fMat& R() const noexcept { return R_; }
	//! number of non-zero elements written
	inline const std::vector<uint64_t>& nnz() const noexcept { return nnz_; }
	//! nnz consistency test, i.e. nnz(R)==nnz(-R)
	inline bool nnzTest() const noexcept {
		for (auto i=R_.ccBegin(), e=R_.ccEnd(); i!=e; ++i) {
			const auto j = std::lower_bound(R_.ccBegin(),R_.ccEnd(),-(*i));
			if (j!=R_.ccEnd() && nnz_[size_t(i)]!=nnz_[size_t(j)])
				return false;
		}
		return true;
	}
	//! end of file, i.e. all Hamiltonian blocks written
	inline bool eof() const noexcept { return r_==R_.ccEnd(); }


	/** @name current block information
	 */
	//! returns the index of the current block
	inline size_t c_ind() const noexcept { return size_t(r_); }
	//! returns the current R vector
	inline const lm__::fCol c_R() const noexcept { return *r_; }
	//! returns the nnz of the current block
	inline uint64_t c_nnz() const noexcept { return *n_; }
	//! returns the bytes written in the current block
	inline size_t c_bytes() const noexcept { return c_nnz()*bytes_(); }
	//! return a description of the current block
	virtual std::string c_descr() const noexcept=0;


	/** @name data insertion
	 */
	//! write Hamiltonian data
	inline void insert(const hdat& inp) {
		*n_ += inp.size();
		insert_(inp);
	}


	/** @name initiate and finalize
	 */
	//! new Hamiltonian block initialization
	virtual void newBlock()=0;
	//! flush current Hamiltonian block function
	inline void flush() {
		flush_();
		++r_, ++n_;
	}


protected:
	/** @name internals
	 */
	//! the number of bytes for each Hamiltonian entry
	virtual size_t bytes_() const noexcept=0;
	/** write data blob function
	 * @param inp vector of Hamiltonian entries
	 */
	virtual void insert_(const hdat& inp)=0;
	//! flush current block Hamiltonian function
	virtual void flush_()=0;
	
	
protected:
	/** @name member variables
	 */
	lm__::fMat R_;				//!< R vectors
	lm__::c_fColItr r_;			//!< iterator to current R vector
	std::vector<uint64_t> nnz_;		//!< number of entries written vector
	std::vector<uint64_t>::iterator n_;	//!< iterator to current number of entries
};


//! write to memory writer
template <class MT>
class ll_writerMEM final: public ll_writer {
public:
	/** @name constructors
	 */
	/** constructors from R vectors and vector of matrices
	 * @param R R vectors
	 * @param H vector of Hamiltonian matrices to write into
	 */
	inline explicit ll_writerMEM(lm__::fMat R, std::vector<MT>& H):
			ll_writer(std::move(R)), m_(H.begin()) {
		assert(H.size()==NR());
	}
	/** constructor from R_H struct
	 * @param hr R_H struct
	 */
	inline explicit ll_writerMEM(ll__::R_H<MT>& hr):
			ll_writer(hr.R()), m_(hr.begin()) {
		assert(!hr.empty());
	}


	/** @name current block information
	 */
	//! description of the current block
	inline std::string c_descr() const noexcept { return lm__::T(this->c_R()).print(0); }


	/** @name initiate and finalize
	 */
	//! new Hamiltonian block initialization
	inline void newBlock() {}


private:
	/** @name internals
	 */
	//! the number of bytes for each Hamiltonian entry
	inline size_t bytes_() const noexcept { return 2*sizeof(size_t)+sizeof(ll__::hel); }
	/** write data blob function
	 * @param inp vector of Hamiltonian entries
	 */
	inline void insert_(const hdat& inp) {
		for (const auto& i: inp)
			(*m_)(i.m,i.n) = i.h;
	}
	//! flush current Hamiltonian block function
	inline void flush_() { ++m_; }


private:
	/** member variables
	 */
	typename std::vector<MT>::iterator m_;	//!< iterator to the current matrix
};
//! the number of bytes for each Hamiltonian entry, overload for real matrices
template<> inline size_t ll_writerMEM<ll__::fMat>::bytes_() const noexcept {
	return 2*sizeof(size_t)+sizeof(double);
}
/** write data blob function, overload for real matrices
 * @param inp vector of Hamiltonian entries
 */
template<> inline void ll_writerMEM<lm__::fMat>::insert_(const hdat& inp) {
	for (const auto& i: inp)
		(*m_)(i.m,i.n) = std::real(i.h);
}
/** write data blob function, overload for sparse matrices
 * @param inp vector of Hamiltonian entries
 */
template<> inline void ll_writerMEM<ll_sparse>::insert_(const hdat& inp) {
	m_->vec_.reserve(m_->vec_.size()+inp.size());
	for (const auto& i: inp)
		m_->vec_.push_back(i);	
}


//! write to memory writer, only R=0 version
template <class MT>
class ll_writerR0 final: public ll_writer {
public:
	/** @name constructors
	 */
	/** constructor from matrix
	 * @param H Hamiltonian matrix
	 * @param dim number of dimensions
	 */
	inline explicit ll_writerR0(MT& H, const size_t dim=DIM__):
			ll_writer(lm__::zeros<lm__::fMat>(dim,1)), H_(H) {}


	/** @name current block information
	 */
	//! description of the current block
	inline std::string c_descr() const noexcept { return lm__::T(this->c_R()).print(0); }


	/** @name initiate and finalize
	 */
	//! new Hamiltonian block initialization
	inline void newBlock() {}


private:
	/** @name internals
	 */
	//! the number of bytes for each Hamiltonian entry
	inline size_t bytes_() const noexcept { return 2*sizeof(size_t)+sizeof(ll__::hel); }
	/** write data blob function
	 * @param inp vector of Hamiltonian entries
	 */
	inline void insert_(const hdat& inp) {
		for (const auto& i: inp)
			H_(i.m,i.n) = i.h;
	}
	//! flush current Hamiltonian block function
	inline void flush_() {}

private:
	MT& H_;
};
//! the number of bytes for each Hamiltonian entry, overload for real matrices
template<> inline size_t ll_writerR0<ll__::fMat>::bytes_() const noexcept {
	return 2*sizeof(size_t)+sizeof(double);
}
/** write data blob function, overload for real matrices
 * @param inp vector of Hamiltonian entries
 */
template<> inline void ll_writerR0<lm__::fMat>::insert_(const hdat& inp) {
	for (const auto& i: inp)
		H_(i.m,i.n) = std::real(i.h);
}
/** write data blob function, overload for sparse matrices
 * @param inp vector of Hamiltonian entries
 */
template<> inline void ll_writerR0<ll_sparse>::insert_(const hdat& inp) {
	H_.vec_.reserve(H_.vec_.size()+inp.size());
	for (const auto& i: inp)
		H_.vec_.push_back(i);	
}


//! wannier90 hr style text writer
class ll_writerW90 final: public ll_writer {
public:
	/** @name constructors
	 */
	/** constructor from filename
	 * @param fileName name of the file
	 * @param R R vectors
	 * @param Nw number of Wannier functions
	 * @param prec floating point precision
	 */
	inline ll_writerW90(const std::string& fileName, lm__::fMat R,
					const size_t Nw, const size_t prec=12):
			ll_writer(std::move(R)),
			fileName_(fileName), file_(aux::openFile<std::ofstream>(fileName)),
			prec_(prec) {
	
		// write header
		{
			std::locale::global(std::locale("en_US.utf8"));
			std::time_t t = std::time(nullptr);
			char mbstr[100];
			std::strftime(mbstr, sizeof(mbstr), "%c", std::localtime(&t));
			
			file_ << mbstr << "\n";
		}
		
		// write Nw, Nr
		file_ << std::setw(precision()) << Nw << "\n";
		file_ << std::setw(precision()) << NR() << "\n";
		
		// fake degeneracy
		for (size_t i=1; i<=NR(); ++i)
			file_ << std::setw(5) << 1 << (!(i%15) && i!=NR() ? "\n": "");
		file_ << "\n";
	}

	
	/** @name general information
	 */
	//! floating point precision
	inline size_t precision() const noexcept { return prec_; }
	//! filename
	inline const std::string fileName() const noexcept { return fileName_; }
	

	/** @name current block information
	 */
	//! description of the current block
	inline std::string c_descr() const noexcept { return lm__::T(this->c_R()).print(0); }


	/** @name initiate and finalize
	 */
	//! new Hamiltonian block initialization
	inline void newBlock() {}


private:
	/** @name internals
	 */
	//! the number of bytes for each Hamiltonian entry
	inline size_t bytes_() const noexcept { return 2*sizeof(size_t)+sizeof(ll__::hel); }
	/** write data blob function
	 * @param inp vector of Hamiltonian entries
	 */
	inline void insert_(const hdat& inp) {
		for (const auto& h: inp) {
			for (const int i: c_R())
				file_ << std::setw(8) << i;
			file_ << std::setw(8) << h.m
			      << std::setw(8) << h.n
			      << std::setw(precision()+8) << std::fixed
			      << std::setprecision(precision()) << std::real(h.h)
			      << std::setw(precision()+8) << std::fixed
			      << std::setprecision(precision()) << std::imag(h.h) << "\n";
		}
	}
	//! flush current Hamiltonian block function
	inline void flush_() noexcept {}


private:
	/** @name member variables
	 */
	std::string fileName_;		//!< filename
	std::ofstream file_;		//!< filestream
	const size_t prec_;		//!< floating point precision
};


//! binary all in one file hamiltonian writer
template <class IT, class FT>
class ll_writerBIN final: public ll_writer {
public:
	/** @name types
	 */
	typedef typename ll_writer::hdat hdat;	//!< Hamiltonian data

public:
	/** @name constructors
	 */
	/** constructor from filename
	 * @param fileName name of the file
	 * @param R R vectors
	 * @param Nw number of Wannier functions
	 */
	inline ll_writerBIN(const std::string& fileName, lm__::fMat R, const size_t Nw):
			ll_writer(std::move(R)),
			fileName_(fileName), file_(aux::openFile<std::ofstream>(fileName,std::ios::binary)),
       			spos_(this->R_.N()), s_(spos_.begin()) {
		
		// write header
		file_.write((char*) header_(), 5*sizeof(char));
	
		// write dim, Nw, NR
		uint32_t b32_;
		b32_ = this->R_.M(), file_.write((char*) &b32_, sizeof(uint32_t));
		b32_ = Nw,           file_.write((char*) &b32_, sizeof(uint32_t));
		b32_ = this->R_.N(), file_.write((char*) &b32_, sizeof(uint32_t));
	}
	//! destructor
	inline ~ll_writerBIN() {
		// reopen file
		file_.close();
		file_.open(fileName(), std::ios::binary | std::ios::out | std::ios::in);

		// write real nnz's
		auto n_ = this->nnz_.cbegin();
		for (const auto s: spos_)
			file_.seekp(s), file_.write((char*) &(*n_), sizeof(uint64_t)), ++n_;
	}
	

	/** @name general information
	 */
	//! filename
	inline const std::string fileName() const noexcept { return fileName_; }
	

	/** @name current block information
	 */
	//! description of the current block
	inline std::string c_descr() const noexcept { return lm__::T(this->c_R()).print(0); }


	/** @name initiate and finalize
	 */
	//! new Hamiltonian block initialization
	inline void newBlock() {
		// write current R to file
		for (const auto& r: this->c_R())
			file_.write((char*) &r, sizeof(double));

		// store stream position of nnz
		*s_ = file_.tellp();

		// write nnz placeholder
		const uint64_t plh = 0;
		file_.write((char*) &plh, sizeof(uint64_t));
	}

private:
	/** @name internals
	 */
	//! the number of bytes for each Hamiltonian entry
	inline size_t bytes_() const noexcept { return 2*sizeof(IT)+sizeof(FT); }
	/** write data blob function
	 * @param inp vector of Hamiltonian entries
	 */
	inline void insert_(const hdat& inp) { assert(false); }
	//! flush current Hamiltonian block function
	inline void flush_() noexcept { ++s_; }
	//! header string
	static inline const char* header_() noexcept { return "inval"; }

private:
	/** @name member variables
	 */
	std::string fileName_;				//!< filename
	std::ofstream file_;				//!< filestream
	std::vector<std::streampos> spos_;		//!< stream positions
	std::vector<std::streampos>::iterator s_;	//!< current stream position
};

//! insert specialization for 16bit indices and single precision Hamiltonian entries
template<> inline void ll_writerBIN<uint16_t,float>::insert_(const hdat& inp) {
	for (const auto& i: inp) {
		const uint16_t nm[] = {uint16_t(i.m),uint16_t(i.n)};
		file_.write((char*) &nm, 2*sizeof(uint16_t));
		const float h = std::real(i.h);
		file_.write((char*) &h, sizeof(float));
	}
}
//! insert specialization for 32bit indices and single precision Hamiltonian entries
template<> inline void ll_writerBIN<uint32_t,float>::insert_(const hdat& inp) {
	for (const auto& i: inp) {
		const uint32_t nm[] = {uint32_t(i.m),uint32_t(i.n)};
		file_.write((char*) &nm, 2*sizeof(uint32_t));
		const float h = std::real(i.h);
		file_.write((char*) &h, sizeof(float));
	}
}
//! insert specialization for 16bit indices and single precision complex Hamiltonian entries
template<> inline void ll_writerBIN<uint16_t,std::complex<float>>::insert_(const hdat& inp) {
	for (const auto& i: inp) {
		const uint16_t nm[] = {uint16_t(i.m),uint16_t(i.n)};
		file_.write((char*) &nm, 2*sizeof(uint16_t));
		const std::complex<float> h = {float(std::real(i.h)),float(std::imag(i.h))};
		file_.write((char*) &h, sizeof(std::complex<float>));
	}
}
//! insert specialization for 32bit indices and single precision complex Hamiltonian entries
template<> inline void ll_writerBIN<uint32_t,std::complex<float>>::insert_(const hdat& inp) {
	for (const auto& i: inp) {
		const uint32_t nm[] = {uint32_t(i.m),uint32_t(i.n)};
		file_.write((char*) &nm, 2*sizeof(uint32_t));
		const std::complex<float> h = {float(std::real(i.h)),float(std::imag(i.h))};
		file_.write((char*) &h, sizeof(std::complex<float>));
	}
}
//! insert specialization for 16bit indices and double precision Hamiltonian entries
template<> inline void ll_writerBIN<uint16_t,double>::insert_(const hdat& inp) {
	for (const auto& i: inp) {
		const uint16_t nm[] = {uint16_t(i.m),uint16_t(i.n)};
		file_.write((char*) &nm, 2*sizeof(uint16_t));
		const double h = std::real(i.h);
		file_.write((char*) &h, sizeof(double));
	}
}
//! insert specialization for 32bit indices and double precision Hamiltonian entries
template<> inline void ll_writerBIN<uint32_t,double>::insert_(const hdat& inp) {
	for (const auto& i: inp) {
		const uint32_t nm[] = {uint32_t(i.m),uint32_t(i.n)};
		file_.write((char*) &nm, 2*sizeof(uint32_t));
		const double h = std::real(i.h);
		file_.write((char*) &h, sizeof(double));
	}
}
//! insert specialization for 16bit indices and double precision complex Hamiltonian entries
template<> inline void ll_writerBIN<uint16_t,std::complex<double>>::insert_(const hdat& inp) {
	for (const auto& i: inp) {
		const uint16_t nm[] = {uint16_t(i.m),uint16_t(i.n)};
		file_.write((char*) &nm, 2*sizeof(uint16_t));
		file_.write((char*) &i.h, sizeof(std::complex<double>));
	}
}
//! insert specialization for 32bit indices and double precision complex Hamiltonian entries
template<> inline void ll_writerBIN<uint32_t,std::complex<double>>::insert_(const hdat& inp) {
	for (const auto& i: inp) {
		const uint32_t nm[] = {uint32_t(i.m),uint32_t(i.n)};
		file_.write((char*) &nm, 2*sizeof(uint32_t));
		file_.write((char*) &i.h, sizeof(std::complex<double>));
	}
}

//! header specializations for 16bit indices and single precision Hamiltonian entries
template<> inline const char* ll_writerBIN<uint16_t,float>::header_() noexcept { return "hrssr"; }
//! header specializations for 32bit indices and single precision Hamiltonian entries
template<> inline const char* ll_writerBIN<uint32_t,float>::header_() noexcept { return "hrlsr"; }
//! header specializations for 16bit indices and double precision Hamiltonian entries
template<> inline const char* ll_writerBIN<uint16_t,double>::header_() noexcept { return "hrsdr"; }
//! header specializations for 32bit indices and double precision Hamiltonian entries
template<> inline const char* ll_writerBIN<uint32_t,double>::header_() noexcept { return "hrldr"; }
//! header specializations for 16bit indices and single precision complex Hamiltonian entries
template<> inline const char* ll_writerBIN<uint16_t,std::complex<float>>::header_() noexcept { return "hrssc"; }
//! header specializations for 32bit indices and single precision complex Hamiltonian entries
template<> inline const char* ll_writerBIN<uint32_t,std::complex<float>>::header_() noexcept { return "hrlsc"; }
//! header specializations for 16bit indices and double precision complex Hamiltonian entries
template<> inline const char* ll_writerBIN<uint16_t,std::complex<double>>::header_() noexcept { return "hrsdc"; }
//! header specializations for 32bit indices and double precision complex Hamiltonian entries
template<> inline const char* ll_writerBIN<uint32_t,std::complex<double>>::header_() noexcept { return "hrldc"; }




namespace ll__ {

	/* R grid related
	 */
	/** probe R grid function.
	 * This function attempts to find the R vectors for which Hamiltonian entries exist
	 * by scanning over a set of adjacent unit cells.
	 * @param W wbh (ll_hbonds or ll_hbondss)
	 * @param B basis
	 * @param Ap atomic positions
	 * @param T atomic types
	 * @param r restriction vector
	 * @param IR interaction radius
	 * @param approx switch to enable approximate interactions (true) or enforce exact matching (false)
	 */
	template <class WT>
	lm__::fMat getConnectedGrid(const WT& W, const lm__::fMat& B, const lm__::fMat& Ap, const aTv& T,
				const rv& r, const double IR, const bool approx) noexcept;
	/** probe R grid function, use cell in wbh version.
	 * This function attempts to find the R vectors for which Hamiltonian entries exist
	 * by scanning over a set of adjacent unit cells.
	 * @param W wbh (ll_hbonds or ll_hbondss)
	 * @param r restriction vector
	 * @param IR interaction radius
	 * @param approx switch to enable approximate interactions (true) or enforce exact matching (false)
	 */
	inline lm__::fMat getConnectedGrid(const ll_hbonds& W, const rv& r,
				const double IR, const bool approx) noexcept {
		return getConnectedGrid(W,W.cell().B(),W.cell().Ap(),W.cell().type(),r,IR,approx);
	}


	/* hamiltonian constructors
	 */
	/** construct Hamiltonian matrices.
	 * @param W wbh (ll_hbonds or ll_hbondss)
	 * @param B basis
	 * @param Ap atomic positions
	 * @param T atomic types
	 * @param writer Hamiltonian writer
	 * @param IR interaction radius
	 * @param approx switch to enable approximate interactions (true) or enforce exact matching (false)
	 * @param Nthreads number of threads to use in parallel sections
	 * @param verbosity bitmask controlling the verbosity level
	 * @param os stream to print into
	 */
	template <class WT>
	void hctor(const WT& W, const fMat& B, const fMat& Ap, const aTv& T,
			  ll_writer& writer, const double IR, const bool approx,
			  const size_t Nthreads=0, const size_t verbosity=0, std::ostream& os=std::cout);
	/** construct Hamiltonian matrices, writer rvalue version.
	 * @param W wbh (ll_hbonds or ll_hbondss)
	 * @param B basis
	 * @param Ap atomic positions
	 * @param T atomic types
	 * @param writer Hamiltonian writer
	 * @param IR interaction radius
	 * @param approx switch to enable approximate interactions (true) or enforce exact matching (false)
	 * @param Nthreads number of threads to use in parallel sections
	 * @param verbosity bitmask controlling the verbosity level
	 * @param os stream to print into
	 */
	template <class WT>
	inline void hctor(const WT& W, const fMat& B, const fMat& Ap, const aTv& T,
			ll_writer&& writer, const double IR, const bool approx,
			const size_t Nthreads=0, const size_t verbosity=0, std::ostream& os=std::cout) {
		hctor(W,B,Ap,T,writer,IR,approx,Nthreads,verbosity,os);
	}
	/** construct Hamiltonian matrices, cell from wbh version
	 * @param W wbh (ll_hbonds or ll_hbondss)
	 * @param writer Hamiltonian writer
	 * @param IR interaction radius
	 * @param approx switch to enable approximate interactions (true) or enforce exact matching (false)
	 * @param Nthreads number of threads to use in parallel sections
	 * @param verbosity bitmask controlling the verbosity level
	 * @param os stream to print into
	 */
	inline void hctor(const ll_hbonds& W, ll_writer& writer, const double IR, const bool approx,
			const size_t Nthreads=0, const size_t verbosity=0, std::ostream& os=std::cout) {
		hctor(W,W.cell().B(),W.cell().getcAp(),W.cell().type(),writer,IR,approx,Nthreads,verbosity,os);
	}
	/** construct Hamiltonian matrices, cell from wbh, writer rvalue version
	 * @param W wbh (ll_hbonds or ll_hbondss)
	 * @param writer Hamiltonian writer
	 * @param IR interaction radius
	 * @param approx switch to enable approximate interactions (true) or enforce exact matching (false)
	 * @param Nthreads number of threads to use in parallel sections
	 * @param verbosity bitmask controlling the verbosity level
	 * @param os stream to print into
	 */
	inline void hctor(const ll_hbonds& W, ll_writer&& writer, const double IR, const bool approx,
			const size_t Nthreads=0, const size_t verbosity=0, std::ostream& os=std::cout) {
		hctor(W,W.cell().B(),W.cell().getcAp(),W.cell().type(),writer,IR,approx,Nthreads,verbosity,os);
	}


	//! read hr sparse files as written by ll_writerBIN
	R_H<ll_sparse> readHrSparse(const std::string& fileName);


	/** adapt ll_hbondss to structure.
	 * @param W wbh
	 * @param inp user input
	 * @param Ap atomic positions
	 * @param T atomic types
	 * @param os stream to print into
	 */
	void adaptWBH(ll_hbondss& W, const ll_hbondss_input& inp,
			const fMat& Ap, const aTv& T, std::ostream& os);
	/** adapt ll_hbondss to structure, ll_hbonds version (it does nothing)
	 * @param W wbh
	 * @param inp user input
	 * @param Ap atomic positions
	 * @param T atomic types
	 * @param os stream to print into
	 */
	inline void adaptWBH(ll_hbonds& W, const ll_hbondss_input& inp,
				const fMat& Ap, const aTv& T, std::ostream& os) {}

} // namespace ll__

#endif // _LL_HIO_

/** @}
 */
