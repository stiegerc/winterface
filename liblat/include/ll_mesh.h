// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_MESH_
#define _LL_MESH_

#include "ll_types.h"
#include "lm_types.h"
#include "lm_tVecItr.h"
#include "lm_tCol.h"
#include "lm_tMat.h"
#include "aux_io.h"
#include "ll_lambda.h"
#include <numeric>
#include <array>
#include <iterator>


namespace ll__ {


/** cat vectors and ints into vectors at compile time
 * @{
 */
template<typename... Args> std::vector<size_t> cat(std::vector<size_t>,Args...) noexcept;
template<typename... Args> std::vector<size_t> cat(size_t,Args...) noexcept;
std::vector<size_t> cat(std::vector<size_t>) noexcept;
std::vector<size_t> cat(size_t) noexcept;	
template<typename... Args>
inline std::vector<size_t> cat(std::vector<size_t> inp, Args... args) noexcept {
	for (const auto i: cat(args...))
		inp.push_back(i);
	return inp;
}
template<typename... Args>
inline std::vector<size_t> cat(size_t inp, Args... args) noexcept {
	std::vector<size_t> res = {inp};
	for (const auto i: cat(args...))
		res.push_back(i);
	return res;
}
inline std::vector<size_t> cat(size_t inp) noexcept { return {inp}; }
inline std::vector<size_t> cat(std::vector<size_t> inp) noexcept { return inp; }
/** @} */


/** add, subtract ints from all elements of vectors
 * @{
 */
inline std::vector<size_t> operator+(std::vector<size_t> lhs, size_t rhs) noexcept {
	std::transform(lhs.cbegin(),lhs.cend(),lhs.begin(),
		[rhs](const size_t i)->size_t{return i+rhs;});
	return lhs;
}
inline std::vector<size_t> operator+(size_t lhs, std::vector<size_t> rhs) noexcept {
	std::transform(rhs.cbegin(),rhs.cend(),rhs.begin(),
		[lhs](const size_t i)->size_t{return i+lhs;});
	return rhs;
}
inline std::vector<size_t> operator-(std::vector<size_t> lhs, size_t rhs) noexcept {
	std::transform(lhs.cbegin(),lhs.cend(),lhs.begin(),
		[rhs](const size_t i)->size_t{return i-rhs;});
	return lhs;
}
inline std::vector<size_t> operator-(size_t lhs, std::vector<size_t> rhs) noexcept {
	std::transform(rhs.cbegin(),rhs.cend(),rhs.begin(),
		[lhs](const size_t i)->size_t{return -i+lhs;});
	return rhs;
}
/** @} */


namespace mesh {

/** mesh iterator base class.
 * A mesh iterator traverses a mesh along a 'line' in the mesh,
 * which is implemented by a matrix underneath. This matrix can
 * be traversed column-wise in a manner similar to along columns
 * or rows with respect to the underlyng linear data. Since a
 * matrix is essentially a 2D grid, it can be traversed along
 * rows or columns (or 'lines') by an iterator with a custom step.
 * A mesh iterator works analogous to lm_tItr with a custom step
 * hopping from column to column along 'lines' in the mesh. \n
 * Template arguments are: \n
 * - IT: column iterator
 * - RT: return type
 * - RT_: other return type
 */
template<class IT, class RT, class RT_>
class itr_b: public std::iterator<
		typename std::iterator_traits<IT>::iterator_category,
		typename std::iterator_traits<IT>::value_type,
		typename std::iterator_traits<IT>::difference_type,
		typename std::iterator_traits<IT>::pointer,
		typename std::iterator_traits<IT>::reference> {
public:
	/** @name types
	 */
	typedef typename std::iterator_traits<IT>::value_type val;	//!< value type
	typedef typename std::iterator_traits<IT>::difference_type diff;//!< difference type
	typedef typename std::iterator_traits<IT>::pointer ptr;		//!< pointer type
	typedef typename std::iterator_traits<IT>::reference ref;	//!< reference type


public:
	/** @name constructors
	 */
	//! constructor from iterator and step
	inline explicit itr_b(IT itr, const size_t s) noexcept:
		itr_(std::move(itr)), s_(s) { assert(s_>0); }


	/** assignment
	 */
	//! swap function
	inline friend void swap(itr_b& lhs, itr_b& rhs) noexcept {
		using namespace std;
		swap(lhs.itr_,rhs.itr_);
		swap(lhs.s_,rhs.s_);
	}


	/** @name information
	 */
	//! underlying vector iterator
	const IT& vitr() const noexcept { return itr_; }
	//! step
	ptrdiff_t s() const noexcept { return s_; }


	/** @name casts
	 */
	//! cast to difference type
	inline explicit operator diff() const noexcept { return (diff)vitr(); }
	//! cast to size_t
	inline explicit operator size_t() const noexcept { return (size_t)vitr(); }
	

	/** @name increment, decrement
	 */
	//! post increment
	inline RT& operator++() noexcept { 
		itr_+=s(); return *static_cast<RT*>(this);
	};
	//! pre increment
	inline RT operator++(int) noexcept {
		RT res=*static_cast<RT*>(this); ++(*this); return res;
	}
	//! post decrement
	inline RT& operator--() noexcept {
		itr_-=s(); return *static_cast<RT*>(this);
	}
	//! pre decrement
	inline RT operator--(int) noexcept {
		RT res=*static_cast<RT*>(this); --(*this); return res;
	}
	

	/** @name arithmetic operators
	 */
	//! iterator +=
	inline RT& operator+=(const diff rhs) noexcept {
		itr_+=(rhs*s()); return *static_cast<RT*>(this);
	}
	//! iterator +
	inline RT operator+(const diff rhs) const noexcept {
		RT res=*static_cast<const RT*>(this); res+=rhs; return res;
	}
	//! iterator -=
	inline RT& operator-=(const diff rhs) noexcept {
		itr_-=(rhs*s()); return *static_cast<RT*>(this); 
	}
	//! iterator -
	inline RT operator-(const diff rhs) const noexcept {
		RT res=*static_cast<const RT*>(this); res-=rhs; return res;
	}
	

	/** @name arithmetic as rhs
	 */
	//! itertator +
	inline friend RT operator+(const diff lhs, const RT& rhs) noexcept {
		return rhs+lhs;
	}
	

	/** @name const dereference
	 */
	//! const dereference operator
	inline ref operator*() const noexcept {
		return this->vitr().operator*();
	}
	//! const linear indexing operator
	inline val operator[](const diff i) const noexcept {
		return this->vitr().operator[](i * this->s());
	}
	//! const arrow operator
	inline ptr operator->() const noexcept {
		return this->vitr().operator->();
	}
	

	/** @name comparison operators
	 */
	//! comparison to iterator
	inline bool operator==(const RT& rhs) const noexcept { return vitr()==rhs.vitr(); }
	//! comparison to other iterator
	inline bool operator==(const RT_& rhs) const noexcept { return vitr()==rhs.vitr(); }
	//! comparison to iterator
	inline bool operator!=(const RT& rhs) const noexcept { return vitr()!=rhs.vitr(); }
	//! comparison to other iterator
	inline bool operator!=(const RT_& rhs) const noexcept { return vitr()!=rhs.vitr(); }
	//! comparison to iterator
	inline bool operator<(const RT& rhs) const noexcept { return vitr()<rhs.vitr(); }
	//! comparison to other iterator
	inline bool operator<(const RT_& rhs) const noexcept { return vitr()<rhs.vitr(); }
	//! comparison to iterator
	inline bool operator>(const RT& rhs) const noexcept { return vitr()>rhs.vitr(); }
	//! comparison to other iterator
	inline bool operator>(const RT_& rhs) const noexcept { return vitr()>rhs.vitr(); }
	//! comparison to iterator
	inline bool operator<=(const RT& rhs) const noexcept { return vitr()<=rhs.vitr(); }
	//! comparison to other iterator
	inline bool operator<=(const RT_& rhs) const noexcept { return vitr()<=rhs.vitr(); }
	//! comparison to iterator
	inline bool operator>=(const RT& rhs) const noexcept { return vitr()>=rhs.vitr(); }
	//! comparison to other iterator
	inline bool operator>=(const RT_& rhs) const noexcept { return vitr()>=rhs.vitr(); }
	

	/** @name difference
	 */
	//! iterator differnce
	inline friend diff operator-(const RT& lhs, const RT& rhs) noexcept {
		assert(lhs.s()==rhs.s());
		return ((diff)lhs-(diff)rhs)/lhs.s();
	}
	//! iterator differnce
	inline friend diff operator-(const RT& lhs, const RT_& rhs) noexcept {
		assert(lhs.s()==rhs.s());
		return ((diff)lhs-(diff)rhs)/lhs.s();
	}

	/** @name printing
	 */
	//! streaming operator
	friend std::ostream& operator<<(std::ostream& os, const itr_b& i) noexcept {
		return (os << "mi(" << i.vitr() << "," << (ptrdiff_t)i
			   << "," << i.s() << ")"); 
	}


protected:
	/** @name member variables
	 */
	IT itr_;	//!< underlying vector iterator
	diff s_;	//!< step
};

}
}


namespace ll__ {
	//! mesh nfo struct
	struct mesh_nfo {
		size_t M;			//!< number of elements in each columns
		std::vector<size_t> maj;	//!< data majority vector
		std::vector<size_t> D;		//!< number of grid points along each dimension
	};
}



/** mesh class, implements a mesh akin of the MATLAB commmand 'meshgrid'.
 * A mesh is implemented by a matrix underneath. Each point in the mesh
 * corresponds to a column in this matrix. The 'lines' across the mesh
 * can be laid out in different majority orders in the data, analogous
 * to row or column major in a matrix. Just like for lm_tMat, custom step
 * iterators are provided for traversing the mesh along 'lines'. \n
 * Each point in the mesh is a column in a matrix. Examples are: \n
 * - A kpoint mesh for bandstructures is a mesh where each column
 *   has 3 (usually) entries. The mesh may be spun across a different
 *   number of dimensions however, in case of 2D materials for example.
 * - The corresponding bandstructure. In this case each column in the
 *   matrix has Nb elements, one for each band, whereas the mesh itself
 *   is underwise identical to the kpoints mesh.
 */
template<class MT=lm__::fMat>
class ll_mesh {
public:
	/** @name constructors
	 */
	//! default constructor
	ll_mesh()=default;
	/** constructor from matrix, majority order and number of points
	 * @param base mesh matrix
	 * @param maj majority order in base
	 * @param D number of points in each dimension
	 */
	ll_mesh(MT base, std::vector<size_t> maj, std::vector<size_t> D) noexcept:
			base_(std::move(base)), s_(steps_(D,maj)),
			maj_(std::move(maj)), D_(std::move(D)) {
		assert(std::all_of(D_.cbegin(),D_.cend(),[](const size_t i)->bool{return i;}));
		assert(base_.N()==std::accumulate(D_.cbegin(),D_.cend(),size_t(1),
					std::multiplies<size_t>()));
	}
	/** constuctor from number of dimensions, majority oder and number of points
	 * @param base_dim number of dimensions in the matrix
	 * @param maj majority order in base
	 * @param D number of points in each dimension
	 */
	ll_mesh(const size_t base_dim, std::vector<size_t> maj, std::vector<size_t> D) noexcept:
		base_(base_dim,D.empty() ? 0: std::accumulate(D.cbegin(),D.cend(),
					size_t(1),std::multiplies<size_t>())),
					s_(steps_(D,maj)), maj_(std::move(maj)), D_(std::move(D)) {
		assert(std::all_of(D_.cbegin(),D_.cend(),[](const size_t i)->bool{return i;}));
	}
	/** constructor from mesh info struct
	 * @param S info struct
	 */
	ll_mesh(ll__::mesh_nfo S) noexcept: ll_mesh(S.M,std::move(S.maj),std::move(S.D)) {}
	

	/** @name information
	 */
	//! returns whether the mesh is empty
	inline bool empty() const noexcept { return base_.empty(); }
	//! returns the linear size of the mesh
	inline size_t size() const noexcept { return base_.N(); }
	//! returns the dimension of the mesh
	inline size_t meshDim() const noexcept { return maj().size(); }
	//! returns number of elements in each column of the matrix
	inline size_t M() const noexcept { return base_.M(); }
	//! returns the linear size of the mesh
	inline size_t N() const noexcept { return base_.N(); }
	/** hopping step for each dimension
	 * @param d dimension index
	 */
	inline size_t step(const size_t d) const noexcept {
		assert(d<meshDim());
		return s()[d];
	}
	/** majority layout for each dimension
	 * @param d dimension index
	 */
	inline size_t maj(const size_t d) const noexcept {
		assert(d<meshDim());
		return maj()[d];
	}
	/** number of points aling each dimension
	 * @param d dimension index
	 */
	inline size_t D(const size_t d) const noexcept {
		assert(d<meshDim());
		return D()[d];
	}
	/** position vector to linear index
	 * @param pos position vector
	 */
	inline size_t posToi(const std::vector<size_t>& pos) const noexcept {
		assert(pos.size()==meshDim());
		#ifndef NDEBUG
		for (size_t d=0; d!=pos.size(); ++d)
			assert(pos[d]<D(d));
		#endif

		size_t res=0;
		for (auto js=s().cbegin(),jse=s().cend(),jp=pos.cbegin(); js!=jse; ++js,++jp)
			res += (*jp)*(*js);
		return res;
	}
	/** linear index to position vector
	 * @param i linear index
	 */
	inline std::vector<size_t> iTopos(size_t i) const noexcept {
		assert(i<N());
		
		std::vector<size_t> res(meshDim());
		size_t dmn = std::accumulate(maj().cbegin(),maj().cend()-1,size_t(1),
				[this](const size_t s, const size_t i)->size_t
				{return s*D(i);});
		for (auto j=maj().crbegin(),e=maj().crend()-1; j!=e; ++j)
			res[*j] = i/dmn, i%=dmn, dmn/=D(*(j+1));
		res[maj().front()] = i/dmn;

		return res;
	}


	/** @name data access
	 */
	//! step vector const reference
	inline const std::vector<size_t>& s() const noexcept { return s_; }
	//! majority vector const reference
	inline const std::vector<size_t>& maj() const noexcept { return maj_; }
	//! number of points vector const reference
	inline const std::vector<size_t>& D() const noexcept { return D_; }


	/** @name printing
	 */
	/** write mesh to file
	 * @param fileName name of the file
	 * @param noheader skip header switch
	 */
	void writeToFile(const std::string& fileName, const bool noheader=false) const;


	/** @name public data member and cast to matrix base
	 */
	MT base_;	//!< underlying matrix
	//! cast to underlying matrix
	inline operator MT() const && { return MT(std::move(base_)); }


	/** @name direct mesh access
	 */
	//! direct column const access by linear index
	auto cAt(const size_t i) const noexcept -> const decltype(MT().cAt(0)) { return base_.cAt(i); }
	//! direct column access by linear index
	auto cAt(const size_t i) noexcept { return base_.cAt(i); }
	//! direct column const access by position vector
	auto cAt(const std::vector<size_t>& pos) const noexcept -> const decltype(MT().cAt(0)) {
		return cAt(posToi(pos));
	}
	//! direct column const access by position vector
	auto cAt(const std::vector<size_t>& pos) noexcept {
		return cAt(posToi(pos));
	}
	//! direct row const access
	auto rAt(const size_t i) const noexcept -> const decltype(MT().rAt(0)) { return base_.rAt(i); }
	//! direct row access
	auto rAt(const size_t i) noexcept { return base_.rAt(i); }


	/** @name regular matrix iterators
	 */
	//! column iterator to first column
	auto cBegin() noexcept { return base_.cBegin(); }
	//! past the end column iterator
	auto cEnd() noexcept { return base_.cEnd(); }
	//! column const_iterator to first column
	auto cBegin() const noexcept { return base_.cBegin(); }
	//! past the end column const_iterator
	auto cEnd() const noexcept { return base_.cEnd(); }
	//! column const_iterator to first column
	auto ccBegin() const noexcept { return base_.ccBegin(); }
	//! past the end column const_iterator
	auto ccEnd() const noexcept { return base_.ccEnd(); }
	//! row iterator to the first row
	auto rBegin() noexcept { return base_.rBegin(); }
	//! past the end row iterator
	auto rEnd() noexcept { return base_.rEnd(); }
	//! row const_iterator to the first row
	auto rBegin() const noexcept { return base_.rBegin(); }
	//! past the end row const_iterator
	auto rEnd() const noexcept { return base_.rEnd(); }
	//! row const_iterator to the first row
	auto crBegin() const noexcept { return base_.crBegin(); }
	//! past the end row const_iterator
	auto crEnd() const noexcept { return base_.crEnd(); }


private:
	/** @name helpers
	 */
	//! construct steps from maj
	inline static std::vector<size_t> steps_(const std::vector<size_t>& D,
						 const std::vector<size_t>& maj) noexcept {
		using namespace aux;
		assert(D.size()==maj.size());
		assert(std::all_of(D.cbegin(),D.cend(),[](const size_t i)->bool{return i;}));
		#ifndef NDEBUG
		std::vector<size_t> ck(maj.size());
		std::iota(ck.begin(),ck.end(),0);
		assert(std::is_permutation(ck.cbegin(),ck.cend(),maj.cbegin()));
		#endif

		size_t f=1;
		std::vector<size_t> res(maj.size());
		for (const auto i: maj)
			res[i] = f, f *= D[i];
		return res;
	}


	/** @name member variables
	 */
	std::vector<size_t> s_;		//!< step vector
	std::vector<size_t> maj_;	//!< majority order vector
	std::vector<size_t> D_;		//!< number of points vector

public:
	class itr;
	class c_itr;
	class r_itr;
	class cr_itr;

	/** @name  mesh iterator classes
	 */
	//! mesh iterator
	class itr: public ll__::mesh::itr_b<decltype(base_.cBegin()),itr,c_itr> {
	public:
		/** @name constructors
		 */
		using ll__::mesh::itr_b<decltype(base_.cBegin()),itr,c_itr>::itr_b;
		//! default constructor
		itr()=default;
		//! copy constructor
		inline itr(const itr& inp) noexcept: itr(inp.vitr(),inp.s()) {}
		//! constructor from reverse iterator
		inline itr(const r_itr& inp) noexcept: itr((inp-1).vitr()+1,inp.s()) {}


		/** @name assignment
		 */
		//! assignment from other iterator
		inline itr& operator=(const itr& inp) noexcept {
			this->itr_=inp.vitr(); this->s_=inp.s(); return *this;
		}
		//! assignment from reverse iterator
		inline c_itr& operator=(const r_itr& inp) noexcept {
			this->itr_=(inp-1).vitr()+1; this->s_=inp.s(); return *this;
		}
	};
	//! mesh const_iterator
	class c_itr: public ll__::mesh::itr_b<decltype(base_.ccBegin()),c_itr,itr> {
	public:
		/** @name constructors
		 */
		using ll__::mesh::itr_b<decltype(base_.ccBegin()),c_itr,itr>::itr_b;
		//! default constructor
		c_itr()=default;
		//! copy constructor
		inline c_itr(const c_itr& inp) noexcept: c_itr(inp.vitr(),inp.s()) {}
		//! constructor from iterator
		inline c_itr(const itr& inp) noexcept: c_itr(inp.vitr(),inp.s()) {}
		//! constructor from reverse const_iterator
		inline c_itr(const cr_itr& inp) noexcept: c_itr((inp-1).vitr()+1,inp.s()) {}
		//! constructor from reverse iterator
		inline c_itr(const r_itr& inp) noexcept: c_itr((inp-1).vitr()+1,inp.s()) {}


		/** @name assignment
		 */
		//! assignment from iterator
		inline c_itr& operator=(const itr& inp) noexcept {
			this->itr_=inp.vitr(); this->s_=inp.s(); return *this;
		}
		//! assignment from const_iterator
		inline c_itr& operator=(const c_itr& inp) noexcept {
			this->itr_=inp.vitr(); this->s_=inp.s(); return *this;
		}
		//! assignment from reverse iterator
		inline c_itr& operator=(const r_itr& inp) noexcept {
			this->itr_=(inp-1).vitr()+1; this->s_=inp.s(); return *this;
		}
		//! assignment from reverse const_iterator
		inline c_itr& operator=(const cr_itr& inp) noexcept {
			this->itr_=(inp-1).vitr()+1; this->s_=inp.s(); return *this;
		}
	};
	//! mesh reverse iterator
	class r_itr: public ll__::mesh::itr_b<decltype(base_.rcBegin()),r_itr,cr_itr> {
	public:
		/** @name constructors
		 */
		using ll__::mesh::itr_b<decltype(base_.rcBegin()),r_itr,cr_itr>::itr_b;
		//! default constructor
		r_itr()=default;
		//! constructor from iterator
		inline r_itr(const itr& inp) noexcept: r_itr((inp-1).vitr()+1,inp.s()) {}
		//! copy constructor
		inline r_itr(const r_itr& inp) noexcept: r_itr(inp.vitr(),inp.s()) {}


		/** @name assignment
		 */
		//! assignment from iterator
		inline itr& operator=(const itr& inp) noexcept {
			this->itr_=(inp-1).vitr()+1; this->s_=inp.s(); return *this;
		}
		//! assignment from reverse iterator
		inline c_itr& operator=(const r_itr& inp) noexcept {
			this->itr_=inp.vitr(); this->s_=inp.s(); return *this;
		}
	

		/** @name base
		 */
		//! base iterator
		inline itr base() const noexcept { return itr(*this); }
	};
	//! mesh reverse const_iterator
	class cr_itr: public ll__::mesh::itr_b<decltype(base_.crcBegin()),cr_itr,r_itr> {
	public:
		/** @name constructors
		 */
		using ll__::mesh::itr_b<decltype(base_.crcBegin()),cr_itr,r_itr>::itr_b;
		//! default constructor
		cr_itr()=default;
		//! constructor from const_iterator
		inline cr_itr(const c_itr& inp) noexcept: cr_itr((inp-1).vitr()+1,inp.s()) {}
		//! constructor from iterator
		inline cr_itr(const itr& inp) noexcept: cr_itr((inp-1).vitr()+1,inp.s()) {}
		//! copy constructor
		inline cr_itr(const cr_itr& inp) noexcept: cr_itr(inp.vitr(),inp.s()) {}
		//! constructor from reverse iterator
		inline cr_itr(const r_itr& inp) noexcept: cr_itr(inp.vitr(),inp.s()) {}


		/** @name assignment
		 */
		//! assignment from iterator
		inline c_itr& operator=(const itr& inp) noexcept {
			this->itr_=(inp-1).vitr()+1; this->s_=inp.s(); return *this;
		}
		//! assignment from const_iterator
		inline c_itr& operator=(const c_itr& inp) noexcept {
			this->itr_=(inp-1).vitr()+1; this->s_=inp.s(); return *this;
		}
		//! assignment from reverse iterator
		inline c_itr& operator=(const r_itr& inp) noexcept {
			this->itr_=inp.vitr(); this->s_=inp.s(); return *this;
		}
		//! assignment from reverse const_iterator
		inline c_itr& operator=(const cr_itr& inp) noexcept {
			this->itr_=inp.vitr(); this->s_=inp.s(); return *this;
		}
		

		/** @name base
		 */
		//! base iterator
		inline c_itr base() const noexcept { return c_itr(*this); }
	};


public:
	/** @name mesh iterators
	 */
	/** mesh iterator to a column at the beginning of a 'line'
	 * @param d dimension index
	 * @param pos position vector to a grid point in the 'line'
	 */
	inline itr begin(const size_t d, const std::vector<size_t>& pos) noexcept {
		return itr(base_.cBegin()+posToi(pos),step(d))-pos[d];
	}
	/** mesh iterator to a column past the end of a 'line'
	 * @param d dimension index
	 * @param pos position vector to a grid point in the 'line'
	 */
	inline itr end(const size_t d, const std::vector<size_t>& pos) noexcept {
		assert(pos.size()==meshDim());
		return begin(d,pos)+D(d);
	}
	/** mesh const_iterator to a column at the beginning of a 'line'
	 * @param d dimension index
	 * @param pos position vector to a grid point in the 'line'
	 */
	inline c_itr begin(const size_t d, const std::vector<size_t>& pos) const noexcept {
		return c_itr(base_.cBegin()+posToi(pos),step(d))-pos[d];
	}
	/** mesh const_iterator to a column past the end of a 'line'
	 * @param d dimension index
	 * @param pos position vector to a grid point in the 'line'
	 */
	inline c_itr end(const size_t d, const std::vector<size_t>& pos) const noexcept {
		assert(pos.size()==meshDim());
		return begin(d,pos)+D(d);
	}
	/** mesh const_iterator to a column at the beginning of a 'line'
	 * @param d dimension index
	 * @param pos position vector to a grid point in the 'line'
	 */
	inline c_itr cbegin(const size_t d, const std::vector<size_t>& pos) const noexcept {
		return begin(d,pos);
	}
	/** mesh const_iterator to a column past the end of a 'line'
	 * @param d dimension index
	 * @param pos position vector to a grid point in the 'line'
	 */
	inline c_itr cend(const size_t d, const std::vector<size_t>& pos) const noexcept {
		return end(d,pos);
	}
	/** reverse mesh iterator to a column at the end of a 'line'
	 * @param d dimension index
	 * @param pos position vector to a grid point in the 'line'
	 */
	inline r_itr rbegin(const size_t d, const std::vector<size_t>& pos) noexcept {
		return r_itr(end(d,pos));
	}
	/** reverse mesh iterator to a column past the beginning of a 'line'
	 * @param d dimension index
	 * @param pos position vector to a grid point in the 'line'
	 */
	inline r_itr rend(const size_t d, const std::vector<size_t>& pos) noexcept {
		return r_itr(begin(d,pos));
	}
	/** reverse mesh const_iterator to a column at the end of a 'line'
	 * @param d dimension index
	 * @param pos position vector to a grid point in the 'line'
	 */
	inline cr_itr rbegin(const size_t d, const std::vector<size_t>& pos) const noexcept {
		return cr_itr(end(d,pos));
	}
	/** reverse mesh const_iterator to a column past the beginning of a 'line'
	 * @param d dimension index
	 * @param pos position vector to a grid point in the 'line'
	 */
	inline cr_itr rend(const size_t d, const std::vector<size_t>& pos) const noexcept {
		return cr_itr(begin(d,pos));
	}
	/** reverse mesh const_iterator to a column at the end of a 'line'
	 * @param d dimension index
	 * @param pos position vector to a grid point in the 'line'
	 */
	inline cr_itr crbegin(const size_t d, const std::vector<size_t>& pos) const noexcept {
		return rbegin(d,pos);
	}
	/** reverse mesh const_iterator to a column past the beginning of a 'line'
	 * @param d dimension index
	 * @param pos position vector to a grid point in the 'line'
	 */
	inline cr_itr crend(const size_t d, const std::vector<size_t>& pos) const noexcept {
		return rend(d,pos);
	}
};

//! streaming operator
template<class MT>
inline std::ostream& operator<<(std::ostream& os, const ll_mesh<MT>& inp) noexcept {
	using namespace aux;
	return (os << "size: " << inp.size() << "\n"
		   << "D: " << inp.D() << "\n"
		   << "maj: " << inp.maj() << "\n"
		   << "step: " << inp.s() << "\n"
		   << "cpx: " << inp.base_.cpx());
}


namespace ll__ {

	//! get nfo from mesh
	template<class MT>
	ll__::mesh_nfo size(const ll_mesh<MT>& inp) noexcept {
		return {inp.M(),inp.maj(),inp.D()};
	}


	/* generic mesh generation
	 */
	/** generate mesh from bounds, majority order and number of points
	 * @param bounds upper and lower bounds, size DIMx2
	 * @param maj majority order vector
	 * @param D number of points vector
	 */
	ll_mesh<> genMesh(const fMat& bounds, std::vector<size_t> maj, std::vector<size_t> D) noexcept;
	/** generate mesh from upper and lower bound, majority order and number of points
	 * @param lb lower bound for all dimensions
	 * @param ub upper bound for all dimensions
	 * @param D number of points for all dimensions
	 * @param maj majority order vector
	 * @param r restriction vector
	 */
	ll_mesh<> genMesh(const double lb, const double ub, const size_t D,
			std::vector<size_t> maj, const rv& r) noexcept;
	

	/* majority styles
	 */
   	//! default majority {0,1,...}
	inline std::vector<size_t> maj_default(const size_t d) noexcept {
		std::vector<size_t> res(d);
		std::iota(res.begin(),res.end(),size_t(0));
		return res;
	}
	//! MATLAB style majority, i.e. {1,0},{1,0,2},...
	inline std::vector<size_t> maj_MATLAB(const rv& r) noexcept {
		auto res = maj_default(r.size());
		const auto I = ninds(r);
		if (I.size()>1) std::swap(res[I[0]],res[I[1]]);
		return res;
	}
	

	/* mesh generators
	 */
	/** mesh in a cell with given density
	 * @param rho point density
	 * @param B basis
	 * @param lb lower bound
	 * @param ub upper bound
	 * @param maj majority order vector
	 * @param r restriction vector
	 */
	ll_mesh<> genMesh_cell(const double rho, const fMat& B,
		const double lb, const double ub, std::vector<size_t> maj, const rv& r);
	/** mesh in a cell with given density
	 * @param rho point density
	 * @param B basis
	 * @param lb lower bound
	 * @param ub upper bound
	 * @param r restriction vector
	 */
	inline ll_mesh<> genMesh_cell_MATLAB(const double rho, const fMat& B, const double lb,
		const double ub, const rv& r) { return genMesh_cell(rho,B,lb,ub,maj_MATLAB(r),r); }


	/** mesh on an integer grid
	 * @param bounds upper and lower bounds, size DIMx2
	 * @param maj majority order vector
	 */
	ll_mesh<> genMesh_int(const fMat& bounds, std::vector<size_t> maj) noexcept;
	/** mesh on an integer grid, forward majority {0,1,...}
	 * @param bounds upper and lower bounds, size DIMx2
	 */
	inline ll_mesh<> genMesh_int_fw(const fMat& bounds) noexcept {
		std::vector<size_t> maj(bounds.M());
		std::iota(maj.begin(),maj.end(),0);
		return genMesh_int(bounds,std::move(maj));
	}
	/** mesh on an integer grid, backward majority {N,N-1,...}
	 * @param bounds upper and lower bounds, size DIMx2
	 */
	inline ll_mesh<> genMesh_int_bw(const fMat& bounds) noexcept {
		std::vector<size_t> maj(bounds.M());
		std::iota(maj.rbegin(),maj.rend(),0);
		return genMesh_int(bounds,std::move(maj));
	}
}

#endif // _LL_MESH_

/** @}
 */
