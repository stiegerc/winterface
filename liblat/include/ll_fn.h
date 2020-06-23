// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_FN_
#define _LL_FN_

#include "ll_defs.h"
#include "ll_compound.h"
#include "ll_cell.h"
#include "ll_lambda.h"
#include "ll_hbonds.h"
#include "ll_mesh.h"
#include <unordered_map>


struct ll_BStest_input;

namespace ll__ {

	//! next neighbour matrix
	extern const fMat NN;

	/* basis finding
	 */
	/** function to find a basis expansion given a basis and a template
	 * @param B the basis
	 * @param V the template
	 * @param lim the upper limit for expansion coefficients
	 * @param tol tolerance level as when a floating point number can be considered an integer
	 */
	fMat findBasisExpansion(const fMat& B, const fMat& V, const double lim=1000.0, const double tol=1e-2);
	/** function to find a basis given a basis and a template
	 * @param B the basis
	 * @param V the template
	 * @param lim the upper limit for expansion coefficients
	 * @param tol tolerance level as when a floating point number can be considered an integer
	 */
	inline fMat findBasis(const fMat& B, const fMat& V, const double lim=1000.0, const double tol=1e-2) {
		return B.prod(findBasisExpansion(B,V,lim,tol));
	}
	/** function to 'orthogonalize' a basis from a restriction vector. This function will try to make
	 * basis vectors deemed restricted orthogonal to the others
	 * @param B basis
	 * @param r restriction vector
	 */
	fMat orthogonalize(fMat B, const rv& r) noexcept;
	

	/* linear path generation
	 */
	/** generate points along a path including position information
	 * @param P points matrix
	 * @param Nps vector holding the number of points along each path between two points in P
	 */
	p_p genPath(const fMat& P, std::vector<size_t> Nps);
	/** generate points along a path including position information. The number of points
	 * for each segment is determined automatically by the relative distance of the points in P.
	 * @param P points matrix
	 * @param Np total number of points
	 * @param B optional basis, if P is in direct coordinates
	 */
	p_p genPath(fMat P, const size_t Np, const fMat& B={});


	/* metrics and clustering
	 */
	/** generate a grid of integers akin to the meshgrid function in MATLAB
	 * @param bounds upper and lower bounds, size DIMx2
	 * @param maj majority order in the grid, size DIM
	 */
	fMat genGrid(const fMat& bounds, std::vector<size_t> maj) noexcept;
	/** generate a grid of integers akin to the meshgrid function in MATLAB
	 * @param bounds upper and lower bounds, size DIMx2
	 */
	fMat genGrid(const fMat& bounds) noexcept;
	/** generate a grid of integers -1:1:1 akin to the meshgrid function in MATLAB
	 * @param r restriction vector, size DIM
	 * @param maj majority order in the grid DIM
	 */
	fMat genNNmat(const rv& r, std::vector<size_t> maj) noexcept;
	/** generate a grid of integers -1:1:1 akin to the meshgrid function in MATLAB
	 * @param r restriction vector, size DIM
	 */
	fMat genNNmat(const rv& r) noexcept;
	/** generate a grid of integers -1:1:1 akin to the meshgrid function in MATLAB.
	 * This version assumes no spacial restrictions.
	 * @param d dimension of space
	 */
	inline fMat genNNmat(const size_t d=DIM__) noexcept { return genNNmat(rv(d,false)); }
	/** function to find the distance of two points in periodic space
	 * @param B basis
	 * @param p1 the first point in direct coordinates
	 * @param p2 the first point in direct coordinates
	 * @param NN next neighbor matrix
	 */
	double dist(const fMat& B, const fArray& p1, const fArray& p2, fMat NN) noexcept;
	/** function to find the distance of two points in periodic space, tolerant mode.
	 * This version will keep increasing the result from the minimum until a large enough
	 * gap is found. Useful for slightly warped lattices.
	 * @param B basis
	 * @param p1 the first point in direct coordinates
	 * @param p2 the first point in direct coordinates
	 * @param f tolerance factor, f<0: upper limit is |f| times shortest bond length
	 *                            f>0: upper limit is f
	 * @param NN next neighbor matrix
	 */
	double dist(const fMat& B, const fArray& p1, const fArray& p2, const double f, fMat NN) noexcept;	
	/** function to find the distance of two points in periodic space
	 * @param B basis
	 * @param p1 the first point in direct coordinates
	 * @param p2 the first point in direct coordinates
	 * @param NN next neighbor matrix
	 * @return ll__::d_b struct including the bond length and the actual bonds
	 * 		indicating which image is the closest
	 */
	d_b distb(const fMat& B, const fArray& p1, const fArray& p2, fMat NN) noexcept;
	/** function to find the distance of two points in periodic space, tolerant mode.
	 * This version will keep increasing the result from the minimum until a large enough
	 * gap is found. Useful for slightly warped lattices.
	 * @param B basis
	 * @param p1 the first point in direct coordinates
	 * @param p2 the first point in direct coordinates
	 * @param f tolerance factor, f<0: upper limit is |f| times shortest bond length
	 *                            f>0: upper limit is f
	 * @param NN next neighbor matrix
	 * @return ll__::d_b struct including the bond length and the actual bonds
	 * 		indicating which image is the closest
	 */
	d_b distb(const fMat& B, const fArray& p1, const fArray& p2, const double f, fMat NN) noexcept;
	/** function to return a distance matrix of points with themselves.
	 * @param B basis
	 * @param P positions in direct coordinates
	 * @param NN next neighbor matrix
	 */
	fMat genDmat(const fMat& B, const fMat& P, const fMat& NN) noexcept;
	/** function to return a distance matrix of one set of points with another.
	 * @param B basis
	 * @param P1 first set of points in direct coordinates
	 * @param P2 second set of points in direct coordinates
	 * @param NN next neighbor matrix
	 */
	fMat genDmat(const fMat& B, const fMat& P1, const fMat& P2, const fMat& NN) noexcept;
	/** function to return a distance matrix of points with themselves.
	 * This version will keep increasing the distance from the minimum until a large enough
	 * gap is found. Useful for slightly warped lattices.
	 * @param B basis
	 * @param P positions in direct coordinates
	 * @param f tolerance factor, f<0: upper limit is |f| times shortest bond length
	 *                            f>0: upper limit is f
	 * @param NN next neighbor matrix
	 */
	fMat genDmat(const fMat& B, const fMat& P, const double f, const fMat& NN) noexcept;
	/** function to return a distance matrix of one set of points with another.
	 * This version will keep increasing the distance from the minimum until a large enough
	 * gap is found. Useful for slightly warped lattices.
	 * @param B basis
	 * @param P1 first set of points in direct coordinates
	 * @param P2 second set of points in direct coordinates
	 * @param f tolerance factor, f<0: upper limit is |f| times shortest bond length
	 *                            f>0: upper limit is f
	 * @param NN next neighbor matrix
	 */
	fMat genDmat(const fMat& B, const fMat& P1, const fMat& P2, const double f, const fMat& NN) noexcept;
	/** function to find the center of mass in periodic space
	 * @param inp array holding a set of points
	 * @param w optional weights
	 */
	fMat com(const fArray& inp, const fMat& w={}) noexcept;
	/** function to generate a grouping of points into clusters using the dbscan algorithm.
	 * @param D points to investigate
	 * @param minpts minimum number of points to count as a cluster
	 * @param eps small parameters
	 */
	wi dbscan(const fMat& D, const size_t minpts, const double eps);
	

	/* generate Wannier hamiltonian
	 */
	/** generate a Hamiltonian in Wannier representation.
	 * @param inp kpoints and U matrices
	 * @param E eigenvalues
	 * @param R R vectors
	 */
	R_H<> genHam(const k_U& inp, const fMat& E, const fMat& R) noexcept;
	/** generate a Hamiltonian in Wannier representation
	 * @param B basis
	 * @param Ap atomic postions
	 * @param id vector of id strings
	 * @param W a wbh (ll_hbonds or ll_hbondss)
	 * @param R R vectors
	 * @param strict switch to allow only exact matched(true) or also approximate matches(false)
	 */
	template<class MT, class WT> R_H<MT> genHam(const fMat& B, const fMat& Ap,
			const idv& id, const WT& W, fMat R, const bool strict) noexcept;
	/** generate a Hamiltonian in Wannier representation
	 * @param B basis
	 * @param Ap atomic postions
	 * @param id vector of id strings
	 * @param r restriction vector
	 * @param W a wbh (ll_hbonds or ll_hbondss)
	 * @param strict switch to allow only exact matched(true) or also approximate matches(false)
	 */
	template<class MT, class WT> R_H<MT> genHam(const fMat& B, const fMat& Ap,
			const idv& id, const rv& r, const WT& W, const bool strict) noexcept;
	/** generate a Hamiltonian in Wannier representation
	 * @param cell unit cell
	 * @param W a wbh (ll_hbonds or ll_hbondss)
	 * @param R R vectors
	 * @param strict switch to allow only exact matched(true) or also approximate matches(false)
	 */
	template<class MT, class WT> R_H<MT> inline genHam(const ll_cell& cell,
					const WT& W, fMat R, const bool strict) noexcept {
		return genHam<MT>(cell.B(),cell.getcAp(),cell.id(cell.type()),W,R,strict);
	}
	/** generate a Hamiltonian in Wannier representation
	 * @param cell unit cell
	 * @param r restriction vector
	 * @param W a wbh (ll_hbonds or ll_hbondss)
	 * @param strict switch to allow only exact matched(true) or also approximate matches(false)
	 */
	template<class MT, class WT> R_H<MT> inline genHam(const ll_cell& cell, const rv& r,
					const WT& W, const bool strict) noexcept {
		return genHam<MT>(cell.B(),cell.getcAp(),cell.id(cell.type()),r,W,strict);
	}
	/** generate a Hamiltonian in Wannier representation using the cell in the wbh
	 * @param W a wbh (ll_hbonds or ll_hbondss)
	 * @param R R vectors
	 */
	template<class MT, class WT> R_H<MT> inline genHam(const WT& W, fMat R) noexcept
			{ return genHam<MT>(W.cell(),W,R,false); }
	/** generate a Hamiltonian in Wannier representation using the cell in the wbh
	 * @param W a wbh (ll_hbonds or ll_hbondss)
	 */
	template<class MT, class WT> R_H<MT> inline genHam(const WT& W) noexcept
			{ return genHam<MT>(W.cell(),W.r(),W,false); }


	/* Wannier matching
	 */
	/** perform 'wannier matching' by matching Wannier centers to 'atomic' positions
	 * @param B basis
	 * @param Ap 'atomic' positions in direct coordinates
	 * @param Wp Wannier centers
	 * @param NN next neighbor matrix
	 */
	wi matchCenters(const fMat& B, const fMat& Ap, const fMat& Wp, const fMat& NN) noexcept;
	/** perform 'wannier matching' by matching Wannier centers to 'atomic' positions
	 * @param cell unit cell
	 * @param Wp Wannier centers
	 * @param NN next neighbor matrix
	 */
	inline wi matchCenters(const ll_cell& cell, const fMat& Wp, const fMat& NN) noexcept {
		return matchCenters(cell.B(),cell.Ap(),Wp,NN);
	}
	/** function to 'balance' out a Wannier matching. This can be useful for suboptimal
	 * Wannierizations where some Wannier centers can be a bit further away than usual. They
	 * may then get matched to the 'wrong' atom. This function attempts to balance things out
	 * by shifting the least obvious matches such that the same atomic types all have the same
	 * number of Wannier centers matched to them despite the spacial warping of the initial distribution.
	 * @param I 'naive' Wannier matching
	 * @param cell unit cell
	 * @param Wp Wannier centers
	 * @param NN next neighbor matrix
	 */
	void balance(wi& I, const ll_cell& cell, const fMat& Wp, const fMat& NN) noexcept;
	/** function to 'clusterize' Wannier centers. This attempts to find clusters of Wannier centers
	 * and matches them to each other.
	 * @param B basis
	 * @param Wp Wannier centers
	 * @param eps small parameter, paramter for dbscan
	 * @param minpts minimum number of points to count as a cluster, parameter for dbscan
	 * @param NN next neighbor matrix
	 */
	wi clusterize(const fMat& B, const fMat& Wp, const double eps, const size_t minpts, const fMat& NN) noexcept;
	/** function to find 'atomic' positions for 'clusterized' Wannier centers. This uses the 'center
	 * of mass' as the 'atomic' positition with the Wannier spreads as weights.
	 * @param B basis
	 * @param Wp Wannier centers
	 * @param s Wannier spreads
	 * @param I Wannier matching
	 */
	fMat genCenters(const fMat& B, const fMat& Wp, const fMat& s, const wi& I) noexcept;


	/* Wannier tools
	 */
	//! convert Wannier matching to atomic type vector
	aTv wiToT(const wi& I) noexcept;
	//! convert atomic type vector to Wannier matching
	wi TToWi(const aTv& T) noexcept;
	//! check integrity of Wannier matching
	bool checkWi(const wi& I) noexcept;
	/** check integrity of Wannier matching
	 * @param I Wannier matching
	 * @param N number of Wannier functions in I
	 */
	bool checkWi(const wi& I, const size_t N) noexcept;
	//! function to find the largest entries of real and imaginary parts in a Wannier Hamiltonian
	fMat checkHr(const R_H<>& hr) noexcept;


	/* simple bandstructure calculation
	 */
	/** function to find the valence and conduction band edges
	 * @param Ef Fermi energy
	 * @param E band energies
	 */
	vb_cb findBandEdges(const double Ef, const fMat& E);
	/** function to compute the bandstructure on a trace.
	 * @param hr Hamiltonian in Wannier representation
	 * @param k kpoints
	 * @param Nthreads number of threads to use in parallel section
	 */
	template<class MT>
	fMat calcBS(const R_H<MT>& hr, const fMat& k, const size_t Nthreads=NTHREADS__) noexcept;
	/** function to compute the bandstructure on a mesh.
	 * @param hr Hamiltonian in Wannier representation
	 * @param k mesh of kpoints
	 * @param Nthreads number of threads to use in parallel section
	 */
	template<class MT>
	inline ll_mesh<> calcBS(const R_H<MT>& hr, const ll_mesh<>& k, const size_t Nthreads=NTHREADS__) noexcept {
		return ll_mesh<>(calcBS(hr,k.base_,Nthreads),k.maj(),k.D());
	}
	/** function to compute the bandstructure on a trace for a super cell using a folding approach.
	 * @param hr Hamiltonian in Wannier representation
	 * @param k kpoints
	 * @param B basis of the super cell
	 * @param Bp basis of the primitive cell
	 * @param Nthreads number of threads to use in parallel section
	 */
	template<class MT>
	fMat calcFoldedBS(const R_H<MT>& hr, const fMat& k, const fMat& B, const fMat& Bp,
				const size_t Nthreads=NTHREADS__) noexcept;
	/** function to compute the bandstructure on a mesh for a super cell using a folding approach.
	 * @param hr Hamiltonian in Wannier representation
	 * @param k mesh of kpoints
	 * @param B basis of the super cell
	 * @param Bp basis of the primitive cell
	 * @param Nthreads number of threads to use in parallel section
	 */
	template<class MT>
	inline ll_mesh<> calcFoldedBS(const R_H<MT>& hr, const ll_mesh<>& k, const fMat& B, const fMat& Bp,
				const size_t Nthreads=NTHREADS__) noexcept {
		return ll_mesh<>(calcFoldedBS(hr,k.base_,B,Bp,Nthreads),k.maj(),k.D());
	}


	/* bandstructure calculation including gradients and curvature tensors
	 */
	//! energies, gradients and curvature tensors struct
	template<class MT>
	struct egc {
		MT E;		//!< energy
		ll_mesh<> G;	//!< gradient
		ll_mesh<> C;	//!< curvature
		//! Number of bands
		inline size_t Nb() const noexcept { return E.M(); }
		//! number of kpoints
		inline size_t Nk() const noexcept { return E.N(); }
		//! dimension of space
		inline size_t dim() const noexcept { return G.empty() ? 0: G.D(0); }
	};
	/** function to compute the bandstructure including gradient and curvature on a trace.
	 * @param hr Hamiltonian in Wannier representation
	 * @param k kpoints
	 * @param B basis (real space)
	 * @param Nthreads number of threads to use in parallel section
	 * @return ll__::egc struct holding energy, gradient and curvature
	 */
	template<class MT>
	egc<fMat> calcBS_gc(const R_H<MT>& hr, const fMat& k, const fMat& B,
			const size_t Nthreads=NTHREADS__) noexcept;
	/** function to compute the bandstructure including gradient and curvature on a mesh.
	 * @param hr Hamiltonian in Wannier representation
	 * @param k mesh of kpoints
	 * @param B basis (real space)
	 * @param Nthreads number of threads to use in parallel section
	 * @return ll__::egc struct holding energy, gradient and curvature
	 */
	template<class MT>
	inline egc<ll_mesh<>> calcBS_gc(const R_H<MT>& hr, const ll_mesh<>& k, const fMat& B,
			const size_t Nthreads=NTHREADS__) noexcept {
		const auto buff = calcBS_gc(hr,k.base_,B,Nthreads);
		return {ll_mesh<>(std::move(buff.E),k.maj(),k.D()),
			ll_mesh<>(std::move(buff.G.base_),cat(0,k.maj()+1),cat(k.M(),k.D())),
			ll_mesh<>(std::move(buff.C.base_),cat(0,1,k.maj()+2),cat(k.M(),k.M(),k.D()))
		};
	}
	/** function to compute the bandstructure including gradient and curvature on a trace in
	 * a super cell using a folding approach.
	 * @param hr Hamiltonian in Wannier representation
	 * @param k kpoints
	 * @param B basis of the super cell (real space)
	 * @param Bp basis of the primtive cell (real space)
	 * @param Nthreads number of threads to use in parallel section
	 * @return ll__::egc struct holding energy, gradient and curvature
	 */
	template<class MT>
	egc<fMat> calcFoldedBS_gc(const R_H<MT>& hr, const fMat& k, const fMat& B, const fMat& Bp,
			const size_t Nthreads=NTHREADS__) noexcept;
	/** function to compute the bandstructure including gradient and curvature on a mesh in
	 * a super cell using a folding approach.
	 * @param hr Hamiltonian in Wannier representation
	 * @param k mesh of kpoints
	 * @param B basis of the super cell (real space)
	 * @param Bp basis of the primtive cell (real space)
	 * @param Nthreads number of threads to use in parallel section
	 * @return ll__::egc struct holding energy, gradient and curvature
	 */
	template<class MT>
	inline egc<ll_mesh<>> calcFoldedBS_gc(const R_H<MT>& hr, const ll_mesh<>& k,
			const fMat& B, const fMat& Bp, const size_t Nthreads=NTHREADS__) noexcept {
		const auto buff = calcFoldedBS_gc(hr,k.base_,B,Bp,Nthreads);
		return {ll_mesh<>(std::move(buff.E),k.maj(),k.D()),
			ll_mesh<>(std::move(buff.G.base_),cat(0,k.maj()+1),cat(k.M(),k.D())),
			ll_mesh<>(std::move(buff.C.base_),cat(0,1,k.maj()+2),cat(k.M(),k.M(),k.D()))
		};
	}


	/* analyze structure
	 */
	//! structure report struct
	struct struct_report {
		//! union for exact matches
		union iiN {
			struct {
				size_t i1;		//!< starting index
				size_t i2;		//!< ending index
				size_t N = 0;		//!< number of occurrence
			};
			std::array<size_t,3> dat;	//!< array access
			//! streaming operator
			friend inline std::ostream& operator<<(std::ostream& os, const iiN& inp) noexcept
				{ return (os<<"("<<inp.i1<<":"<<inp.i2<<":"<<inp.N<<")"); }
		};
		//! union for approximate matches
		union iijN {
			struct {
				size_t i1;		//!< starting index
				size_t i2;		//!< ending index
				size_t j;		//!< real target
				size_t N = 0;		//!< number of occurrence
			};
			std::array<size_t,4> dat;	//!< array access
			//! streaming operator
			friend inline std::ostream& operator<<(std::ostream& os, const iijN& inp) noexcept
				{ return (os<<"("<<inp.i1<<":"<<inp.i2<<":"<<inp.j<<":"<<inp.N<<")"); }
		};
		std::vector<iiN> exact;		//!< vector for exact matches
		std::vector<iijN> approx;	//!< vector for approximate matches
		std::vector<size_t> unmatched;	//!< vector for failed matches
	};
	/** function to analyze a structure in tandem with a corresponding wbh. What are the exact matches
	 * during queries for bonds? What are the approximate matches? etc.
	 * @param Ap atomic positions
	 * @param id id strings
	 * @param W wbh (ll_hbonds or ll_hbondss)
	 * @param tol tolerance level as to which matches should be included. Matches with interaction
	 * 		matrices with no entry exceeding this are discarded.
	 * @param f factor for the cutoff radius, Rcut = W.radius()*f
	 */
	template<class WT>
	struct_report analyzeStructure(const fMat& Ap, const idv& id,
				       const WT& W, const double tol, const double f=.5) noexcept;
}

#endif // _LL_FN_

/** @}
 */
