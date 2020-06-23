// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_IO_
#define _LL_IO_

#include "ll_defs.h"
#include "ll_compound.h"
#include "lm_fn.h"
#include "ll_cell.h"
#include "ll_hbonds.h"
#include "ll_hbondss.h"
#include "ll_mesh.h"
#include <string>
#include <vector>
#include <cfloat>
#include <fstream>


namespace ll__ {
	
	/* dft related
	 */
	/** read VASP POSCAR file.
	 * @param fileName name of the file
	 * @param d dimension of space
	 * @return ll__::psc struct
	 */
	psc readPOSCAR(const std::string& fileName, const size_t d=DIM__);
	/** print VASP POSCAR file
	 * @param fileName name of the file
	 * @param B basis of the unit cell
	 * @param Ap atomic positions
	 * @param T atomic type vector
	 * @param id id strings vector
	 * @param a0 lattice constant
	 * @param direct switch to print in direct(true) or cartesian(false) coodrinates
	 * @param prec floating point precision
	 * @param header first line in the POSCAR file
	 */
	void printPOSCAR(const std::string& fileName, const fMat& B, const fMat& Ap,
			const aTv& T, const idv& id={},
			const double a0=1.0, const bool direct=true,
			const size_t prec=PPREC__, const std::string& header="");
	/** print VASP POSCAR file
	 * @param fileName name of the file
	 * @param B basis of the unit cell
	 * @param Ap atomic positions
	 * @param id id strings vector
	 * @param a0 lattice constant
	 * @param direct switch to print in direct(true) or cartesian(false) coodrinates
	 * @param prec floating point precision
	 * @param header first line in the POSCAR file
	 */
	void printPOSCAR(const std::string& fileName, const fMat& B, const fMat& Ap,
			const idv& id,
			const double a0=1.0, const bool direct=true,
			const size_t prec=PPREC__, const std::string& header="");
	/** print VASP POSCAR file
	 * @param fileName name of the file
	 * @param inp the unit cell
	 * @param direct switch to print in direct(true) or cartesian(false) coodrinates
	 * @param strip switch to strip indices
	 * @param prec floating point precision
	 * @param header first line in the POSCAR file
	 */
	void printPOSCAR(const std::string& fileName, const ll_cell& inp,
			const bool direct=true, const bool strip=false, 
			const size_t prec=PPREC__, const std::string& header="");
	/** read Fermi energy from VASP OUTCAR file
	 * @param fileName name of the file
	 */
	double readEf(const std::string& fileName=OUTCAR__);


	/* omen related
	 */
	/** generate OMEN style restriction vector
	 * @param d dimension of space
	 * @return restriction vector, 111 for 1D, 110 for 2D, 100 for 3D
	 */
	rv genr(const size_t d) noexcept;
	//! read OMEN layer matrix
	Ap_T readLayerMatrix(const std::string& fileName=LM__);
	/** print OMEN material file
	 * @param fileName name of the file
	 * @param E band edges
	 * @param Norb number of orbitals vectors
	 * @param mass mass vector
	 */
	void printOmf(const std::string& fileName,const vb_cb& E,
		      const std::vector<size_t>& Norb, const std::vector<double>& mass);
	/** print OMEN material file
	 * @param fileName name of the file
	 * @param D ll__::omf OMEN material struct
	 */
	inline void printOmf(const std::string& fileName, const omf& D) {
		printOmf(fileName,D.E,D.Norb,D.mass);
	}
	/** print OMEN lattice file
	 * @param fileName name of the file
	 * @param B basis
	 * @param Ap atomic positions
	 * @param id id strings
	 * @param nn number of next neighbors
	 * @param bl bond length
	 */
	void printOlf(const std::string& fileName, const fMat& B, const fMat& Ap,
		      const idv& id, const size_t nn, const double bl);
	/** print OMEN lattice file
	 * @param fileName name of the file
	 * @param D ll__::olf OMEN lattice struct
	 */
	inline void printOlf(const std::string& fileName, const olf& D) {
		printOlf(fileName,D.B,D.Ap,D.id,D.nn,D.bl);
	}
	/** read OMEN material file
	 * @param fileName name of the file
	 * @return ll__::omf OMEN material struct
	 */
	omf readOmf(const std::string& fileName);
	/** read OMEN lattice file
	 * @param fileName name of the file
	 * @return ll__::olf OMEN lattice struct
	 */
	olf readOlf(const std::string& fileName);
	

	/* wannier90 related
	 */
	//! read dimension restrictions from a wannier90 hr file
	rv hrDim(const std::string& fileName=HR__);
	//! read basis from a wannier90 wout file
	fMat readB(const std::string& fileName=WOUT__);
	/** read atomic positions and id strings from wannier90 wout file.
	 * @param fileName name of the wout file
	 * @param direct switch to read in direct(true) or cartesian(false) coordinates
	 * @param abl atomic blacklist
	 * @return ll__::Ap_id struct
	 */
	Ap_id readAp(const std::string& fileName=WOUT__, const bool direct=true,
		    		const std::vector<size_t>& abl={});
	/** read Wannier centers and spreads from a wannier90 wout file
	 * @param fileName name of the file
	 * @param wbl Wannier blacklist
	 * @return ll__::Wp_s struct
	 */
	Wp_s readWp(const std::string& fileName=WOUT__, const aTv& wbl={});
	/** read operator in Wannier representation
	 * @param strm stream to read from
	 * @param frmt format string
	 * @param Nw number of Wannier functions
	 * @param NR number of R vectors
	 * @param wbl Wannier blacklist
	 * @param keep lambda determining which matrices to keep
	 */
	R_H<> readOp(std::istream& strm, const char* frmt, const size_t Nw, const size_t NR,
			const aTv& wbl, const std::function<bool(const cMat& inp)>& keep);
	/** read Hamiltonian in Wannier representation
	 * @param fileName name of the file
	 * @param wbl Wannier blacklist
	 * @param keep lambda determining which matrices to keep
	 */
	R_H<> readHr(const std::string& fileName=HR__, const aTv& wbl={},
		const std::function<bool(const cMat& inp)>& keep=
			[](const cMat& inp){return true;});
	/** read Hamiltonian in Wannier representation
	 * @param fileName name of the file
	 * @param wbl Wannier blacklist
	 * @param tol tolerance level, matrices with at least one entry above this will be included
	 */
	inline R_H<> readHr(const std::string& fileName, const aTv& wbl, const double tol) {
		return readHr(fileName,wbl,[tol](const cMat& inp)->bool
			{return std::any_of(inp.cbegin(),inp.cend(),
				[tol](const auto i)->bool{return cmph(i,tol);});});
	}
	/** write wannier90 Hamiltonian to file
	 * @param fileName name of the file
	 * @param hr ll__::R_H struct holding the Hamiltonian data and R vectors
	 * @param header header to put into the file
	 */
	template<class MT> void writeHr(const std::string& fileName, const R_H<MT>& hr,
			const std::string& header="");
	//! read wannier90 AMN file
	std::vector<cMat> readAMN(const std::string& fileName);
	/** read wannier position operator
	 * @param fileName name of the file
	 * @param d dimension of space
	 * @param wbl Wannier blacklist
	 * @param keep lambda determining which matrices to keep
	 */
	R_H<> readXYZ(const std::string& fileName=R__, const size_t d=0, const aTv& wbl={},
		const std::function<bool(const cMat& inp)>& keep=[](const cMat& inp){return true;});
	/** read x component of wannier position operator
	 * @param fileName name of the file
	 * @param wbl Wannier blacklist
	 * @param keep lambda determining which matrices to keep
	 */
	inline R_H<> readX(const std::string& fileName=R__, const aTv& wbl={},
		const std::function<bool(const cMat& inp)>& keep=[](const cMat& inp){return true;})
			{ return readXYZ(fileName,0,wbl,keep); }
	/** read y component of wannier position operator
	 * @param fileName name of the file
	 * @param wbl Wannier blacklist
	 * @param keep lambda determining which matrices to keep
	 */
	inline R_H<> readY(const std::string& fileName=R__, const aTv& wbl={},
		const std::function<bool(const cMat& inp)>& keep=[](const cMat& inp){return true;})
			{ return readXYZ(fileName,1,wbl,keep); }
	/** read z component of wannier position operator
	 * @param fileName name of the file
	 * @param wbl Wannier blacklist
	 * @param keep lambda determining which matrices to keep
	 */
	inline R_H<> readZ(const std::string& fileName=R__, const aTv& wbl={},
		const std::function<bool(const cMat& inp)>& keep=[](const cMat& inp){return true;})
			{ return readXYZ(fileName,2,wbl,keep); }
	//! read eigenvalues from wannier90 eig file
	fMat readEig(const std::string& filename=WEIG__);
	//! read Wannier transformation, only k and U
	k_U readWannierTransf(const std::string& fileName);
	/** read Wannier transformation, k, U and Udis
	 * @param fileNameU name of the U matrices file
	 * @param fileNameUdis name of the Udis matrices file
	 */
	k_U readWannierTransf(const std::string& fileNameU, const std::string& fileNameUdis);
	//! read Wannier transformation from wannier90 CHK file
	k_U readChk(const std::string& fileName=CHK__);
	/** print Wannier matching to file
	 * @param fileName name of the file
	 * @param I Wannier matching indices
	 */
	void printWiToFile(const std::string& fileName, const wi& I);


	/** write gnu plot files
	 * @param xdat data on the x axis
	 * @param ydat data on the y axis
	 * @param seed seed to use for output files, files are seed.gnu and seed.dat
	 * @param prefix prefix to put in front of filenames
	 * @param window plot window
	 * @param xticks ticks on the x axis
	 * @param xlabels labels for the ticks on the x axis
	 * @param yticks ticks on the y axis
	 * @param ylabels labels for the ticks on the y axis
	 */
	void writeGnuPlot(const fMat& xdat, const fMat& ydat,
		const std::string& seed, const std::string& prefix="./", std::vector<double> window={},
		const std::vector<double>& xticks={}, const idv& xlabels={},
		const std::vector<double>& yticks={}, const idv& ylabels={});
	

	struct struct_report;
	/** print structural information for a wbh (ll_hbonds or ll_hbondss)
	 * @param rep report as generated by ll__::analyzeStructure
	 * @param mask bitmask, ...0001 prints exact matches, ...0010 prints approximate matches
	 * @param W the wbh
	 * @param os stream to print into
	 */
	template<class WT>
	std::ostream& printStructureReport(const struct_report& rep, const size_t mask,
			const WT& W, std::ostream& os) noexcept;


	//! extract unit cell from a wbh
	ll_cell extractCell(const std::string& fileName);
}

#endif// _LL_IO_

/** @}
 */
