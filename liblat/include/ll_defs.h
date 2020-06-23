// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _LL_DEFS_
#define _LL_DEFS_

#include "lm_defs.h"


/* general
 */
#ifndef DIM__
#define DIM__ 3				//!< dimension of space
#endif
#ifndef TXTPREC__
#define TXTPREC__ 6			//!< default precision for printing floats to textfiles
#endif


/* wannier bonds related
 */
#ifndef WBH__
#define WBH__ "wannier90.wbh"		//!< wbh default name
#endif


/* dft related
 */
#ifndef POSCAR__
#define POSCAR__ "POSCAR"               //!< Vasp POSCAR default name
#endif
#ifndef OUTCAR__
#define OUTCAR__ "OUTCAR"               //!< default filename for VASP OUTCAR file
#endif
#ifndef QEI__
#define QEI__ "scf.in"                  //!< Quantum Espresso input file default name
#endif
#ifndef QEO__
#define QEO__ "scf.out"                 //!< default filename for Quantum Espresso out file
#endif


/* omen related
 */
#ifndef LM__
#define LM__ "Layer_Matrix.dat"         //!< default filename for OMEN layer matrix
#endif
#ifndef OLF__
#define OLF__ "lattice_dat"             //!< OMEN lattice file default name
#endif
#ifndef OMF__
#define OMF__ "ph_mat_par"              //!< OMEN material file default name
#endif
#ifndef EPC__
#define EPC__ 1000000000		//!< entries per column in OMEN material file
#endif
#ifndef OIS__
#define OIS__ "stump.cmd"		//!< default filename for OMEN input stump
#endif


/* wannier90 related
 */
#ifndef WOUT__
#define WOUT__ "wannier90.wout"		//!< default filename for wannier90 wout file
#endif
#ifndef HR__
#define HR__ "wannier90_hr.dat"		//!< default filename for wannier90 hamiltonian data
#endif
#ifndef R__
#define R__ "wannier90_r.dat"		//!< default filename for wannier90 position operator data
#endif
#ifndef CHK__
#define CHK__ "wannier90.chk.fmt"	//!< default filename for wannier90 chk data
#endif
#ifndef WEIG__
#define WEIG__ "wannier90.eig"		//!< default filename for wannier90 eig
#endif
#ifndef UMAT__
#define UMAT__ "wannier90_u.mat"	//!< default filename for wannier90 u mats
#endif
#ifndef WTOL__
#define WTOL__ 1e-4			//!< floating point tolerance to used with wannier data
#endif					//!< dictated by wannier90 output precision


/* default threads
 */
#ifndef NTHREADS__
#define NTHREADS__ 4			//!< default number of threads in omp sections
#endif


#endif // _LL_DEFS_

/** @}
 */
