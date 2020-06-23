// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @defgroup liblat liblat
 * This library provides an infrastructure for dealing with unit cells,
 * bonds and Wannier interactions. It also provides a range of
 * functions to read and write to files provided by the dft code VASP,
 * the wannierization tool Wannier90 and the quantum transport tool OMEN.
 * The main functionality is upscaling Hamiltonians expressed in terms
 * of a small unit cell to larger structures, i.e. generate Hamiltonians
 * for quantum transport simulations. Functionality for manipulating
 * atomic lattices and compute bandstructures is provided as well.
 * All of the above makes heavy use of the libmat library.
 */

/** @namespace ll__
 * The main namespace for liblat
 */

/** @namespace aux
 * Namespace for auxiliary functionality
 */
