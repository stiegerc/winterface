// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @defgroup libmat libmat
 * This library provides real and complex matrices using the FORTRAN
 * libraries BLAS/LAPACK underneath for speed. The main feature is
 * the use of row and column iterators allowing for convenient use of
 * the C++ standard library. Due to the column major ordering of the
 * memory, matrices should be viewed as vectors of columns.
 */

/** @namespace lm__
 * The main namespace for libmat
 */
