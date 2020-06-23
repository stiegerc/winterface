// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup libmat
 * @{
 */

#ifndef _LM_DEFS_
#define _LM_DEFS_

#include <complex>
#include <string>

// default setting
#ifndef SINGLE__
#ifndef DOUBLE__
#define DOUBLE__			//!< compile for double precision
#endif
#endif

// precision
#ifdef SINGLE__
#define RE__ float			//!< real type for single precision
#define CPX__ std::complex<float>	//!< complex type for single precision
#define MTOL__ 1e-2f			//!< default tolerance for single precision
#endif
#ifdef DOUBLE__
#define RE__ double			//!< real type for double precision
#define CPX__ std::complex<double>	//!< complex type for double precision
#define MTOL__ 1e-8			//!< default tolerance for double precision
#endif

// off limits constant
#ifndef NPOS__
#define NPOS__ std::string::npos	//!< npos value
#endif

// printing precision
#ifndef PPREC__
#define PPREC__ 12			//!< default floating point precision
#endif

#endif // _LM_DEFS_

/** @}
 */
