// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#ifndef _BLAS_
#define _BLAS_

#ifndef FORTRAN_NAME
    #if defined FORTRAN_ADD_
        #define FORTRAN_NAME(X,Y) (X ## _)
    #elif defined FORTRAN_NOCHANGE
        #define FORTRAN_NAME(X,Y) (X)
    #elif defined FORTRAN_UPCASE
        #define FORTRAN_NAME(X,Y) (Y)
    #else
        #define FORTRAN_NAME(X,Y) (X ## _)
    #endif
#endif

#include <complex>
typedef std::complex<float> scpx;
typedef std::complex<double> dcpx;


/* FORTRAN FUNCTIONS */
#ifdef __cplusplus
extern "C" {
#endif



// Blas level 1
// DOT
//float FORTRAN_NAME(sdot,SDOT)(const int* n, const float* sx, const int* incx, const float* sy, const int* incy);
float FORTRAN_NAME(sdot,SDOT)(const int* n, const float* sx, const int* incx, const float* sy, const int* incy);
void FORTRAN_NAME(cdotc,CDOTC)(scpx* res, const int* n, const scpx* cx, const int* incx,
				const scpx* cy, const int* incy);
void FORTRAN_NAME(cdotu,CDOTU)(scpx* res, const int* n, const scpx* cx, const int* incx,
				const scpx* cy, const int* incy);
double FORTRAN_NAME(ddot,DDOT)(const int* n, const double* dx, const int* incx, const double* dy, const int* incy);
void FORTRAN_NAME(zdotc,ZDOTC)(dcpx* res, const int* n, const dcpx* zx, const int* incx,
				const dcpx* zy, const int* incy);
void FORTRAN_NAME(zdotu,ZDOTU)(dcpx* res, const int* n, const dcpx* zx, const int* incx,
				const dcpx* zy, const int* incy);
void FORTRAN_NAME(saxpy,SAXPY)(const int* n, const float* sa, const float* sx, const int* incx,
				const float* sy, const int* incy);
void FORTRAN_NAME(caxpy,CAXPY)(const int* n, const scpx* ca, const scpx* cx, const int* incx,
				const scpx* cy, const int* incy);
void FORTRAN_NAME(daxpy,DAXPY)(const int* n, const double* da, const double* dx, const int* incx,
				const double* dy, const int* incy);
void FORTRAN_NAME(zaxpy,ZAXPY)(const int* n, const dcpx* za, const dcpx* zx, const int* incx,
				const dcpx* zy, const int* incy);

// Blas level 2
// GEMV
void FORTRAN_NAME(sgemv,SGEMV)(const char* trans, const int* m, const int* n, const float* alpha, const float* a,
                               const int* lda, const float* x, const int* incx, const float* beta, float* y,
                               const int* incy);
void FORTRAN_NAME(cgemv,CGEMV)(const char* trans, const int* m, const int* n, const scpx* alpha, const scpx* a,
                               const int* lda, const scpx* x, const int* incx, const scpx* beta, scpx* y,
                               const int* incy);
void FORTRAN_NAME(dgemv,DGEMV)(const char* trans, const int* m, const int* n, const double* alpha, const double* a,
                               const int* lda, const double* x, const int* incx, const double* beta, double* y,
                               const int* incy);	
void FORTRAN_NAME(zgemv,ZGEMV)(const char* trans, const int* m, const int* n, const dcpx* alpha, const dcpx* a,
                               const int* lda, const dcpx* x, const int* incx, const dcpx* beta, dcpx* y,
                               const int* incy);					   

// Blas level 3
// GEMM
void FORTRAN_NAME(sgemm,SGEMM)(const char* transa, const char* transb, const int* m, const int* n, const int* k,
                               const float* alpha, const float* a, const int* lda, const float* b, const int* ldb,
                               const float* beta, float *c, const int* ldc);							   
void FORTRAN_NAME(cgemm,CGEMM)(const char* transa, const char* transb, const int* m, const int* n, const int* k,
			       const scpx* alpha, const scpx* a, const int* lda, const scpx* b, const int* ldb,
			       const scpx* beta, scpx* c, const int* ldc);
void FORTRAN_NAME(dgemm,DGEMM)(const char* transa, const char* transb, const int* m, const int* n, const int* k,
                               const double* alpha, const double* a, const int* lda, const double* b, const int* ldb,
                               const double* beta, double* c, const int* ldc);							   
void FORTRAN_NAME(zgemm,ZGEMM)(const char* transa, const char* transb, const int* m, const int* n, const int* k,
			       const dcpx* alpha, const dcpx* a, const int* lda, const dcpx* b, const int* ldb,
			       const dcpx* beta, dcpx* c, const int* ldc);

// Lapack
// GETRF
void FORTRAN_NAME(sgetrf,SGETRF)(const int* m, const int* n, float* a, const int* lda, int* ipiv, int* info);
void FORTRAN_NAME(cgetrf,CGETRF)(const int* m, const int* n, scpx* a, const int* lda, int* ipiv, int* info);
void FORTRAN_NAME(dgetrf,DGETRF)(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
void FORTRAN_NAME(zgetrf,ZGETRF)(const int* m, const int* n, dcpx* a, const int* lda, int* ipiv, int* info);


// GETRI
void FORTRAN_NAME(sgetri,SGETRI)(const int* n, float* a, const int* lda, const int* ipiv, float* work,
                                const int* lwork, int* info);
void FORTRAN_NAME(cgetri,CGETRI)(const int *n, scpx* a, const int *lda, const int *ipiv, scpx* work,
                                const int* lwork, int* info);
void FORTRAN_NAME(dgetri,DGETRI)(const int *n, double* a, const int *lda, const int *ipiv, double* work,
                                const int* lwork, int* info);
void FORTRAN_NAME(zgetri,ZGETRI)(const int *n, dcpx* a, const int *lda, const int *ipiv, dcpx* work,
                                const int* lwork, int* info);	

// GEQRF
void FORTRAN_NAME(sgeqrf,SGEQRF)(const int* m, const int* n, float* a, const int* lda, float* tau, float* work,
				const int* lwork, int* info);
void FORTRAN_NAME(cgeqrf,CGEQRF)(const int* m, const int* n, scpx* a, const int* lda, scpx* tau, scpx* work,
				const int* lwork, int* info);
void FORTRAN_NAME(dgeqrf,DGEQRF)(const int* m, const int* n, double* a, const int* lda, double* tau, double* work,
				const int* lwork, int* info);
void FORTRAN_NAME(zgeqrf,ZGEQRF)(const int* m, const int* n, dcpx* a, const int* lda, dcpx* tau, dcpx* work,
				const int* lwork, int* info);

// ORGQR,UNGQR
void FORTRAN_NAME(sorgqr,SORGQR)(const int* m, const int* n, const int* k, float* a, const int* lda, float* tau,
				float* work, const int* lwork, int* info);
void FORTRAN_NAME(cungqr,CUNGQR)(const int* m, const int* n, const int* k, scpx* a, const int* lda, scpx* tau,
				scpx* work, const int* lwork, int* info);
void FORTRAN_NAME(dorgqr,DORGQR)(const int* m, const int* n, const int* k, double* a, const int* lda, double* tau,
				double* work, const int* lwork, int* info);
void FORTRAN_NAME(zungqr,ZUNGQR)(const int* m, const int* n, const int* k, dcpx* a, const int* lda, dcpx* tau,
				dcpx* work, const int* lwork, int* info);

// ORMQR,UNMQR
void FORTRAN_NAME(sormqr,SORMQR)(const char* side, const char* trans, const int* m, const int* n, const int* k,
				float* a, const int* lda, float* tau, float* c, const int* ldc, float* work,
				const int* lwork, int* info);
void FORTRAN_NAME(cunmqr,CUNMQR)(const char* side, const char* trans, const int* m, const int* n, const int* k,
				scpx* a, const int* lda, scpx* tau, scpx* c, const int* ldc, scpx* work,
				const int* lwork, int* info);
void FORTRAN_NAME(dormqr,DORMQR)(const char* side, const char* trans, const int* m, const int* n, const int* k,
				double* a, const int* lda, double* tau, double* c, const int* ldc, double* work,
				const int* lwork, int* info);
void FORTRAN_NAME(zunmqr,ZUNMQR)(const char* side, const char* trans, const int* m, const int* n, const int* k,
				dcpx* a, const int* lda, dcpx* tau, dcpx* c, const int* ldc, dcpx* work,
				const int* lwork, int* info);

// POTRF
void FORTRAN_NAME(spotrf,SPOTRF)(const char* uplo, const int* n, float* a, const int* lda, int* info);
void FORTRAN_NAME(cpotrf,CPOTRF)(const char* uplo, const int* n, scpx* a, const int* lda, int* info);
void FORTRAN_NAME(dpotrf,DPOTRF)(const char* uplo, const int* n, double* a, const int* lda, int* info);
void FORTRAN_NAME(zpotrf,ZPOTRF)(const char* uplo, const int* n, dcpx* a, const int* lda, int* info);

// POTRS
void FORTRAN_NAME(spotrs,SPOTRS)(const char* uplo, const int* n, const int* nrhs, const float* a, const int* lda,
                                 float* b, const int* ldb, int* info);
void FORTRAN_NAME(cpotrs,CPOTRS)(const char* uplo, const int* n, const int* nrhs, const scpx* a, const int* lda,
                                 scpx* b, const int* ldb, int* info);
void FORTRAN_NAME(dpotrs,DPOTRS)(const char* uplo, const int* n, const int* nrhs, const double* a, const int* lda,
                                 double* b, const int* ldb, int* info);
void FORTRAN_NAME(zpotrs,ZPOTRS)(const char* uplo, const int* n, const int* nrhs, const dcpx* a, const int* lda,
                                 dcpx* b, const int* ldb, int* info);

// GESV
void FORTRAN_NAME(sgesv,SGESV)(const int* n, const int* nrhs, float* a, const int* lda, int* ipiv,
                               float* b, const int* ldb, int* info);
void FORTRAN_NAME(cgesv,CGESV)(const int* n, const int* nrhs, scpx* a, const int* lda, int* ipiv,
                               scpx* b, const int* ldb, int* info);
void FORTRAN_NAME(dgesv,DGESV)(const int* n, const int* nrhs, double* a, const int* lda, int* ipiv,
                               double* b, const int* ldb, int* info);
void FORTRAN_NAME(zgesv,ZGESV)(const int* n, const int* nrhs, dcpx* a, const int* lda, int* ipiv,
                               dcpx* b, const int* ldb, int* info);

// GESVD																					
void FORTRAN_NAME(sgesvd,SGESVD)(const char* jobu, const char* jobvt, const int* m, const int* n, float* a,
                                 const int* lda, float* s, float* u, const int* ldu, float* vt, const int* ldvt,
                                 float* work, const int* lwork, int* info);
void FORTRAN_NAME(cgesvd,CGESVD)(const char* jobu, const char* jobvt, const int* m, const int* n, scpx* a,
                                 const int* lda, float* s, scpx* u, const int* ldu, scpx* vt, const int* ldvt,
                                 scpx* work, const int* lwork, float* rwork, int* info);																					
void FORTRAN_NAME(dgesvd,DGESVD)(const char* jobu, const char* jobvt, const int* m, const int* n, double* a,
                                 const int* lda, double* s, double* u, const int* ldu, double* vt, const int* ldvt,
                                 double* work, const int* lwork, int* info);
void FORTRAN_NAME(zgesvd,ZGESVD)(const char* jobu, const char* jobvt, const int* m, const int* n, dcpx* a,
                                 const int* lda, double* s, dcpx* u, const int* ldu, dcpx* vt, const int* ldvt,
                                 dcpx* work, const int* lwork, double* rwork, int* info);

// ILAENV
int FORTRAN_NAME(ilaenv,ILAENV)(const int* ispec, const char* name, const char* opts, const int* n1, const int* n2,
                                const int* n3, const int* n4);
								
// GEEV
void FORTRAN_NAME(sgeev,SGEEV)(const char* jobvl, const char* jobvr, const int* n, float* a, const int* lda,
				float* wr, float* wi, float* vl, const int* ldvl, float* vr,
				const int* ldvr, float* work, const int* lwork, int* info);
void FORTRAN_NAME(cgeev,CGEEV)(const char* jobvl, const char* jobvr, const int* n, scpx* a, const int* lda,
				scpx* wr, scpx* vl, const int* ldvl, scpx* vr, const int* ldvr,
				scpx* work, const int* lwork, float* rwork, int* info);
void FORTRAN_NAME(dgeev,DGEEV)(const char* jobvl, const char* jobvr, const int* n, double* a, const int* lda,
				double* wr, double* wi, double* vl, const int* ldvl, double* vr,
				const int* ldvr, double* work, const int* lwork, int* info);
void FORTRAN_NAME(zgeev,ZGEEV)(const char* jobvl, const char* jobvr, const int* n, dcpx* a, const int* lda,
				dcpx* wr, dcpx* vl, const int* ldvl, dcpx* vr, const int* ldvr,
				dcpx* work, const int* lwork, double* rwork, int* info);
									
// SYEV
void FORTRAN_NAME(ssyev,SSYEV)(const char* jobz, const char* uplo, const int* n, float* a, const int* lda,
				float* w, float* work, const int* lwork, int* info);
void FORTRAN_NAME(dsyev,DSYEV)(const char* jobz, const char* uplo, const int* n, double* a, const int* lda,
				double* w, double* work, const int* lwork, int* info);
					
// HEEV
void FORTRAN_NAME(cheev,CHEEV)(const char* jobz, const char* uplo, const int* n, scpx* a, const int* lda,
				float* w, scpx* work, const int* lwork, float* rwork, int* info);
void FORTRAN_NAME(zheev,ZHEEV)(const char* jobz, const char* uplo, const int* n, dcpx* a, const int* lda,
				double* w, dcpx* work, const int* lwork, double* rwork, int* info);		

#ifdef __cplusplus
}
#endif





// C INTERFACE

// Blas level 1
// DOT
inline float c_xdot(const int n, const float* sx, const int incx, const float* sy, const int incy) {
	return FORTRAN_NAME(sdot,SDOT)(&n,sx,&incx,sy,&incy);
}
inline scpx c_xdotc(const int n, const scpx* cx, const int incx, const scpx* cy, const int incy) {
	
	scpx res{.0,.0};
	for (auto cxe=cx+n*incx; cx!=cxe; cx+=incx,cy+=incy)
		res += std::conj(*cx)*(*cy);
	return res;
	
//	doesn't work with the netlib version for some reason
//	scpx res;
//	FORTRAN_NAME(cdotc,CDOTC)(&res,&n,cx,&incx,cy,&incy);
//	return res;
}
inline scpx c_xdotu(const int n, const scpx* cx, const int incx, const scpx* cy, const int incy) {
	
	scpx res{.0,.0};
	for (auto cxe=cx+n*incx; cx!=cxe; cx+=incx,cy+=incy)
		res += (*cx)*(*cy);
	return res;
	
//	doesn't work with the netlib version for some reason
//	scpx res;
//	FORTRAN_NAME(cdotu,CDOTU)(&res,&n,cx,&incx,cy,&incy);
//	return res;
}
inline double c_xdot(const int n, const double* dx, const int incx, const double* dy, const int incy) {
	return FORTRAN_NAME(ddot,DDOT)(&n,dx,&incx,dy,&incy);
}
inline dcpx c_xdotc(const int n, const dcpx* zx, const int incx, const dcpx* zy, const int incy) {

	dcpx res{.0,.0};
	for (auto zxe=zx+n*incx; zx!=zxe; zx+=incx,zy+=incy)
		res += std::conj(*zx)*(*zy);
	return res;
	
//	doesn't work with the netlib version for some reason
//	dcpx res;
//	FORTRAN_NAME(zdotc,ZDOTC)(&res,&n,zx,&incx,zy,&incy);
//	return res;
}
inline dcpx c_xdotu(const int n, const dcpx* zx, const int incx, const dcpx* zy, const int incy) {
	
	dcpx res{.0,.0};
	for (auto zxe=zx+n*incx; zx!=zxe; zx+=incx,zy+=incy)
		res += (*zx)*(*zy);
	return res;

//	doesn't work with the netlib version for some reason
//	dcpx res;
//	FORTRAN_NAME(zdotu,ZDOTU)(&res,&n,zx,&incx,zy,&incy);
//	return res;
}
// AXPY
inline void c_xaxpy(const int n, const float sa, const float* sx, const int incx, const float* sy, const int incy) {
	FORTRAN_NAME(saxpy,SAXPY)(&n,&sa,sx,&incx,sy,&incy);
}
inline void c_xaxpy(const int n, const scpx& ca, const scpx* cx, const int incx, const scpx* cy, const int incy) {
	FORTRAN_NAME(caxpy,CAXPY)(&n,&ca,cx,&incx,cy,&incy);
}
inline void c_xaxpy(const int n, const double da, const double* dx, const int incx, const double* dy, const int incy) {
	FORTRAN_NAME(daxpy,DAXPY)(&n,&da,dx,&incx,dy,&incy);
}
inline void c_xaxpy(const int n, const dcpx& za, const dcpx* zx, const int incx, const dcpx* zy, const int incy) {
	FORTRAN_NAME(zaxpy,ZAXPY)(&n,&za,zx,&incx,zy,&incy);
}

// Blas level 2
// GEMV
inline void c_xgemv(const char trans, const int m, const int n, const float alpha, const float* a, const int lda,
                    const float* x, const int incx, const float beta, float* y, const int incy) {
	FORTRAN_NAME(sgemv,SGEMV)(&trans,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}
inline void c_xgemv(const char trans, const int m, const int n, const scpx alpha, const scpx* a, const int lda,
                    const scpx* x, const int incx, const scpx beta, scpx* y, const int incy) {
	FORTRAN_NAME(cgemv,CGEMV)(&trans,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}
inline void c_xgemv(const char trans, const int m, const int n, const double alpha, const double* a, const int lda,
                    const double* x, const int incx, const double beta, double* y, const int incy) {
	FORTRAN_NAME(dgemv,DGEMV)(&trans,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}
inline void c_xgemv(const char trans, const int m, const int n, const dcpx alpha, const dcpx* a, const int lda,
                    const dcpx* x, const int incx, const dcpx beta, dcpx* y, const int incy) {
	FORTRAN_NAME(zgemv,ZGEMV)(&trans,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}

// Blas level 3
// GEMM
inline void c_xgemm(const char transa, const char transb, const int m, const int n, const int k, const float alpha,
                    const float* a, const int lda, const float* b, const int ldb, const float beta, float* c,
                    const int ldc) {
	FORTRAN_NAME(sgemm,SGEMM)(&transa,&transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
}
inline void c_xgemm(const char transa, const char transb, const int m, const int n, const int k, const scpx alpha,
                    const scpx* a, const int lda, const scpx* b, const int ldb, const scpx beta, scpx* c,
                    const int ldc) {
	FORTRAN_NAME(cgemm,CGEMM)(&transa,&transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
}
inline void c_xgemm(const char transa, const char transb, const int m, const int n, const int k, const double alpha,
                    const double* a, const int lda, const double* b, const int ldb, const double beta, double* c,
                    const int ldc) {
	FORTRAN_NAME(dgemm,DGEMM)(&transa,&transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
}
inline void c_xgemm(const char transa, const char transb, const int m, const int n, const int k, const dcpx alpha,
                    const dcpx* a, const int lda, const dcpx* b, const int ldb, const dcpx beta, dcpx* c,
                    const int ldc) {
	FORTRAN_NAME(zgemm,ZGEMM)(&transa,&transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
}

// Lapack
// GETRF
inline void c_xgetrf(const int m, const int n, float* a, const int lda, int* ipiv, int* info) {
	FORTRAN_NAME(sgetrf,SGETRF)(&m,&n,a,&lda,ipiv,info);
}
inline void c_xgetrf(const int m, const int n, scpx* a, const int lda, int* ipiv, int* info) {
	FORTRAN_NAME(cgetrf,CGETRF)(&m,&n,a,&lda,ipiv,info);
}
inline void c_xgetrf(const int m, const int n, double* a, const int lda, int* ipiv, int* info) {
	FORTRAN_NAME(dgetrf,DGETRF)(&m,&n,a,&lda,ipiv,info);
}
inline void c_xgetrf(const int m, const int n, dcpx* a, const int lda, int* ipiv, int* info) {
	FORTRAN_NAME(zgetrf,ZGETRF)(&m,&n,a,&lda,ipiv,info);
}

// GETRI
inline void c_xgetri(const int n, float* a, const int lda, const int* ipiv, float* work,
                     const int lwork, int* info) {
	FORTRAN_NAME(sgetri,SGETRI)(&n,a,&lda,ipiv,work,&lwork,info);
}
inline void c_xgetri(const int n, scpx* a, const int lda, const int* ipiv, scpx* work,
                     const int lwork, int* info) {
	FORTRAN_NAME(cgetri,CGETRI)(&n,a,&lda,ipiv,work,&lwork,info);
}
inline void c_xgetri(const int n, double* a, const int lda, const int* ipiv, double* work,
                     const int lwork, int* info) {
	FORTRAN_NAME(dgetri,DGETRI)(&n,a,&lda,ipiv,work,&lwork,info);
}
inline void c_xgetri(const int n, dcpx* a, const int lda, const int* ipiv, dcpx* work,
                     const int lwork, int* info) {
	FORTRAN_NAME(zgetri,ZGETRI)(&n,a,&lda,ipiv,work,&lwork,info);
}

// GEQRF
inline void c_xgeqrf(const int m, const int n, float* a, const int lda, float* tau, float* work,
		     const int lwork, int* info) {
	FORTRAN_NAME(sgeqrf,SGEQRF)(&m,&n,a,&lda,tau,work,&lwork,info);
}
inline void c_xgeqrf(const int m, const int n, scpx* a, const int lda, scpx* tau, scpx* work,
		     const int lwork, int* info) {
	FORTRAN_NAME(cgeqrf,CGEQRF)(&m,&n,a,&lda,tau,work,&lwork,info);
}
inline void c_xgeqrf(const int m, const int n, double* a, const int lda, double* tau, double* work,
		     const int lwork, int* info) {
	FORTRAN_NAME(dgeqrf,DGEQRF)(&m,&n,a,&lda,tau,work,&lwork,info);
}
inline void c_xgeqrf(const int m, const int n, dcpx* a, const int lda, dcpx* tau, dcpx* work,
		     const int lwork, int* info) {
	FORTRAN_NAME(zgeqrf,ZGEQRF)(&m,&n,a,&lda,tau,work,&lwork,info);
}

// ORGQR,UNGQR
inline void c_xungqr(const int m, const int n, const int k, float* a, const int lda, float* tau,
		     float* work, const int lwork, int* info) {
	FORTRAN_NAME(sorgqr,SORGQR)(&m,&n,&k,a,&lda,tau,work,&lwork,info);
}
inline void c_xungqr(const int m, const int n, const int k, scpx* a, const int lda, scpx* tau,
		     scpx* work, const int lwork, int* info) {
	FORTRAN_NAME(cungqr,CUNGQR)(&m,&n,&k,a,&lda,tau,work,&lwork,info);
}
inline void c_xungqr(const int m, const int n, const int k, double* a, const int lda, double* tau,
		     double* work, const int lwork, int* info) {
	FORTRAN_NAME(dorgqr,DORGQR)(&m,&n,&k,a,&lda,tau,work,&lwork,info);
}
inline void c_xungqr(const int m, const int n, const int k, dcpx* a, const int lda, dcpx* tau,
		     dcpx* work, const int lwork, int* info) {
	FORTRAN_NAME(zungqr,ZUNGQR)(&m,&n,&k,a,&lda,tau,work,&lwork,info);
}

// ORMQR,UNMQR
inline void c_xunmqr(const char side, const char trans, const int m, const int n, const int k,
		     float* a, const int lda, float* tau, float* c, const int ldc, float* work,
		     const int lwork, int* info) {
	FORTRAN_NAME(sormqr,SORMQR)(&side,&trans,&m,&n,&k,a,&lda,tau,c,&ldc,work,&lwork,info);
}
inline void c_xunmqr(const char side, const char trans, const int m, const int n, const int k,
		     scpx* a, const int lda, scpx* tau, scpx* c, const int ldc, scpx* work,
		     const int lwork, int* info) {
	FORTRAN_NAME(cunmqr,CUNMQR)(&side,&trans,&m,&n,&k,a,&lda,tau,c,&ldc,work,&lwork,info);
}
inline void c_xunmqr(const char side, const char trans, const int m, const int n, const int k,
		     double* a, const int lda, double* tau, double* c, const int ldc, double* work,
		     const int lwork, int* info) {
	FORTRAN_NAME(dormqr,DORMQR)(&side,&trans,&m,&n,&k,a,&lda,tau,c,&ldc,work,&lwork,info);
}
inline void c_xunmqr(const char side, const char trans, const int m, const int n, const int k,
		     dcpx* a, const int lda, dcpx* tau, dcpx* c, const int ldc, dcpx* work,
		     const int lwork, int* info) {
	FORTRAN_NAME(zunmqr,ZUNMQR)(&side,&trans,&m,&n,&k,a,&lda,tau,c,&ldc,work,&lwork,info);
}

// POTRF
inline void c_xpotrf(const char uplo, const int n, float* a, const int lda, int* info) {
	FORTRAN_NAME(spotrf,SPOTRF)(&uplo,&n,a,&lda,info);
}
inline void c_xpotrf(const char uplo, const int n, scpx* a, const int lda, int* info) {
	FORTRAN_NAME(cpotrf,CPOTRF)(&uplo,&n,a,&lda,info);
}
inline void c_xpotrf(const char uplo, const int n, double* a, const int lda, int* info) {
	FORTRAN_NAME(dpotrf,DPOTRF)(&uplo,&n,a,&lda,info);
}
inline void c_xpotrf(const char uplo, const int n, dcpx* a, const int lda, int* info) {
	FORTRAN_NAME(zpotrf,ZPOTRF)(&uplo,&n,a,&lda,info);
}

// POTRS
inline void c_xpotrs(const char uplo, const int n, const int nrhs, const float* a, const int lda, float* b,
                     const int ldb, int* info) {
	FORTRAN_NAME(spotrs,SPOTRS)(&uplo,&n,&nrhs,a,&lda,b,&ldb,info);
}
inline void c_xpotrs(const char uplo, const int n, const int nrhs, const scpx* a, const int lda, scpx* b,
                     const int ldb, int* info) {
	FORTRAN_NAME(cpotrs,CPOTRS)(&uplo,&n,&nrhs,a,&lda,b,&ldb,info);
}
inline void c_xpotrs(const char uplo, const int n, const int nrhs, const double* a, const int lda, double* b,
                     const int ldb, int* info) {
	FORTRAN_NAME(dpotrs,DPOTRS)(&uplo,&n,&nrhs,a,&lda,b,&ldb,info);
}
inline void c_xpotrs(const char uplo, const int n, const int nrhs, const dcpx* a, const int lda, dcpx* b,
                     const int ldb, int* info) {
	FORTRAN_NAME(zpotrs,ZPOTRS)(&uplo,&n,&nrhs,a,&lda,b,&ldb,info);
}

// GESV
inline void c_xgesv(const int n, const int nrhs, float* a, const int lda, int* ipiv, float* b,
                    const int ldb, int* info) {
	FORTRAN_NAME(sgesv,SGESV)(&n,&nrhs,a,&lda,ipiv,b,&ldb,info);
}
inline void c_xgesv(const int n, const int nrhs, scpx* a, const int lda, int* ipiv, scpx* b,
                    const int ldb, int* info) {
	FORTRAN_NAME(cgesv,CGESV)(&n,&nrhs,a,&lda,ipiv,b,&ldb,info);
}
inline void c_xgesv(const int n, const int nrhs, double* a, const int lda, int* ipiv, double* b,
                    const int ldb, int* info) {
	FORTRAN_NAME(dgesv,DGESV)(&n,&nrhs,a,&lda,ipiv,b,&ldb,info);
}
inline void c_xgesv(const int n, const int nrhs, dcpx* a, const int lda, int* ipiv, dcpx* b,
                    const int ldb, int* info) {
	FORTRAN_NAME(zgesv,ZGESV)(&n,&nrhs,a,&lda,ipiv,b,&ldb,info);
}

// GESVD
inline void c_xgesvd(const char jobu, const char jobvt, const int m, const int n, float* a, const int lda,
                     float* s, float* u, const int ldu, float* vt, const int ldvt, float* work,
                     const int lwork, int* info) {
	FORTRAN_NAME(sgesvd,SGESVD)(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,info);
}
inline void c_xgesvd(const char jobu, const char jobvt, const int m, const int n, scpx* a, const int lda,
                     float* s, scpx* u, const int ldu, scpx* vt, const int ldvt, scpx* work,
                     const int lwork, float* rwork, int* info) {
	FORTRAN_NAME(cgesvd,CGESVD)(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,rwork,info);
}
inline void c_xgesvd(const char jobu, const char jobvt, const int m, const int n, double* a, const int lda,
                     double* s, double* u, const int ldu, double* vt, const int ldvt, double* work,
                     const int lwork, int* info) {
	FORTRAN_NAME(dgesvd,DGESVD)(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,info);
}
inline void c_xgesvd(const char jobu, const char jobvt, const int m, const int n, dcpx* a, const int lda,
                     double* s, dcpx* u, const int ldu, dcpx* vt, const int ldvt, dcpx* work,
                     const int lwork, double* rwork, int* info) {
	FORTRAN_NAME(zgesvd,ZGESVD)(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,rwork,info);
}

// ILAENV
inline int c_ilaenv(const int ispec, const char* name, const char* opts, const int n1, const int n2,
                    const int n3, const int n4) {
	return FORTRAN_NAME(ilaenv,ILAENV)(&ispec,name,opts,&n1,&n2,&n3,&n4);
}

// GEEV
inline void c_xgeev(const char jobvl, const char jobvr, const int n, float* a, const int lda,
					float* wr, float* wi, float* vl, const int ldvl, float* vr,
					const int ldvr, float* work, const int lwork, int* info) {
	FORTRAN_NAME(sgeev,SGEEV)(&jobvl,&jobvr,&n,a,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,info);		
}
inline void c_xgeev(const char jobvl, const char jobvr, const int n, scpx* a, const int lda,
					scpx* w, scpx* vl, const int ldvl, scpx* vr, const int ldvr,
					scpx* work, const int lwork, float* rwork, int* info) {
	FORTRAN_NAME(cgeev,CGEEV)(&jobvl,&jobvr,&n,a,&lda,w,vl,&ldvl,vr,&ldvr,work,&lwork,rwork,info);
}
inline void c_xgeev(const char jobvl, const char jobvr, const int n, double* a, const int lda,
					double* wr, double* wi, double* vl, const int ldvl, double* vr,
					const int ldvr, double* work, const int lwork, int* info) {
	FORTRAN_NAME(dgeev,DGEEV)(&jobvl,&jobvr,&n,a,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,info);
}
inline void c_xgeev(const char jobvl, const char jobvr, const int n, dcpx* a, const int lda,
					dcpx* w, dcpx* vl, const int ldvl, dcpx* vr, const int ldvr,
					dcpx* work, const int lwork, double* rwork, int* info) {
	FORTRAN_NAME(zgeev,ZGEEV)(&jobvl,&jobvr,&n,a,&lda,w,vl,&ldvl,vr,&ldvr,work,&lwork,rwork,info);
}

// SYEV
inline void c_xsyev(const char jobz, const char uplo, const int n, float* a, const int lda,
									float* w, float* work, const int lwork, int* info) {
	FORTRAN_NAME(ssyev,SSYEV)(&jobz,&uplo,&n,a,&lda,w,work,&lwork,info);
}
inline void c_xsyev(const char jobz, const char uplo, const int n, double* a, const int lda,
									double* w, double* work, const int lwork, int* info) {
	FORTRAN_NAME(dsyev,DSYEV)(&jobz,&uplo,&n,a,&lda,w,work,&lwork,info);
}

// HEEV
inline void c_xheev(const char jobz, const char uplo, const int n, scpx* a, const int lda,
									float* w, scpx* work, const int lwork, float* rwork, int* info) {
	FORTRAN_NAME(cheev,CHEEV)(&jobz,&uplo,&n,a,&lda,w,work,&lwork,rwork,info);
}
inline void c_xheev(const char jobz, const char uplo, const int n, dcpx* a, const int lda,
									double* w, dcpx* work, const int lwork, double* rwork, int* info) {
	FORTRAN_NAME(zheev,ZHEEV)(&jobz,&uplo,&n,a,&lda,w,work,&lwork,rwork,info);
}

#endif // _BLAS_
