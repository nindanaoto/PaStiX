/**
 *
 * @file core_zrqrrt.c
 *
 * PaStiX Rank-revealing QR kernel based on randomization technique with rotation.
 *
 * @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Esragul Korkmaz
 * @date 2019-11-12
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <cblas.h>
#include <lapacke.h>
#include "common/frobeniusupdate.h"
#include "blend/solver.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"
#include "z_nan_check.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static pastix_complex64_t zone  = 1.0;
static pastix_complex64_t zzero = 0.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief Compute a randomized QR factorization with rotation technique.
 *
 * This kernel implements the second method (he did not write a pseudocode for
 * the second method) described in :
 *
 * Blocked rank-revealing QR factorizations: How randomized sampling can be used
 * to avoid single-vector pivoting. P. G. Martinsson
 *
 * Note that we only use the trailing matrix for resampling. We don't use power
 * of it for getting better results, since it would be too costly.
 *
 * Difference from randomized QRCP is :
 *   1) resampling is used instead of sample update formula
 *   2) Instead of pivoting A, rotation is applied on it
 *   3) Instead of working with B and omega, directly their transposes are handled
 *    for simplicity
 *
 * The main difference in this algorithm compared to figure 5 in the Martinsson's
 * paper:
 *   1) First, we stop the iterative process based on a tolerance criterion
 *   2) Second, we do not apply SVD for gathering the mass on the diagonal
 *      blocks, so our method is less costly but we expect it to be less
 *      close to SVD result
 *   3) Third, we do not apply the power iteration for resampling for a closer
 *      result to SVD, since it is too costly
 *
 *******************************************************************************
 *
 * @param[in] tol
 *          The relative tolerance criterion. Computations are stopped when the
 *          frobenius norm of the residual matrix is lower than tol.
 *          If tol < 0, then maxrank reflectors are computed.
 *
 * @param[in] maxrank
 *          Maximum number of reflectors computed. Computations are stopped when
 *          the rank exceeds maxrank. If maxrank < 0, all reflectors are computed
 *          or up to the tolerance criterion.
 *
 * @param[in] nb
 *          Tuning parameter for the GEMM blocking size. if nb < 0, nb is set to
 *          32.
 *
 * @param[in] m
 *          Number of rows of the matrix A.
 *
 * @param[in] n
 *          Number of columns of the matrix A.
 *
 * @param[inout] A
 *          The matrix of dimension lda-by-n that needs to be compressed.
 *          On output, the A matrix is computed up to the founded
 *          rank k (k first columns and rows). Everything else, should be dropped.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1, m).
 *
 * @param[out] tau
 *          Contains scalar factors of the elementary reflectors for the matrix
 *          A.
 *
 * @param[out] B
 *          The matrix of dimension ldb-by-maxrank that will holds the partial
 *          projection of the matrix A.
 *          On output, each block of 32 columns can be used to computed the
 *          reverse rotation of the R part of A.
 *
 * @param[in] ldb
 *          The leading dimension of the matrix B. ldb >= max(1, n).
 *
 * @param[out] tau_b
 *          Contains scalar factors of the elementary reflectors for the matrix
 *          B.
 *
 * @param[in] work
 *          Workspace array of size lwork.
 *
 * @param[in] lwork
 *          The dimension of the work area. lwork >= (nb * max(n,n))
 *          If lwork == -1, the function returns immediately and work[0]
 *          contains the optimal size of work.
 *
 * @param[in] normA
 *          The norm of the input matrixA. If negative, the norm is computed by
 *          the kernel.
 *
 *******************************************************************************
 *
 * @return This routine will return the rank of A (>=0) or -1 if it didn't
 *         manage to compress within the margins of tolerance and maximum rank.
 *
 *******************************************************************************/
int
core_zrqrrt( double tol, pastix_int_t maxrank, pastix_int_t nb,
             pastix_int_t m, pastix_int_t n,
             pastix_complex64_t *A, pastix_int_t lda, pastix_complex64_t *tau,
             pastix_complex64_t *B, pastix_int_t ldb, pastix_complex64_t *tau_b,
             pastix_complex64_t *work, pastix_int_t lwork, double normA )
{
    int                 SEED[4] = {26, 67, 52, 197};
    int                 ret, i;
    pastix_int_t        bp = ( nb < 0 ) ? 32 : nb;
    pastix_int_t        d, ib, loop = 1;
    pastix_int_t        ldo = m;
    pastix_int_t        size_O = ldo * bp;
    pastix_int_t        rk, minMN, lwkopt;
    pastix_int_t        sublw = n * bp;
    pastix_complex64_t *omega = work;
    pastix_complex64_t *subw  = work;
    double normR;

    sublw = pastix_imax( sublw, size_O );

    char trans;
#if defined(PRECISION_c) || defined(PRECISION_z)
    trans = 'C';
#else
    trans = 'T';
#endif

    lwkopt = sublw;
    if ( lwork == -1 ) {
        work[0] = (pastix_complex64_t)lwkopt;
        return 0;
    }
#if !defined(NDEBUG)
    if (m < 0) {
        return -1;
    }
    if (n < 0) {
        return -2;
    }
    if (lda < pastix_imax(1, m)) {
        return -4;
    }
    if( lwork < lwkopt ) {
        return -8;
    }
#endif

    minMN = pastix_imin(m, n);
    if ( maxrank < 0 ) {
        maxrank = minMN;
    }
    maxrank = pastix_imin( maxrank, minMN );

    /* Compute the norm of A if not provided by the user */
    if ( normA < 0 ) {
        normA = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                     A, lda, NULL );
    }

    /**
     * If maximum rank is 0, then either the matrix norm is below the tolerance,
     * and we can return a null rank matrix, or it is not and we need to return
     * a full rank matrix.
     */
    if ( maxrank == 0 ) {
        if ( tol < 0. ) {
            return 0;
        }
        if ( normA < tol ) {
            return 0;
        }
        return -1;
    }

    /* Quick exit if A is null rank for the given tolerance */
    if ( normA < tol ) {
        return 0;
    }

#if defined(PASTIX_DEBUG_LR)
    omega = malloc( size_O * sizeof(pastix_complex64_t) );
    subw  = malloc( sublw  * sizeof(pastix_complex64_t) );
#endif

    /* Computation of the Gaussian matrix */
    LAPACKE_zlarnv_work(3, SEED, size_O, omega);

    rk = 0;
    while ( (rk < maxrank) && loop )
    {
        /*
         * Note that we can use maxrank instead of minMN to compute ib, as it is
         * useless to compute extra columns with rotation. The residual will
         * tell us if we managed to compress or not
         */
        ib = pastix_imin( bp, maxrank-rk );
        d = ib;

        /* Computation of the projection matrix B = A_{0:m,k:n}^H * omega */
        cblas_zgemm( CblasColMajor, CblasConjTrans, CblasNoTrans,
                     n-rk, d, m-rk,
                     CBLAS_SADDR(zone),  A + rk*lda + rk, lda,
                                         omega,           ldo,
                     CBLAS_SADDR(zzero), B + rk*ldb + rk, ldb );

        /* Try to do some power iteration to refine the projection */
        if (0)
        {
            for(int l=0; l<2; l++)
            {
                cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                             m-rk, d, n-rk,
                             CBLAS_SADDR(zone),  A + rk*lda + rk, lda,
                                                 B + rk*ldb + rk, ldb,
                             CBLAS_SADDR(zzero), omega,           ldo );

                cblas_zgemm( CblasColMajor, CblasConjTrans, CblasNoTrans,
                             n-rk, d, m-rk,
                             CBLAS_SADDR(zone),  A + rk*lda + rk, lda,
                                                 omega,           ldo,
                             CBLAS_SADDR(zzero), B + rk*ldb + rk, ldb );
            }

            /* Computation of the Gaussian matrix */
            LAPACKE_zlarnv_work(3, SEED, size_O, omega);
        }

        /*
         * QR factorization of the sample matrix B = Q_{B} R_{B}.
         * At the end, the householders will be stored at the lower part of the matrix
         */
        ret = LAPACKE_zgeqrf_work( LAPACK_COL_MAJOR, n-rk, d,
                                   B + rk*ldb + rk, ldb, tau_b+rk,
                                   subw, sublw );
        assert(ret == 0);

        /*
         *  A_{0:m,k:n} = A_{0:m,k:n} Q_{B} for rotational version
         */
        ret = LAPACKE_zunmqr_work( LAPACK_COL_MAJOR, 'R', 'N',
                                   m - rk, n - rk, d,
                                   B + rk*ldb + rk, ldb, tau_b+rk,
                                   A + rk*lda + rk, lda,
                                   subw, sublw );
        assert(ret == 0);

        /*
         * Factorize d columns of A_{k:m,k:k+d} without pivoting
         */
        ret = LAPACKE_zgeqrf_work( LAPACK_COL_MAJOR, m-rk, d,
                                   A + rk*lda + rk, lda, tau + rk,
                                   subw, sublw );
        assert(ret == 0);

        if ( rk+d < n ) {
            /*
             * Update trailing submatrix: A_{k:m,k+d:n} <- Q^h A_{k:m,k+d:n}
             */
            ret = LAPACKE_zunmqr_work( LAPACK_COL_MAJOR, 'L', trans,
                                       m-rk, n-rk-d, d,
                                       A +  rk   *lda + rk, lda, tau + rk,
                                       A + (rk+d)*lda + rk, lda,
                                       subw, sublw );
            assert(ret == 0);
        }

        /* Let's update the residual norm */
        normR = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m-rk-d, n-rk-d, A + (rk+d) * (lda+1), lda, NULL );
        if ( normR < tol ) {
            /* Let's refine the rank: r <= rk+d */
            double ssq = 1.;
            double scl = normR;

            loop = 0;

            for( i=d-1; i>=0; i-- ) {
                double normRk = cblas_dznrm2( n-rk-i, A + (rk+i) * lda + (rk+i), lda );

                frobenius_update( 1., &scl, &ssq, &normRk );
                normRk = scl * sqrt( ssq );

                if ( normRk > tol ) {
                    /*
                     * The actual rank is i (the i^th column has just been
                     * removed from the selection), and we need to be below the
                     * threshold tol, so we need the i from the previous
                     * iteration (+1)
                     */
                    d = i+1;
                    break;
                }
            }
        }
        rk += d;
    }

#if defined(PASTIX_DEBUG_LR)
    free( omega );
    free( subw  );
#endif

    (void)ret;
    if ( (loop && (rk < minMN)) || (rk > maxrank) ) {
        return -1;
    }
    else {
        return rk;
    }
}

/**
 *******************************************************************************
 *
 * @brief Convert a full rank matrix in a low rank matrix, using RQRRT.
 *
 *******************************************************************************
 *
 * @param[in] use_reltol
 *          Defines if the kernel should use relative tolerance (tol *||A||), or
 *          absolute tolerance (tol).
 *
 * @param[in] tol
 *          The tolerance used as a criterion to eliminate information from the
 *          full rank matrix
 *
 * @param[in] rklimit
 *          The maximum rank to store the matrix in low-rank format. If
 *          -1, set to min(m, n) / PASTIX_LR_MINRATIO.
 *
 * @param[in] m
 *          Number of rows of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] n
 *          Number of columns of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] A
 *          The matrix of dimension lda-by-n that needs to be compressed
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1, m)
 *
 * @param[out] Alr
 *          The low rank matrix structure that will store the low rank
 *          representation of A
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zge2lr_rqrrt( int use_reltol, pastix_fixdbl_t tol, pastix_int_t rklimit,
                   pastix_int_t m, pastix_int_t n,
                   const void *A, pastix_int_t lda,
                   pastix_lrblock_t *Alr )
{
    return core_zge2lr_qrrt( core_zrqrrt, use_reltol, tol, rklimit,
                             m, n, A, lda, Alr );
}
