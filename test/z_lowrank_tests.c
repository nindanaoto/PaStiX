/**
 *
 * @file z_lowrank_tests.c
 *
 * Test functions for the low-rank kernels.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2019-12-05
 *
 * @precisions normal z -> z c d s
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdio.h>
#include <cblas.h>
#include <assert.h>
#include <pastix.h>
#include "common/common.h"
#include <lapacke.h>
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"
#include "kernels/pastix_zlrcores.h"
#include "kernels/pastix_lowrank.h"
#include "tests.h"
#include "z_tests.h"

#define PI 3.14159265358979323846

static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t zzero =  0.0;

static int ISEED[4] = { 42, 15, 314, 666 }; /* initial seed for zlarnv() */

/**
 *******************************************************************************
 *
 * @brief Generate random matrices with given singular values in dense format.
 *
 *******************************************************************************
 *
 * @param[in] mode
 *          Defines the distribution of the singular values of the generated
 *          matrix.
 *          There are three experimental cases issued from the paper:
 *          "Martinsson, P. G. (2015). “Blocked rank-revealing QR factorizations
 *          : How randomized sampling can be used to avoid single-vector
 *          pivoting” where we test different singular value distributions based
 *          on mode value:
 *          0) A is an i.i.d normalized Gaussian matrix, so singular values
 *             should decrease slowly and drop at the end
 *          2) A = U D V^H where U,V are random orthonormal matrices. The
 *             singular values are flat up to a point, after decrease rapidly
 *             and then level out again to threshold. (Uses cubic spline
 *             interpolation to describe the curve)
 *          3) A = U D V^H where U,V are random orthonormal matrices. The
 *             singular values are flat up to a point, after decrease rapidly
 *             and then level out again to threshold= tolerance^2. (Uses cosinus
 *             to describe the curve)
 *          -) A = U D V^H where U,V are random orthonormal matrices, and
 *             D(i,i) = max( exp(log(tolerance) / rank) ^ i, threshold )
 *
 * @param[in] tolerance
 *          Define the value of the middle point in the singular values
 *          distribution.
 *
 * @param[in] threshold
 *          Defines the value of the plateau value in the singular values
 *          distribution.
 *
 * @param[in,out] tA
 *          The test matrix to generate.
 *          On entry, m, n, ld, and rk must be defined.
 *          On output, tA->fr is filled with the generated random matrix of size
 *          m-by-n.
 *
 *******************************************************************************/
int
z_lowrank_genmat( int            mode,
                  double         tolerance,
                  double         threshold,
                  test_matrix_t *tA )
{
    pastix_int_t        rank = tA->rk;
    pastix_int_t        m    = tA->m;
    pastix_int_t        n    = tA->n;
    pastix_complex64_t *A    = tA->fr;
    pastix_int_t        lda  = tA->ld;
    pastix_int_t        minMN = pastix_imin( m, n );
    double              rcond = (double)minMN;
    double              dmax  = 1.0;
    pastix_complex64_t *work;
    double *            S;
    int                 i, rc;
    int                 SEED[4] = { 26, 67, 52, 197 };

    if ( m < 0 ) {
        fprintf( stderr, "Invalid m parameter\n" );
        return -4;
    }
    if ( n < 0 ) {
        fprintf( stderr, "Invalid n parameter\n" );
        return -5;
    }
    if ( lda < m ) {
        fprintf( stderr, "Invalid lda parameter\n" );
        return -6;
    }
    if ( rank > pastix_imin( m, n ) ) {
        fprintf( stderr, "Invalid rank parameter\n" );
        return -3;
    }

    if ( rank == 0 ) {
        LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', m, n,
                             0., 0., A, lda );
        tA->norm = 0.;
        return 0;
    }

    S    = malloc( minMN * sizeof( double ) );
    work = malloc( 3 * pastix_imax( m, n ) * sizeof( pastix_complex64_t ) );

    if ( ( !S ) || ( !work ) ) {
        fprintf( stderr, "Out of Memory\n" );
        free( S );
        free( work );
        return -2;
    }

    /*
     * There are three experimental cases issued from the paper: "Martinsson,
     * P. G. (2015). “Blocked rank-revealing QR factorizations : How randomized
     * sampling can be used to avoid single-vector pivoting”
     * where we test different singular value distributions based on mode value:
     *    0) A is an i.i.d normalized Gaussian matrix, so singular values should
     *       decrease slowly and drop at the end
     *    2) A = U D V^H where U,V are random orthonormal matrices. The singular
     *       values are flat up to a point, after decrease rapidly and then
     *       level out again to threshold. (Uses cubic spline interpolation to
     *       describe the curve)
     *    3) A = U D V^H where U,V are random orthonormal matrices. The singular
     *       values are flat up to a point, after decrease rapidly and then
     *       level out again to threshold= tolerance^2. (Uses cosinus to
     *       describe the curve)
     *    -) A = U D V^H where U,V are random orthonormal matrices, and
     *           D(i,i) = max( exp(log(tolerance) / rank) ^ i, threshold )
     */
    switch ( mode ) {
        case 0: {
            /* Generate  a random matrix of rank rank */
            pastix_complex64_t *O1     = malloc( m * rank * sizeof( pastix_complex64_t ) );
            pastix_complex64_t *O2     = malloc( n * rank * sizeof( pastix_complex64_t ) );
            pastix_complex64_t *A2     = malloc( n * lda * sizeof( pastix_complex64_t ) );
            double *            superb = malloc( minMN * sizeof( double ) );
            double              alpha;

            LAPACKE_zlarnv_work( 3, SEED, m * rank, O1 );
            LAPACKE_zlarnv_work( 3, SEED, n * rank, O2 );
            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         m, n, rank,
                         CBLAS_SADDR( zone ),  O1, m,
                                               O2, pastix_imax(rank, 1),
                         CBLAS_SADDR( zzero ), A,  lda );

            LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, A2, lda );
            LAPACKE_zgesvd(
                LAPACK_COL_MAJOR, 'n', 'n', m, n, A2, lda, S, NULL, 1, NULL, 1, superb );

            LAPACKE_zlascl_work( LAPACK_COL_MAJOR, 'g', -1, -1, 1., S[0], m, n, A, lda );
            alpha = 1. / S[0];
            cblas_dscal( minMN, alpha, S, 1 );

            free( superb );
            free( A2 );
            free( O1 );
            free( O2 );
        }
        break;

        case 2: {
            int    s1 = pastix_imax( 0, (rank / 2) - 1 );
            int    s2 = pastix_imin( pastix_imax( rank + ( rank / 2 ) - 1, rank ),
                                     minMN );
            double p1 = log( 1. / tolerance );
            double p2 = log( threshold / tolerance );

            /**
             * a3( p1, p2 ) = abs(p1) <= abs(p2) ? -     p1           :  3  * p2 + 2. * p1
             * a2( p1, p2 ) = abs(p1) <= abs(p2) ? -3. * p1           :  6. * p2 + 3. * p1
             * a1( p1, p2 ) = abs(p1) <= abs(p2) ? -3. * p1           :  3. * p2
             * b3( p1, p2 ) = abs(p1) <= abs(p2) ? -3  * p1 - 2. * p2 :       p2
             * b2( p1, p2 ) = abs(p1) <= abs(p2) ?  6. * p1 + 3. * p2 : -3. * p2
             * b1( p1, p2 ) = abs(p1) <= abs(p2) ? -3. * p1           :  3. * p2
             *
             * A(p1, p2, t) = a3(p1,p2) * t*t*t + a2(p1,p2) *t*t + a1(p1, p2) * t
             * B(p1, p2, t) = b3(p1,p2) * t*t*t + b2(p1,p2) *t*t + b1(p1, p2) * t
             *
             * t(x) = (x-100) / 50
             * f(p1, p2, x) = x < 100 ? A(p1, p2, t(x)) : B(p1, p2, t(x))
             */
            double a1, a2, a3, b1, b2, b3;
            if ( fabs( p1 ) < fabs( p2 ) ) {
                a3 = -p1;
                a2 = -3. * p1;
                a1 = -3. * p1;
                b3 = -3 * p1 - 2. * p2;
                b2 = 6. * p1 + 3. * p2;
                b1 = -3. * p1;
            }
            else {
                a3 = 3 * p2 + 2. * p1;
                a2 = 6. * p2 + 3. * p1;
                a1 = 3. * p2;
                b3 = p2;
                b2 = -3. * p2;
                b1 = 3. * p2;
            }

            for ( i = 0; i <= s1; i++ ) {
                S[i] = 1.;
            }

            for ( i = s1; i <= rank-1; i++ ) {
                double x  = ( (double)( 2 * i - s1 - s2 ) / (double)( s2 - s1 ) );
                double Ax = x * ( a1 + x * ( a2 + a3 * x ) );
                S[i]      = tolerance * exp( Ax );
            }

            for ( ; i < s2; i++ ) {
                double x  = ( (double)( 2 * i - s1 - s2 ) / (double)( s2 - s1 ) );
                double Bx = x * ( b1 + x * ( b2 + b3 * x ) );
                S[i]      = tolerance * exp( Bx );
            }

            for ( ; i < minMN; i++ ) {
                S[i] = threshold;
            }
        }
        break;

        case 3: {
            int    s1    = rank / 2;
            int    s2    = pastix_imin( rank + ( rank / 2 ), minMN );
            double tol2  = tolerance * tolerance;
            double alpha = exp( 2. * log( tolerance ) / rank );

            for ( i = 0; i <= s1; i++ ) {
                S[i] = 1.;
            }

            for ( ; i <= pastix_imin(s2,minMN-1); i++ ) {
                S[i] = S[i - 1] * alpha;
            }

            for ( i = s1; i < s2; i++ ) {
                double x    = ( (double)( i - s1 ) / (double)rank );
                double cosx = -cos( x * PI );
                S[i]        = tolerance * exp( log( tolerance ) * cosx );
            }

            for ( ; i < minMN; i++ ) {
                S[i] = tol2;
            }
        }
        break;

        default: {
            double alpha;
            if ( rank == 0 ) {
                S[0]  = 0.;
                alpha = 1;
            }
            else {
                alpha = exp( log( tolerance ) / rank );
                S[0] = alpha;
            }
            for ( i = 1; i < minMN; i++ ) {
                if ( rank == i ) {
                    alpha = exp( 2 * log( threshold/tolerance ) / rank );
                }
                S[i] = S[i - 1] * alpha;
                if ( S[i] <= threshold ) {
                    alpha = 1.;
                }
            }
        }
    }

    /* Initialize A */
    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, m, n, 'U', ISEED, 'N', S, 0,
                         rcond, dmax, m, n, 'N', A, lda, work );

    tA->norm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n, A, lda, NULL );

    /* Let's store the singular values in a file */
    {
        char *filename;
        FILE *f;

        rc = asprintf( &filename,
                       "singular_values_N%d_mode%d_rank%d_tol%e.txt",
                       (int)minMN, (int)mode, (int)rank, tolerance );
        f = fopen( filename, "w" );
        free( filename );

        for ( i = 0; i < minMN; i++ ) {
            fprintf( f, "%4d %e %e\n", (int)(i+1), S[i], cblas_dnrm2( minMN - i, S + i, 1 ) / (tA->norm) );
        }
        fclose( f );
    }

    free( S );
    free( work );
    (void)rc;
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Generate a random matrix with given singular values in both dense and
 * low-rank format.
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The data structure that defines the kernels used for the
 *          compression. It also defines, the tolerance of the low-rank
 *          representation and if absolute or relative tolerance is applied.
 *
 * @param[in] mode
 *          Defines the distribution of the singular values of the generated
 *          matrix.
 *          There are three experimental cases issued from the paper:
 *          "Martinsson, P. G. (2015). “Blocked rank-revealing QR factorizations
 *          : How randomized sampling can be used to avoid single-vector
 *          pivoting” where we test different singular value distributions based
 *          on mode value:
 *          0) A is an i.i.d normalized Gaussian matrix, so singular values
 *             should decrease slowly and drop at the end
 *          2) A = U D V^H where U,V are random orthonormal matrices. The
 *             singular values are flat up to a point, after decrease rapidly
 *             and then level out again to threshold. (Uses cubic spline
 *             interpolation to describe the curve)
 *          3) A = U D V^H where U,V are random orthonormal matrices. The
 *             singular values are flat up to a point, after decrease rapidly
 *             and then level out again to threshold= tolerance^2. (Uses cosinus
 *             to describe the curve)
 *          -) A = U D V^H where U,V are random orthonormal matrices, and
 *             D(i,i) = max( exp(log(tolerance) / rank) ^ i, threshold )
 *
 * @param[in] threshold
 *          Defines the value of the plateau value in the singular values
 *          distribution.
 *
 * @param[in,out] A
 *          The test matrix to generate.
 *          On entry, m, n, ld, and rk must be defined.
 *          On output, A->fr is filled with the generated random matrix of size
 *          m-by-n, A->lr is the low-rank representation of A->fr with the
 *          compression parameters from lowrank.
 *
 *******************************************************************************/
void
z_lowrank_genmat_comp( const pastix_lr_t *lowrank,
                       int                mode,
                       double             threshold,
                       test_matrix_t     *A )
{
    A->fr = malloc( A->ld * A->n * sizeof(pastix_complex64_t) );
    z_lowrank_genmat( mode, lowrank->tolerance, threshold, A );
    lowrank->core_ge2lr( lowrank->use_reltol, lowrank->tolerance, pastix_imin( A->m, A->n ),
                         A->m, A->n, A->fr, A->ld, &(A->lr) );
}

/**
 *******************************************************************************
 *
 * @brief Check that a compression kernel is correct.
 *
 * This routine compress a dense matrix into a low-rank form according to the
 * lowrank parameters, and check that the compressed form of the matrix respects
 * the relation:
 * \[ || A - c(A) || < \tau \], with
 *     \[ \tau = \eps \] in absolute compression, and
 *     \[ \tau = \eps * ||A|| \] in relative compression.
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The data structure that defines the kernels used for the
 *          compression. It also defines, the tolerance of the low-rank
 *          representation and if absolute or relative tolerance is applied.
 *
 * @param[in,out] A
 *          The matrix to test.
 *          On entry, m, n, ld, rk, and fr must be defined.
 *          On output, A is unchanged, but A->lr is used during the test to
 *          temporarily store the compressed form of the matrix.
 *
 *******************************************************************************
 *
 * @retval 0 on success
 * @retval <0, if one of the parameter is incorrect
 * @retval >0, if one or more of the tests failed.
 *
 *******************************************************************************/
int
z_lowrank_check_ge2lr( pastix_lr_t   *lowrank,
                       test_matrix_t *A )
{
    pastix_int_t        m  = A->m;
    pastix_int_t        n  = A->n;
    pastix_int_t        ld = A->ld;
    pastix_complex64_t *A2;
    pastix_int_t        minMN = pastix_imin( m, n );
    Clock               timer, total_timer = 0.;
    double              resid, errbound, norm_residual;
    pastix_int_t        i, rankA, nb_tests = 1;
    int                 rc = 0;

    if ( m < 0 ) {
        fprintf( stderr, "Invalid m parameter\n" );
        return -4;
    }
    if ( n < 0 ) {
        fprintf( stderr, "Invalid n parameter\n" );
        return -5;
    }
    if ( ld < m ) {
        fprintf( stderr, "Invalid ld parameter\n" );
        return -6;
    }
    A2 = malloc( n * ld * sizeof( pastix_complex64_t ) );

    for ( i = 0; i < nb_tests; i++ ) {
        /*
         * Compress and then uncompress
         */
        timer = clockGetLocal();
        lowrank->core_ge2lr( lowrank->use_reltol, lowrank->tolerance, minMN,
                             m, n, A->fr, ld, &(A->lr) );
        timer = clockGetLocal() - timer;
        assert( timer >= 0. );
        total_timer += timer;

        /*
         * Let's check the result
         */
        core_zlr2ge( PastixNoTrans, m, n, &(A->lr), A2, ld );

        core_zgeadd( PastixNoTrans, m, n,
                     -1., A->fr, ld,
                      1., A2,    ld );

        norm_residual = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                             A2, ld, NULL );

        errbound = lowrank->tolerance;
        if ( (lowrank->use_reltol) && (A->norm > 0.) ) {
            errbound *= A->norm;
        }
        if ( A->lr.rk > 0 ) {
            errbound *= (double)(A->lr.rk);
        }
        resid = norm_residual / errbound;
        if ( resid > 10. ) {
            rc += 1;
        }

        rankA = A->lr.rk;
        core_zlrfree( &(A->lr) );
    }

    total_timer /= nb_tests;

    fprintf( stdout, "%7s %4d %e %e %e %e %7s\n",
             compmeth_shnames[lowrank->compress_method], (int)rankA, total_timer, A->norm, norm_residual, resid,
             (rc != 0) ? "FAILED" : "SUCCESS" );

    free( A2 );

    return rc;
}

/**
 *******************************************************************************
 *
 * @brief Check the operation \[ B = \alpha A + B \]
 *
 * This routine checks the addition of two low-rank matrices together.
 *
 *  Let's have \[ A = c(A) + r(A) \], such that \[ || A - c(A) || < \tau_A \],
 *             \[ B = c(B) + r(B) \], such that \[ || B - c(B) || < \tau_B \]
 *
 * We want to test the correctness of \[ c(C) = \alpha c(A) + c(B) \] with respect to
 * \[ C = \alpha A + B \]
 *
 * For that let's check:
 *    \[ || C - c(C) || = || \alpha A + B - \alpha c(A) - c(B) || = || \alpha r(A) + r(B) || \]
 * So, we consider that the test is correct if:
 *   \[ || C - c(C) || < |\alpha| * ||r(A)|| + ||r(B)|| = \alpha * \tau_A + \tau_B \]
 *
 * Note that if we consider absolute compression, \[ \tau_A = \eps \],
 * otherwise \[ \tau_A = \eps * ||A|| \], and respectively for B.
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The data structure that defines the kernels used for the
 *          compression. It also defines, the tolerance of the low-rank
 *          representation and if absolute or relative tolerance is applied.
 *
 * @param[in] offx
 *          The row offset of the matrix A into the matrix B.
 *          offx >= 0 and A->m + offx <= B->m.
 *
 * @param[in] offy
 *          The column offset of the matrix A into the matrix B.
 *          offy >= 0 and A->n + offy <= B->n.
 *
 * @param[in] alpha
 *          The scalar alpha that multiplies A.
 *
 * @param[in] A
 *          The input matrix A to test.
 *          On entry, m, n, ld, rk, fr, and lr must be defined.
 *          A-lr must be the lowrank representation of A->lr.
 *
 * @param[in] B
 *          The input matrix B to test.
 *          On entry, m, n, ld, rk, fr, and lr must be defined.
 *          B-lr must be the lowrank representation of B->lr.
 *
 * @param[in,out] C
 *          The matrix to test.
 *          On entry, m, n, ld, rk, and fr must be defined. C->fr must be set
 *          as the result of \alpha A + B in the full-rank form.
 *          During the test, C->lr is modified to store the low-rank form of
 *          alpha c(A) + c(B), and is freed on exit.
 *
 *******************************************************************************
 *
 * @retval 0 on success
 * @retval <0, if one of the parameter is incorrect
 * @retval >0, if one or more of the tests failed.
 *
 *******************************************************************************/
int
z_lowrank_check_rradd( pastix_lr_t         *lowrank,
                       pastix_int_t         offx,
                       pastix_int_t         offy,
                       pastix_complex64_t   alpha,
                       const test_matrix_t *A,
                       const test_matrix_t *B,
                       test_matrix_t       *C )
{
    pastix_complex64_t *Clr;
    double              norm_diff, res;
    Clock               timer;
    pastix_int_t        rkCmax;
    pastix_int_t        rankmax = core_get_rklimit( B->m, B->n );
    int                 rc = 0;

    assert( C->m == B->m );
    assert( C->n == B->n );
    assert( (A->m + offx) <= B->m );
    assert( (A->n + offy) <= B->n );

    rkCmax = A->lr.rk + B->lr.rk;

    if ( ( A->lr.rk == -1 ) || ( B->lr.rk == -1 ) ) {
        printf( "Operation not supported\n" );
        return 0;
    }

    /* Init lrC */
    memset( &(C->lr), 0, sizeof( pastix_lrblock_t ) );

    /* Backup B into C */
    core_zlrcpy( NULL, PastixNoTrans, 1.,
                 B->m, B->n, &(B->lr),
                 B->m, B->n, &(C->lr), 0, 0 );

    /* Perform C = B - A in LR format */
    timer = clockGetLocal();
    lowrank->core_rradd( lowrank, PastixNoTrans, &alpha,
                         A->m, A->n, &(A->lr),
                         C->m, C->n, &(C->lr),
                         offx, offy );
    timer = clockGetLocal() - timer;

    /*
     * Check || (\alpha * A+B) - c( \alpha * c(A) + c(B) ) || < tol * (|\alpha| + 1)         (Absolute)
     *    or || (\alpha * A+B) - c( \alpha * c(A) + c(B) ) || < tol * (|\alpha| ||A||+||B||) (Relative)
     */
    Clr = malloc( C->m * C->n * sizeof( pastix_complex64_t ) );

    /* Uncompress the low-rank sum */
    core_zlr2ge( PastixNoTrans, C->m, C->n,
                 &(C->lr), Clr, C->m );

    /* Compute the diff */
    core_zgeadd( PastixNoTrans, C->m, C->n,
                 -1., C->fr, C->ld,
                  1., Clr,   C->m );

    norm_diff = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', C->m, C->n,
                                     Clr, C->m, NULL );

    if ( ( A->lr.rk != 0 ) || ( B->lr.rk != 0 ) ) {
        double errbound = lowrank->use_reltol ? ( cabs(alpha) * A->norm + B->norm ) : cabs(alpha) + 1.;
        res = norm_diff / ( lowrank->tolerance * errbound );
    }
    else {
        res = norm_diff;
    }

    free( Clr );

    /* Check the correctness of the result */
    if ( res > 10.0 ) {
        rc += 1;
    }

    /* Check that final matrix is not full rank, we should have exited before */
    if ( ( C->lr.rk == -1 ) && ( rkCmax <= rankmax ) ) {
        rc += 2;
    }

    /* Check that final rank does not exceed the sum of the ranks */
    if ( C->lr.rk > rkCmax ) {
        rc += 4;
    }

    /* Check that final rank is larger than the minimal rank (rA, rB, abs(rA-rb)) */
    if ( ( C->lr.rk != -1 ) &&
         ( C->lr.rk < pastix_imin( pastix_imin( A->lr.rk, B->lr.rk ), abs( A->lr.rk - B->lr.rk ) ) ) )
    {
        rc += 8;
    }

    fprintf( stdout, "%7s %4d %e %e %e %e ",
             compmeth_shnames[lowrank->compress_method],
             (int)C->lr.rk, clockVal(timer), C->norm, norm_diff, res );

    if ( rc == 0 ) {
        fprintf( stdout, "SUCCESS\n" );
    }
    else {
        fprintf( stdout, "FAILED(%d)\n", rc );
    }

    core_zlrfree( &(C->lr) );
    return ( rc > 0 ) ? 1 : 0;
}

/**
 *******************************************************************************
 *
 * @brief Check the operation \[ C = \alpha A B + \beta C \]
 *
 * This routine checks the addition of two low-rank matrices together.
 *
 *  Let's have \[ A = c(A) + r(A) \], such that \[ || A - c(A) || < \tau_A \]
 *
 * We want to test the correctness of \[ \alpha c(A) c(B) + \beta c(C) \] with respect to
 * \[ \alpha A B + \beta C \]
 *
 * For that let's check:
 *    \[ || C - c(C) || = \alpha * [ c(A)r(B) + r(A)c(B) + rA)r(B) ] + \beta r(C) \]
 * So, we consider that the test is correct if:
 *   \[ || C - c(C) || < |\alpha| * ... + \beta \tau_C \]
 *
 * Note that if we consider absolute compression, \[ \tau_A = \eps \],
 * otherwise \[ \tau_A = \eps * ||A|| \].
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The data structure that defines the kernels used for the
 *          compression. It also defines, the tolerance of the low-rank
 *          representation and if absolute or relative tolerance is applied.
 *
 * @param[in] offx
 *          The row offset of the matrix A into the matrix B.
 *          offx >= 0 and A->m + offx <= B->m.
 *
 * @param[in] offy
 *          The column offset of the matrix A into the matrix B.
 *          offy >= 0 and A->n + offy <= B->n.
 *
 * @param[in] alpha
 *          The scalar alpha that multiplies A B.
 *
 * @param[in] A
 *          The input matrix A to test.
 *          On entry, m, n, ld, rk, fr, and lr must be defined.
 *          A-lr must be the lowrank representation of A->lr.
 *
 * @param[in] B
 *          The input matrix B to test.
 *          On entry, m, n, ld, rk, fr, and lr must be defined.
 *          B-lr must be the lowrank representation of B->lr.
 *
 * @param[in] beta
 *          The scalar beta that multiplies C.
 *
 * @param[in] C
 *          The matrix to test.
 *          On entry, m, n, ld, rk, and fr must be defined. C->fr must be set
 *          as the result of \alpha A B + \beta C in the full-rank form.
 *          C->lr must store the low-rank form of the matrix C before any
 *          computation.
 *
 *******************************************************************************
 *
 * @retval 0 on success
 * @retval <0, if one of the parameter is incorrect
 * @retval >0, if one or more of the tests failed.
 *
 *******************************************************************************/
int
z_lowrank_check_lrmm( pastix_lr_t         *lowrank,
                      pastix_int_t         offx,
                      pastix_int_t         offy,
                      pastix_complex64_t   alpha,
                      const test_matrix_t *A,
                      const test_matrix_t *B,
                      pastix_complex64_t   beta,
                      const test_matrix_t *C )
{
    pastix_lrblock_t    lrC2;
    pastix_complex64_t *Clr;
    double              norm_diff, res;
    Clock               timer;
    int                 rc = 0;

    /* Init lrC */
    memset( &lrC2, 0, sizeof( pastix_lrblock_t ) );

    /* Backup C into C2 */
    core_zlrcpy( NULL, PastixNoTrans, 1., C->m, C->n, &(C->lr), C->m, C->n, &lrC2, 0, 0 );

    /* Compute the low-rank matrix-matrix */
    {
        pastix_atomic_lock_t lock = PASTIX_ATOMIC_UNLOCKED;
        core_zlrmm_t         zlrmm_params;
        zlrmm_params.lowrank = lowrank;
        zlrmm_params.transA  = PastixNoTrans;
        zlrmm_params.transB  = PastixConjTrans;
        zlrmm_params.M       = A->m;
        zlrmm_params.N       = B->m;
        zlrmm_params.K       = A->n;
        zlrmm_params.Cm      = C->m;
        zlrmm_params.Cn      = C->n;
        zlrmm_params.offx    = offx;
        zlrmm_params.offy    = offy;
        zlrmm_params.alpha   = alpha;
        zlrmm_params.A       = &(A->lr);
        zlrmm_params.B       = &(B->lr);
        zlrmm_params.beta    = beta;
        zlrmm_params.C       = &lrC2;
        zlrmm_params.work    = NULL;
        zlrmm_params.lwork   = -1;
        zlrmm_params.lwused  = -1;
        zlrmm_params.lock    = &lock;

        timer = clockGetLocal();
        core_zlrmm( &zlrmm_params );
        timer = clockGetLocal() - timer;
    }

    /*
     * Check || (A+B) - (c(A)+c(B)) || < tol * || A+B ||
     * What we really want is || (A+B) - c(A+B) || < tol * || A+B ||
     */
    Clr = malloc( C->m * C->n * sizeof( pastix_complex64_t ) );

    /* Uncompress the low-rank sum */
    core_zlr2ge( PastixNoTrans, C->m, C->n,
                 &lrC2, Clr, C->ld );

    /* Compute the diff */
    core_zgeadd( PastixNoTrans, C->m, C->n,
                 -1., C->fr, C->ld,
                  1., Clr,   C->ld );

    norm_diff = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', C->m, C->n,
                                     Clr, C->ld, NULL );

    if ( ( C->lr.rk != 0.0 ) && ( ( A->lr.rk + B->lr.rk ) != 0 ) ) {
        double errbound = lowrank->use_reltol ? C->norm : 1.;
        res = norm_diff / ( lowrank->tolerance * errbound );
    }
    else {
        res = norm_diff;
    }

    fprintf( stdout, "%7s %4d %e %e %e %e ",
             compmeth_shnames[lowrank->compress_method],
             (int)lrC2.rk, clockVal(timer), C->norm, norm_diff, res );

    free( Clr );
    core_zlrfree( &lrC2 );

    /* Check the correctness of the result */
    if ( res > 10.0 ) {
        rc += 1;
    }

    if ( rc == 0 ) {
        fprintf( stdout, "SUCCESS\n" );
    }
    else {
        fprintf( stdout, "FAILED(%d)\n", rc );
    }

    return rc;
}
