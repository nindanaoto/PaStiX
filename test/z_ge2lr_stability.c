/**
 *
 * @file z_ge2lr_stability.c
 *
 * Tests and validate the Xge2lr routine.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Esragul Korkmaz
 * @date 2019-11-12
 *
 * @precisions normal z -> z c d s
 *
 **/

#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <pastix.h>
#include "common/common.h"
#include <lapacke.h>
#include <cblas.h>
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"
#include "kernels/pastix_zlrcores.h"
#include "kernels/pastix_lowrank.h"
#include "flops.h"
#include "z_tests.h"
#include "tests.h"

static pastix_complex64_t mzone = -1.0;

/**
 *******************************************************************************
 *
 * @brief Compress a dense matrix with all the vectors to print the decrease of
 * the residual norm.
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The data structure that defines the kernels used for the
 *          compression. It also defines, the tolerance of the low-rank
 *          representation and if absolute or relative tolerance is applied.
 *
 * @param[in] A
 *          The test matrix to study.
 *          On entry, m, n, ld, rk, and fr must be defined.
 *
 * @param[in,out] output
 *          The output file to which the residuals norm are dumped into.
 *
 *******************************************************************************
 *
 * @retval 0 on success
 * @retval <0, if one of the parameter is incorrect
 * @retval >0, if one or more of the tests failed.
 *
 *******************************************************************************/
int
z_lowrank_stability_ge2lr( const pastix_lr_t   *lowrank,
                           const test_matrix_t *A,
                           FILE                *output )
{
    pastix_lrblock_t    lrA;
    pastix_complex64_t *A2;
    pastix_complex64_t *u, *v;
    pastix_int_t m     = A->m;
    pastix_int_t n     = A->n;
    pastix_int_t lda   = A->ld;
    pastix_int_t minMN = pastix_imin(m, n);
    pastix_int_t i, ldu, ldv;
    double norm_residual;

    if (m < 0) {
        fprintf(stderr, "Invalid m parameter\n");
        return -4;
    }
    if (n < 0) {
        fprintf(stderr, "Invalid n parameter\n");
        return -5;
    }
    if (lda < m) {
        fprintf(stderr, "Invalid lda parameter\n");
        return -6;
    }

    /* Backup A into A2 */
    A2 = malloc( m * n * sizeof(pastix_complex64_t));
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                         A->fr, lda, A2, m );

    /*
     * Fully compress the matrix A
     */
    lowrank->core_ge2lr( 1, -1., minMN, m, n, A->fr, lda, &lrA );

    /* Let's check we have the maximal rank */
    assert( lrA.rk == minMN );

    /*
     * Let's compute the frobenius norm of A - U[:,1:i] * V[:,1:i]^T for i in [1:minMN]
     */
    u   = lrA.u;
    v   = lrA.v;
    ldu = m;
    ldv = lrA.rkmax;
    for(i=0; i<minMN; i++) {
        norm_residual = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                             A2, m, NULL );

        fprintf( output, "%4d %e\n", (int)(i+1), norm_residual/A->norm );

	/* Subtract the i^th outer product */
        cblas_zgerc( CblasColMajor, m, n,
                     CBLAS_SADDR(mzone),
                     u + ldu * i, 1,
                     v + i,       ldv,
                     A2,          m );
    }

    core_zlrfree(&lrA);
    free(A2);

    return 0;
}

int main( int argc, char **argv )
{
    test_matrix_t A;
    pastix_int_t n;
    int mode, p, i, rc;
    double tolerance, threshold;
    test_param_t params;
    double eps = LAPACKE_dlamch_work('e');
    pastix_lr_t lowrank;

    testGetOptions( argc, argv, &params, eps );

    lowrank.compress_when       = PastixCompressWhenEnd;
    lowrank.compress_method     = PastixCompressMethodPQRCP;
    lowrank.compress_min_width  = 0;
    lowrank.compress_min_height = 0;
    lowrank.use_reltol          = 1;
    lowrank.tolerance           = params.tol_gen;
    lowrank.core_ge2lr          = core_zge2lr_svd;
    lowrank.core_rradd          = core_zrradd_svd;

    tolerance = params.tol_gen;
    threshold = params.threshold;
    for (n=params.n[0]; n<=params.n[1]; n+=params.n[2]) {
        A.m  = n;
        A.n  = n;
        A.ld = n;
        A.fr = malloc( A.ld * A.n * sizeof(pastix_complex64_t) );

        for (p=params.prank[0]; p<=params.prank[1]; p+=params.prank[2]) {
            A.rk = (p * n) / 100;

            for (mode=params.mode[0]; mode<=params.mode[1]; mode+=params.mode[2]) {
                printf( "   -- Test GE2LR Tol=%e M=N=LDA=%ld R=%ld MODE=%d\n",
                        tolerance, (long)A.n, (long)A.rk, mode );

                /*
                 * Generate a matrix of a given rank for the prescribed tolerance
                 */
                z_lowrank_genmat( mode, tolerance, threshold, &A );

                /* Let's test all methods we have */
                for(i=params.method[0]; i<=params.method[1]; i+=params.method[2])
                {
                    FILE *f;
                    char *filename;

                    lowrank.compress_method = i;
                    lowrank.core_ge2lr = ge2lrMethods[i][PastixComplex64-2];
                    lowrank.core_rradd = rraddMethods[i][PastixComplex64-2];

                    rc = asprintf( &filename,
                                   "stability_m%d_n%d_mode%d_tol%e_prank%d_%s.log",
                                   (int)A.m, (int)A.n, mode, tolerance, (int)A.rk,
                                   compmeth_shnames[i]);
                    f = fopen( filename, "w" );
                    free(filename);

                    fprintf( f,
                             "# method = %s\n"
                             "# mode = %d\n"
                             "# M = %d\n"
                             "# N = %d\n"
                             "# tol = %e\n"
                             "# ||A|| = %e\n",
                             compmeth_lgnames[i], mode,
                             (int)A.m, (int)A.n, tolerance, A.norm );

                    z_lowrank_stability_ge2lr( &lowrank, &A, f );
                    fclose(f);
                }
            }
        }

        free(A.fr);
    }
    (void)rc;
}
