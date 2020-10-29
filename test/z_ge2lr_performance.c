/**
 *
 * @file z_ge2lr_performance.c
 *
 * Tests and validate the Xge2lr routine.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Mathieu Faverge
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

/**
 *******************************************************************************
 *
 * @brief Compress a dense matrix with all the vectors to print the decrease of
 * the residual norm.
 *
 *******************************************************************************
 *
 * @param[in,out] f
 *          The output file to which the results are printed.
 *
 * @param[in] nbruns
 *          The number of times the compression is performed ot get an average
 *          result. nbruns > 1
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
 *******************************************************************************
 *
 * @retval 0 on success
 * @retval <0, if one of the parameter is incorrect
 * @retval >0, if one or more of the tests failed.
 *
 *******************************************************************************/
int
z_lowrank_ge2lr_performance( FILE *f, int nbruns,
                             const pastix_lr_t   *lowrank,
                             const test_matrix_t *A )
{
    pastix_complex64_t *A2;
    pastix_lrblock_t    lrA;
    pastix_int_t m     = A->m;
    pastix_int_t n     = A->n;
    pastix_int_t lda   = A->ld;
    pastix_int_t minMN = pastix_imin(m, n);
    pastix_fixdbl_t flops, gflops;
    double resid, normR;
    Clock timer, total_timer = 0.;
    int i;

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

    /* Backup A in A2 */
    A2 = malloc( m * n * sizeof(pastix_complex64_t));
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                         A->fr, lda, A2, m );

    memset( &lrA, 0, sizeof(pastix_lrblock_t) );

    /* Compress A */
    nbruns = pastix_imax( nbruns, 1 );
    for (i=0; i<nbruns; i++) {
        core_zlrfree(&lrA);

        clockStart(timer);
        lowrank->core_ge2lr( lowrank->use_reltol, lowrank->tolerance,
                             minMN, m, n, A->fr, lda, &lrA );
        clockStop(timer);
        assert( timer >= 0. );
        total_timer += timer;
    }

    flops = FLOPS_ZGEQRF( m, A->rk ) +
        FLOPS_ZUNMQR( m, n - A->rk, A->rk, PastixLeft ) +
        FLOPS_ZUNGQR( m, A->rk, A->rk);
    gflops = (nbruns * flops * 1.e-9) / total_timer;

    /*
     * Let's check the last result
     */
    core_zlr2ge( PastixNoTrans, m, n,
                 &lrA, A2, lda );

    core_zgeadd( PastixNoTrans, m, n,
                 -1., A->fr, lda,
                  1., A2,    lda );

    /* Frobenius norm of ||A - (U_i *V_i)|| */
    normR = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n, A2, lda, NULL );
    resid = lowrank->tolerance;
    if ( (lowrank->use_reltol) && (A->norm > 0.) ) {
        resid *= A->norm;
    }
    if ( lrA.rk > 0 ) {
        resid *= (double)(lrA.rk);
    }
    resid = normR / resid;

    if ( f == stdout ) {
        fprintf( f, " %5d %e %e %e %e %s\n",
                 lrA.rk, total_timer/nbruns, gflops, normR, resid,
                 (resid > 10.) ? "FAILED" : "SUCCESS" );
    }
    else {
        fprintf( f, "%d;%e;%e;%e;%e;%s\n",
                 lrA.rk, total_timer/nbruns, gflops, normR, resid,
                 (resid > 10.) ? "FAILED" : "SUCCESS" );
    }

    free(A2);
    core_zlrfree(&lrA);
    return (resid > 10.);
}

int main( int argc, char **argv )
{
    test_matrix_t A;
    pastix_int_t n;
    int mode, p, i, ret, rc = 0;
    test_param_t params;
    double eps = LAPACKE_dlamch_work('e');
    pastix_lr_t lowrank;

    testGetOptions( argc, argv, &params, eps );

    if ( params.output == stdout ) {
        fprintf( params.output,
                 "%7s %5s %4s %12s %6s %5s %12s %12s %5s %12s %12s %12s %12s\n",
                 "Method", "PRank", "Mode", "TolGen",
                 "N", "Rank", "NormA",
                 "TolCmp", "CRank", "Time", "GFlops", "||A-UVt||_f", "||A-UVt||_f/(||A||_f * eps)" );
    }
    else {
        fprintf( params.output,
                 "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;Check\n",
                 "Method", "PRank", "Mode", "TolGen",
                 "N", "Rank", "NormA",
                 "TolCmp", "CRank", "Time", "GFlops", "||A-UVt||_f", "||A-UVt||_f/(||A||_f * eps)" );
    }

    lowrank.compress_when       = PastixCompressWhenEnd;
    lowrank.compress_method     = PastixCompressMethodPQRCP;
    lowrank.compress_min_width  = 0;
    lowrank.compress_min_height = 0;
    lowrank.use_reltol          = params.use_reltol;
    lowrank.tolerance           = params.tol_cmp;
    lowrank.core_ge2lr          = core_zge2lr_svd;
    lowrank.core_rradd          = core_zrradd_svd;

    for (n=params.n[0]; n<=params.n[1]; n+=params.n[2]) {
        A.m  = n;
        A.n  = n;
        A.ld = n;
        A.fr = malloc( A.ld * A.n * sizeof(pastix_complex64_t) );

        for (p=params.prank[0]; p<=params.prank[1]; p+=params.prank[2]) {
            A.rk = (p * n) / 100;

            for (mode=params.mode[0]; mode<=params.mode[1]; mode+=params.mode[2])
            {
                /*
                 * Generate a matrix of a given rank for the prescribed tolerance
                 */
                z_lowrank_genmat( mode, params.tol_gen,
                                  params.threshold, &A );

                /* Let's test all methods we have */
                for(i=params.method[0]; i<=params.method[1]; i+=params.method[2])
                {
                    lowrank.compress_method = i;
                    lowrank.core_ge2lr = ge2lrMethods[i][PastixComplex64-2];
                    lowrank.core_rradd = rraddMethods[i][PastixComplex64-2];

                    if ( params.output == stdout ) {
                        fprintf( params.output, "%7s %5d %4d %e %6d %5d %e %e",
                                 compmeth_shnames[i], p, mode, params.tol_gen,
                                 A.n, A.rk, A.norm, params.tol_cmp );
                    }
                    else {
                        fprintf( params.output, "%s;%d;%d;%e;%d;%d;%e;%e;",
                                 compmeth_shnames[i], p, mode, params.tol_gen,
                                 A.n, A.rk, A.norm, params.tol_cmp );
                    }

                    ret = z_lowrank_ge2lr_performance( params.output, params.nb_runs,
                                                       &lowrank, &A );
                    rc += (ret ? 1 : 0 );
                }
            }
        }
        free(A.fr);
    }

    if ( params.output != stdout ) {
        fclose( params.output );
    }

    if( rc == 0 ) {
        printf( " -- All tests PASSED --\n" );
        return EXIT_SUCCESS;
    }
    else
    {
        printf( " -- %d tests FAILED --\n", rc );
        return EXIT_FAILURE;
    }
}
