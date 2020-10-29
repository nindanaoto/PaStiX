/**
 *
 * @file z_ge2lr_tests.c
 *
 * Tests and validate the Xge2lr routine.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Gregoire Pichon
 * @date 2019-11-12
 *
 * @precisions normal z -> z c d s
 *
 **/
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
#include "z_tests.h"
#include "tests.h"
#include "kernels/pastix_zlrcores.h"

int main( int argc, char **argv )
{
    test_matrix_t A;
    pastix_int_t n;
    int mode, p, i, ret, rc = 0;
    test_param_t params;
    double eps = LAPACKE_dlamch_work('e');
    pastix_lr_t lowrank;

    testGetOptions( argc, argv, &params, eps );

    fprintf( stdout, "%7s %4s %12s %12s %12s %12s\n",
             "Method", "Rank", "Time", "||A||_f", "||A-UVt||_f",
             "||A-UVt||_f/(||A||_f * eps)" );

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
                printf( "   -- Test GE2LR TolGen=%e TolCmp=%e M=N=LDA=%ld R=%ld MODE=%d\n",
                        params.tol_gen, lowrank.tolerance, (long)A.n, (long)A.rk, mode );

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

                    ret = z_lowrank_check_ge2lr( &lowrank, &A );
                    rc += (ret ? 1 : 0 );
                }
            }
        }
        free(A.fr);
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
