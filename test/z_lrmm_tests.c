/**
 *
 * @file z_lrmm_tests.c
 *
 * Tests and validate the core_zlrmm() routine.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Gregoire Pichon
 * @date 2019-11-12
 *
 * @precisions normal z -> c d s
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
#include <cblas.h>
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"
#include "kernels/pastix_zlrcores.h"
#include "z_tests.h"

#define PRINT_RES(_ret_)                        \
    if(_ret_ == -1) {                           \
        printf("UNDEFINED\n");                  \
    }                                           \
    else if(_ret_ > 0) {                        \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

int main( int argc, char **argv )
{
    pastix_int_t        m, n, k, Cm, Cn, offx, offy, use_reltol;
    double              eps = LAPACKE_dlamch_work('e');
    double              tolerance = sqrt(eps);
    double              threshold = tolerance * tolerance;
    test_matrix_t       A, B, C;
    pastix_complex64_t *Cfr;
    pastix_complex64_t  alpha = -3.48;
    pastix_complex64_t  beta  = 1.;
    pastix_lr_t         lowrank = {
        .compress_when       = PastixCompressWhenEnd,
        .compress_method     = PastixCompressMethodPQRCP,
        .compress_min_width  = 0,
        .compress_min_height = 0,
        .use_reltol          = 0,
        .tolerance           = tolerance,
        .core_ge2lr          = core_zge2lr_svd,
        .core_rradd          = core_zrradd_svd,
    };
    int ranks[3], r, s, i, j, l, meth;
    int err  = 0;
    int ret  = 0;
    int mode = 0;

    core_get_rklimit = core_get_rklimit_max;

    for(use_reltol=0; use_reltol<2; use_reltol++) {
        lowrank.use_reltol = use_reltol;
        for (s=100; s<=200; s = 2*s) {
            ranks[0] = s + 1;
            ranks[1] = 16;
            ranks[2] = 2;

            m = s / 2;
            n = s / 4;
            k = s / 6;

            offx = 1;
            offy = 2;

            Cm = s;
            Cn = s;

            /* Matrix A */
            for (i=0; i<3; i++) {
                A.m  = m;
                A.n  = k;
                A.ld = m;
                A.rk = pastix_imin( m, k ) / ranks[i];
                z_lowrank_genmat_comp( &lowrank, mode, threshold, &A );

                /* Matrix B */
                for (j=0; j<3; j++) {
                    B.m  = n;
                    B.n  = k;
                    B.ld = n;
                    B.rk = pastix_imin( n, k ) / ranks[j];
                    z_lowrank_genmat_comp( &lowrank, mode, threshold, &B );

                    /* Matrix C */
                    for (l=0; l<3; l++) {
                        C.m  = Cm;
                        C.n  = Cn;
                        C.ld = Cm;
                        C.rk = pastix_imin( Cm, Cn ) / ranks[l];
                        z_lowrank_genmat_comp( &lowrank, mode, threshold, &C );

                        printf( "  -- Test LRMM Cm=%ld, Cn=%ld, m=%ld, n=%ld, k=%ld, rA=%ld, rB=%ld, rC=%ld (%s)\n",
                                (long)Cm, (long)Cn, (long)m, (long)n, (long)k,
                                (long)(A.lr.rk), (long)(B.lr.rk), (long)(C.lr.rk),
                                use_reltol ? "Relative": "Absolute" );

                        /* Compute the full rank GEMM */
                        Cfr = C.fr;
                        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasConjTrans,
                                     m, n, k,
                                     CBLAS_SADDR( alpha ), A.fr, A.ld,
                                                           B.fr, B.ld,
                                     CBLAS_SADDR( beta ),  Cfr + C.ld * offy + offx, C.ld );

                        /* Update the norm of C with the norm of the result */
                        C.norm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', Cm, Cn,
                                                      Cfr, C.ld, NULL );

                        fprintf( stdout, "%7s %4s %12s %12s %12s %12s\n",
                                 "Method", "Rank", "Time", "||C||_f", "||c(C)-C||_f",
                                 "||c(C)-C||_f/(||C||_f * eps)" );

                        ret = 0;

                        /* Let's test all methods we have */
                        for(meth=0; meth<PastixCompressMethodNbr; meth++)
                        {
                            lowrank.compress_method = meth;
                            lowrank.core_ge2lr = ge2lrMethods[meth][PastixComplex64-2];
                            lowrank.core_rradd = rraddMethods[meth][PastixComplex64-2];

                            r = z_lowrank_check_lrmm( &lowrank, offx, offy,
                                                      alpha, &A, &B, beta, &C );
                            ret += r * (1 << meth);
                        }
                        lowrank.compress_method = PastixCompressMethodSVD;
                        lowrank.core_ge2lr = ge2lrMethods[PastixCompressMethodSVD][PastixComplex64-2];
                        lowrank.core_rradd = rraddMethods[PastixCompressMethodSVD][PastixComplex64-2];

                        core_zlrfree( &(C.lr) );
                        free( C.fr );
                        PRINT_RES( ret );
                    }

                    core_zlrfree( &(B.lr) );
                    free( B.fr );
                }
                core_zlrfree( &(A.lr) );
                free( A.fr );
            }
        }
    }

    if( err == 0 ) {
        printf(" -- All tests PASSED --\n");
        return EXIT_SUCCESS;
    }
    else
    {
        printf(" -- %d tests FAILED --\n", err);
        return EXIT_FAILURE;
    }

    (void) argc;
    (void) argv;
}
