/**
 *
 * @file bvec_gemv_tests.c
 *
 * Tests for preformance of gemv in parallel and sequential versions.
 * Times is computed on the average of 50 calls to gemv function.
 * Size of the matrice is given by the first two dimensions of the Laplacian
 * driver.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Vincent Bridonneau
 * @date 2020-01-29
 *
 **/
#include <pastix.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <cblas.h>
#include <lapacke.h>
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */

#include "common/common.h"
#include "common/flops.h"
#include "bcsc/bcsc.h"
#include "bcsc/bcsc_z.h"
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"
#include "refinement/z_refine_functions.h"

#include "z_tests.h"
#include "c_tests.h"
#include "d_tests.h"
#include "s_tests.h"

#define PRINT_RES(_ret_)                        \
    if(_ret_) {                                 \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* schednames[] = { "Sequential", "Static" };

int main ( int argc, char **argv )
{
    pastix_data_t   *pastix_data = NULL;
    pastix_int_t     iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    pastix_fixdbl_t  dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    spm_driver_t     driver = SpmDriverRSA;
    char            *filename = NULL;
    int              check = 1;
    pastix_int_t     m, n, s, t;
    pastix_int_t     dim3, dof;
    spm_coeftype_t   flttype;
    pastix_fixdbl_t  alpha, beta;
    int err = 0;
    int rc = 0;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );
    iparm[IPARM_VERBOSE] = 0;

    /**
     * Get options from command line
     * Prevent from failing if no arguments is given
     */
    if ( argc > 1 ) {
        pastixGetOptions( argc, argv,
                          iparm, dparm,
                          &check, &driver, &filename );
    }

    /**
     * Startup pastix to start the scheduler
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    if ( driver != SpmDriverLaplacian ) {
        flttype = SpmDouble;
        m = 10000;
        n = iparm[IPARM_GMRES_IM];
    }
    else {
        spmParseLaplacianInfo( filename, &flttype, &m, &n, &dim3, &alpha, &beta, &dof );
    }

    if ( check ) {
        for( s=PastixSchedSequential; s<=PastixSchedStatic; s++ )
        {
            iparm[IPARM_SCHEDULER] = s;
            for( t=SpmFloat; t<=SpmComplex64; t++ )
            {
                printf("   Case %s - %s: ",
                       fltnames[t], schednames[s] );
                switch( t ){
                case SpmComplex64:
                    rc = z_bvec_gemv_check( pastix_data, check, m, n );
                    break;

                case SpmComplex32:
                    rc = c_bvec_gemv_check( pastix_data, check, m, n );
                    break;

                case SpmFloat:
                    rc = s_bvec_gemv_check( pastix_data, check, m, n );
                    break;

                case SpmDouble:
                default:
                    rc = d_bvec_gemv_check( pastix_data, check, m, n );
                }
                PRINT_RES(rc);
            }
        }
        if( err == 0 ) {
            printf(" -- All tests PASSED --\n");
            rc = EXIT_SUCCESS;
        }
        else
        {
            printf(" -- %d tests FAILED --\n", err);
            rc = EXIT_FAILURE;
        }
    }
    else {
        switch( flttype ){
        case SpmComplex64:
            rc = z_bvec_gemv_check( pastix_data, check, m, n );
            break;

        case SpmComplex32:
            rc = c_bvec_gemv_check( pastix_data, check, m, n );
            break;

        case SpmFloat:
            rc = s_bvec_gemv_check( pastix_data, check, m, n );
            break;

        case SpmDouble:
        default:
            rc = d_bvec_gemv_check( pastix_data, check, m, n );
        }
    }

    pastixFinalize( &pastix_data );

    return rc;
}
