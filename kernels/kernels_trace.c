/**
 *
 * @file kernels_trace.c
 *
 * @copyright 2004-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * PaStiX trace and modelling routines
 *
 * @version 6.0.3
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2019-12-19
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "common.h"
#include "bcsc.h"
#include "solver.h"
#include "kernels_trace.h"

/**
 * @brief Compute the maximal rank accepted for a given matrix size. The pointer is set according to the low-rank strategy used.
 * @param[in] M The number of rows of the matrix
 * @param[in] N The number of columns of the matrix
 * @return The maximal rank accepted for this matrix size.
 */
pastix_int_t (*core_get_rklimit)( pastix_int_t, pastix_int_t ) = core_get_rklimit_end;

volatile double  kernels_flops[PastixKernelLvl1Nbr];

volatile int32_t kernels_trace_started = 0;

#if defined(PASTIX_WITH_EZTRACE)

int pastix_eztrace_level = 1;

#endif

#if defined(PASTIX_GENERATE_MODEL)

pastix_model_entry_t *model_entries     = NULL;
volatile int32_t      model_entries_nbr = -1;
int32_t               model_size        = 0;

#endif

pastix_atomic_lock_t lock_flops = PASTIX_ATOMIC_UNLOCKED;
double overall_flops[3] = { 0.0, 0.0, 0.0 };

/**
 *******************************************************************************
 *
 * @brief Start the trace module
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure of the problem to give input information
 *          to the different trace modes.
 *
 *******************************************************************************/
void
kernelsTraceStart( const pastix_data_t *pastix_data )
{
    const SolverMatrix *solvmtx = pastix_data->solvmatr;
    int32_t nbstart;

    pastix_atomic_lock( &lock_flops );
    nbstart = pastix_atomic_inc_32b( &(kernels_trace_started) );
    if ( nbstart > 1 ) {
        pastix_atomic_unlock( &lock_flops );
        return;
    }

#if defined(PASTIX_WITH_EZTRACE)
    {
        char *level = pastix_getenv("PASTIX_EZTRACE_LEVEL");
        if (level != NULL) {
            pastix_eztrace_level = atoi(level);
            pastix_cleanenv(level);
        }

        if ( pastix_data->dir_local != NULL ) {
            pastix_setenv( "EZTRACE_TRACE_DIR", pastix_data->dir_local, 1 );
        }
        eztrace_start ();
    }
#endif /* defined(PASTIX_WITH_EZTRACE) */

#if defined(PASTIX_GENERATE_MODEL)
    {
        pastix_int_t cblknbr   = solvmtx->cblknbr;
        pastix_int_t cblkmin2d = solvmtx->cblkmin2d;
        pastix_int_t total_number_of_tasks = 0;
        pastix_int_t nbfact, nbtrsm, nbgemm;
        pastix_int_t cblknum;
        SolverCblk  *cblk;

        /* Factorization kernels */
        nbfact = cblknbr;

        /* TRSM kernels */
        nbtrsm = cblkmin2d + (cblknbr - cblkmin2d) * solvmtx->cblkmaxblk;
        if ( solvmtx->factotype == PastixFactLU ) {
            nbtrsm *= 2;
        }

        /* GEMM kernels */
        nbgemm = solvmtx->bloknbr - cblknbr;
        if ( solvmtx->factotype == PastixFactLU ) {
            nbgemm *= 2;
        }

        cblk = solvmtx->cblktab+cblkmin2d;
        for(cblknum = cblkmin2d; cblknum < cblknbr; cblknum++, cblk++ ) {
            pastix_int_t nbodb = (cblk[1].fblokptr - cblk[0].fblokptr) - 1;

            if ( solvmtx->factotype == PastixFactLU ) {
                nbgemm += nbodb * nbodb;
            }
            else {
                nbgemm += (nbodb * (nbodb-1)) / 2;
            }
        }

        total_number_of_tasks = nbfact + nbtrsm + nbgemm;
        model_entries = malloc( total_number_of_tasks * sizeof(pastix_model_entry_t) );
        model_size = total_number_of_tasks;
    }
#endif

    memset( (void*)kernels_flops, 0, PastixKernelLvl1Nbr * sizeof(double) );

    overall_flops[0] = 0.0;
    overall_flops[1] = 0.0;
    overall_flops[2] = 0.0;
    kernels_trace_started = 1;

    (void)solvmtx;
    pastix_atomic_unlock( &lock_flops );
    return;
}

/**
 *******************************************************************************
 *
 * @brief Stop the trace module
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure of the problem to get input information
 *          for the different trace modes, and store output statistics.
 *
 *******************************************************************************/
double
kernelsTraceStop( const pastix_data_t *pastix_data )
{
    double total_flops = 0.0;
    int32_t nbstart;

    assert( kernels_trace_started > 0 );
    pastix_atomic_lock( &lock_flops );
    nbstart = pastix_atomic_dec_32b( &(kernels_trace_started) );
    if ( nbstart > 0 ) {
        pastix_atomic_unlock( &lock_flops );
        return total_flops;
    }

#if defined(PASTIX_WITH_EZTRACE)
    eztrace_stop ();
#endif

#if defined(PASTIX_GENERATE_MODEL)
    {
        char *prec_names[4] = {
            "s - single real", "d - double real",
            "c - single complex", "z - double complex"
        };
        pastix_model_entry_t *entry = model_entries;
        pastix_int_t i, gpucase;
        FILE *f;

        f = fopen( "model.csv", "w" );
        if ( f == NULL ) {
            goto end_model;
        }

        gpucase = pastix_data->iparm[IPARM_GPU_NBR];
        if ( gpucase ) {
            fprintf(f, "# GPU Model data\n");
        }
        else {
            fprintf(f, "# CPU Model data\n");
        }

        fprintf( f, "# Precision: %d - %s\n", pastix_data->bcsc->flttype - 2, prec_names[ pastix_data->bcsc->flttype - 2 ] );
        fprintf( f, "Kernel;M;N;K;Time\n" );

        for(i=0; i <= model_entries_nbr; i++, entry++ ) {
            switch( entry->ktype ) {
            case PastixKernelGETRF:        pastix_attr_fallthrough;
            case PastixKernelHETRF:        pastix_attr_fallthrough;
            case PastixKernelPOTRF:        pastix_attr_fallthrough;
            case PastixKernelPXTRF:        pastix_attr_fallthrough;
            case PastixKernelSYTRF:        pastix_attr_fallthrough;
            case PastixKernelSCALOCblk:    pastix_attr_fallthrough;
            case PastixKernelSCALOBlok:    pastix_attr_fallthrough;
            case PastixKernelTRSMCblk1d:   pastix_attr_fallthrough;
            case PastixKernelTRSMCblk2d:   pastix_attr_fallthrough;
            case PastixKernelTRSMCblkLR:   pastix_attr_fallthrough;
            case PastixKernelTRSMBlokLR:   pastix_attr_fallthrough;
            case PastixKernelGEMMCblk1d1d: pastix_attr_fallthrough;
            case PastixKernelGEMMCblkFRLR: pastix_attr_fallthrough;
            case PastixKernelGEMMCblkLRLR: pastix_attr_fallthrough;
            case PastixKernelGEMMBlokLRLR:
                if ( gpucase ) {
                    continue;
                }

                pastix_attr_fallthrough;
            default:
                fprintf( f, "%d;%d;%d;%d;%e\n",
                         entry->ktype, entry->m, entry->n, entry->k, entry->time );
            }
        }

        fclose( f );

        free( model_entries );

        /* Reinitialize values */
        model_entries     = NULL;
        model_entries_nbr = -1;
        model_size        = 0;
    }
  end_model:
#endif

    /* Update the real number of Flops performed */
    pastix_data->dparm[DPARM_FACT_RLFLOPS] = overall_flops[0] + overall_flops[1] + overall_flops[2];

#if defined(PASTIX_SUPERNODE_STATS)
    if (pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot) {
        fprintf( stdout,
                 "    Details of the number of operations:\n"
                 "      - POTRF(A11) + TRSM(A11, A21): %6.2lf %cFlops\n"
                 "      - HERK(A21, A22)             : %6.2lf %cFlops\n"
                 "      - POTRF(A22)                 : %6.2lf %cFlops\n"
                 "      Total                        : %6.2lf %cFlops\n",
                 pastix_print_value( overall_flops[0] ), pastix_print_unit( overall_flops[0] ),
                 pastix_print_value( overall_flops[1] ), pastix_print_unit( overall_flops[1] ),
                 pastix_print_value( overall_flops[2] ), pastix_print_unit( overall_flops[2] ),
                 pastix_print_value( pastix_data->dparm[DPARM_FACT_RLFLOPS] ),
                 pastix_print_unit( pastix_data->dparm[DPARM_FACT_RLFLOPS] ) );
    }
#endif /* defined(PASTIX_SUPERNODE_STATS) */

    kernels_trace_started = 0;
    pastix_atomic_unlock( &lock_flops );
    (void)pastix_data;
    return total_flops;
}
