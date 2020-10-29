/**
 *
 * @file sequential_ztrsm.c
 *
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2020-01-26
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include "bcsc/bcsc.h"
#include "bcsc/bcsc_z.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_STARPU)
#include "starpu/pastix_zstarpu.h"
#endif

#if defined(PASTIX_WITH_MPI)
#include "coeftab.h"
#endif

void
sequential_ztrsm( pastix_data_t *pastix_data, int side, int uplo, int trans, int diag,
                  sopalin_data_t *sopalin_data,
                  int nrhs, pastix_complex64_t *b, int ldb )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk   *cblk;
    pastix_int_t i, cblknbr;

    pastix_solv_mode_t mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];

    /* Backward like */
    if ( ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
         ( (side == PastixLeft)  && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
         ( (side == PastixRight) && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
         ( (side == PastixRight) && (uplo == PastixLower) && (trans == PastixNoTrans) ) )
    {
        cblknbr = (mode == PastixSolvModeLocal) ? datacode->cblkschur : datacode->cblknbr;

        cblk = datacode->cblktab + cblknbr - 1;
        for (i=0; i<cblknbr; i++, cblk--){
            if( cblk->cblktype & CBLK_RECV ){
                cpucblk_zsend_rhs_backward( datacode, cblk, b );
                continue;
            }

            if( cblk->cblktype & CBLK_FANIN ){
                cpucblk_zrecv_rhs_backward( datacode, cblk, b );
            }

            solve_cblk_ztrsmsp_backward( mode, side, uplo, trans, diag,
                                         datacode, cblk, nrhs, b, ldb );
        }

        bvec_znullify_remote( pastix_data, b );
        bvec_zallreduce( pastix_data, b );
    }
    /* Forward like */
    else
        /**
         * ( (side == PastixRight) && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
         * ( (side == PastixRight) && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
         * ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
         * ( (side == PastixLeft)  && (uplo == PastixLower) && (trans == PastixNoTrans) )
         */
    {
        pastix_complex64_t *work;
        MALLOC_INTERN( work, datacode->colmax, pastix_complex64_t );

        bvec_znullify_remote( pastix_data, b );

        cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;
        cblk = datacode->cblktab;
        for (i=0; i<cblknbr; i++, cblk++){
            if( cblk->cblktype & CBLK_FANIN ){
                cpucblk_zsend_rhs_forward( datacode, cblk, b );
                continue;
            }

            if( cblk->cblktype & CBLK_RECV ) {
                cpucblk_zrecv_rhs_forward( datacode, cblk, work, nrhs, b, ldb );
                continue;
            }

            solve_cblk_ztrsmsp_forward( mode, side, uplo, trans, diag,
                                        datacode, cblk, nrhs, b, ldb );
        }

        memFree_null(work);
    }
}

struct args_ztrsm_t
{
    pastix_data_t  *pastix_data;
    int side, uplo, trans, diag;
    sopalin_data_t *sopalin_data;
    int nrhs;
    pastix_complex64_t *b;
    int ldb;
};

void
thread_ztrsm_static( isched_thread_t *ctx, void *args )
{
    struct args_ztrsm_t *arg = (struct args_ztrsm_t*)args;
    pastix_data_t      *pastix_data  = arg->pastix_data;
    sopalin_data_t     *sopalin_data = arg->sopalin_data;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    pastix_complex64_t *b = arg->b;
    int side  = arg->side;
    int uplo  = arg->uplo;
    int trans = arg->trans;
    int diag  = arg->diag;
    int nrhs  = arg->nrhs;
    int ldb   = arg->ldb;
    SolverCblk *cblk;
    Task       *t;
    pastix_int_t i, ii;
    pastix_int_t tasknbr, *tasktab;
    pastix_solv_mode_t mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];
    int rank = ctx->rank;

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    /* Backward like */
    if ( ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
         ( (side == PastixLeft)  && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
         ( (side == PastixRight) && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
         ( (side == PastixRight) && (uplo == PastixLower) && (trans == PastixNoTrans) ) )
    {
        /* Init ctrbcnt in parallel */
        for (ii=0; ii<tasknbr; ii++) {
            i = tasktab[ii];
            t = datacode->tasktab + i;
            cblk = datacode->cblktab + t->cblknum;

            if ( (cblk->cblktype & CBLK_IN_SCHUR) && (mode != PastixSolvModeSchur) ) {
                cblk->ctrbcnt = 0;
            }
            else {
                cblk->ctrbcnt = cblk[1].fblokptr - cblk[0].fblokptr - 1;
            }
        }
        isched_barrier_wait( &(ctx->global_ctx->barrier) );

        for (ii=tasknbr-1; ii>=0; ii--) {
            i = tasktab[ii];
            t = datacode->tasktab + i;
            cblk = datacode->cblktab + t->cblknum;

            /* Wait */
            do { } while( cblk->ctrbcnt );

            solve_cblk_ztrsmsp_backward( mode, side, uplo, trans, diag,
                                         datacode, cblk, nrhs, b, ldb );
        }
    }
    /* Forward like */
    else
        /**
         * ( (side == PastixRight) && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
         * ( (side == PastixRight) && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
         * ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
         * ( (side == PastixLeft)  && (uplo == PastixLower) && (trans == PastixNoTrans) )
         */
    {
        /* Init ctrbcnt in parallel */
        for (ii=0; ii<tasknbr; ii++) {
            i = tasktab[ii];
            t = datacode->tasktab + i;
            cblk = datacode->cblktab + t->cblknum;
            cblk->ctrbcnt = cblk[1].brownum - cblk[0].brownum;
        }
        isched_barrier_wait( &(ctx->global_ctx->barrier) );

        for (ii=0; ii<tasknbr; ii++) {
            i = tasktab[ii];
            t = datacode->tasktab + i;
            cblk = datacode->cblktab + t->cblknum;

            if ( (cblk->cblktype & CBLK_IN_SCHUR) && (mode != PastixSolvModeSchur) )
                continue;

            /* Wait */
            do { } while( cblk->ctrbcnt );

            solve_cblk_ztrsmsp_forward( mode, side, uplo, trans, diag,
                                        datacode, cblk, nrhs, b, ldb );
        }
    }
}

void
static_ztrsm( pastix_data_t *pastix_data, int side, int uplo, int trans, int diag,
              sopalin_data_t *sopalin_data,
              int nrhs, pastix_complex64_t *b, int ldb )
{
    struct args_ztrsm_t args_ztrsm = {pastix_data, side, uplo, trans, diag, sopalin_data, nrhs, b, ldb};
    isched_parallel_call( pastix_data->isched, thread_ztrsm_static, &args_ztrsm );
}

#if defined(PASTIX_WITH_MPI)
void
runtime_ztrsm( pastix_data_t *pastix_data, int side, int uplo, int trans, int diag,
               sopalin_data_t *sopalin_data,
               int nrhs, pastix_complex64_t *b, int ldb )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk *cblk;
    pastix_int_t i, cblknbr;

    /* Collect the matrix on node 0 */
    coeftab_gather( datacode, datacode->solv_comm, 0, PastixComplex64 );

    if ( sopalin_data->solvmtx->clustnum == 0 ) {
        pastix_solv_mode_t mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];

        /* Backward like */
        if ( ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
             ( (side == PastixLeft)  && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
             ( (side == PastixRight) && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
             ( (side == PastixRight) && (uplo == PastixLower) && (trans == PastixNoTrans) ) )
        {
            cblknbr = (mode == PastixSolvModeLocal) ? datacode->cblkschur : datacode->cblknbr;

            cblk = datacode->cblktab + cblknbr - 1;
            for ( i=0; i<cblknbr; i++, cblk-- ) {
                assert( !(cblk->cblktype & (CBLK_FANIN | CBLK_RECV)) );
                solve_cblk_ztrsmsp_backward( mode, side, uplo, trans, diag,
                                             datacode, cblk, nrhs, b, ldb );
            }
        }
        /* Forward like */
        else
            /**
             * ( (side == PastixRight) && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
             * ( (side == PastixRight) && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
             * ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
             * ( (side == PastixLeft)  && (uplo == PastixLower) && (trans == PastixNoTrans) )
             */
        {
            cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;
            cblk = datacode->cblktab;
            for (i=0; i<cblknbr; i++, cblk++){
                solve_cblk_ztrsmsp_forward( mode, side, uplo, trans, diag,
                                            datacode, cblk, nrhs, b, ldb );
            }
        }

        /* Free the gathered coefficients of the matrix */
        coeftab_nullify( datacode );
    }
    else {
        memset( b, 0, ldb * nrhs * sizeof(pastix_complex64_t) );
    }

    bvec_zallreduce( pastix_data, b );
}
#endif

static void (*ztrsm_table[5])(pastix_data_t *, int, int, int, int, sopalin_data_t *,
                              int, pastix_complex64_t *, int) = {
    sequential_ztrsm,
    static_ztrsm,
#if defined(PASTIX_WITH_PARSEC)
    NULL, /* parsec_ztrsm not yet implemented */
#else
    NULL,
#endif
#if defined(PASTIX_WITH_STARPU)
    starpu_ztrsm,
#else
    NULL,
#endif
    static_ztrsm
};

void
sopalin_ztrsm( pastix_data_t *pastix_data, int side, int uplo, int trans, int diag,
               sopalin_data_t *sopalin_data,
               int nrhs, pastix_complex64_t *b, int ldb )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*ztrsm)(pastix_data_t *, int, int, int, int, sopalin_data_t *,
                  int, pastix_complex64_t *, int) = ztrsm_table[ sched ];

    if (ztrsm == NULL) {
        ztrsm = static_ztrsm;
    }

#if defined(PASTIX_WITH_MPI)
    if( pastix_data->inter_node_procnbr > 1 ) {
        if( (sched == PastixSchedStarPU) || (sched == PastixSchedParsec) ) {
            ztrsm = runtime_ztrsm;
        }
        else {
            /* Force sequential if MPI as no other runtime is supported yet */
            ztrsm = sequential_ztrsm;
        }
    }
#endif

    ztrsm( pastix_data, side, uplo, trans, diag, sopalin_data, nrhs, b, ldb );

#if defined(PASTIX_WITH_MPI)
   MPI_Barrier( pastix_data->inter_node_comm );
#endif
}
