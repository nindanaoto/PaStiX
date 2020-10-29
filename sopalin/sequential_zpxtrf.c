/**
 *
 * @file sequential_zpxtrf.c
 *
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2020-01-29
 *
 * @precisions normal z -> z c
 *
 **/
#include "common.h"
#include "isched.h"
#include "solver.h"
#include "sopalin_data.h"
#include "sopalin/coeftab_z.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_PARSEC)
#include "parsec/pastix_zparsec.h"
#endif

#if defined(PASTIX_WITH_STARPU)
#include "starpu/pastix_zstarpu.h"
#endif

void
sequential_zpxtrf( pastix_data_t  *pastix_data,
                   sopalin_data_t *sopalin_data )
{
    SolverMatrix       *datacode = pastix_data->solvmatr;
    SolverCblk         *cblk;
    pastix_complex64_t *work;
    pastix_int_t  i, lwork;
    (void)sopalin_data;

    lwork = datacode->gemmmax;
    if ( datacode->lowrank.compress_when == PastixCompressWhenBegin ) {
        lwork = pastix_imax( lwork, 2 * datacode->blokmax );
    }
    MALLOC_INTERN( work, lwork, pastix_complex64_t );

    cblk = datacode->cblktab;
    for (i=0; i<datacode->cblknbr; i++, cblk++){
        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            break;
        }

        /* Wait for incoming dependencies */
        if ( cpucblk_zincoming_deps( 0, PastixLCoef,
                                     datacode, cblk ) )
        {
            continue;
        }

        /* Compute */
        cpucblk_zpxtrfsp1d( datacode, cblk,
                            work, lwork );
    }

    memFree_null( work );
}

void
thread_zpxtrf_static( isched_thread_t *ctx, void *args )
{
    sopalin_data_t     *sopalin_data = (sopalin_data_t*)args;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    Task               *t;
    pastix_complex64_t *work;
    pastix_int_t i, ii, lwork;
    pastix_int_t tasknbr, *tasktab;
    int rank = ctx->rank;

    lwork = datacode->gemmmax;
    if ( datacode->lowrank.compress_when == PastixCompressWhenBegin ) {
        lwork = pastix_imax( lwork, 2 * datacode->blokmax );
    }
    MALLOC_INTERN( work, lwork, pastix_complex64_t );

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            continue;
        }

        /* Wait for incoming dependencies */
        if ( cpucblk_zincoming_deps( rank, PastixLCoef,
                                     datacode, cblk ) )
        {
            continue;
        }

        /* Compute */
        cpucblk_zpxtrfsp1d( datacode, cblk,
                            work, lwork );
    }

    memFree_null( work );
}

void
static_zpxtrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    isched_parallel_call( pastix_data->isched, thread_zpxtrf_static, sopalin_data );
}

struct args_zpxtrf_t
{
    sopalin_data_t     *sopalin_data;
    volatile int32_t    taskcnt;
};

void
thread_zpxtrf_dynamic( isched_thread_t *ctx, void *args )
{
    struct args_zpxtrf_t *arg = (struct args_zpxtrf_t*)args;
    sopalin_data_t       *sopalin_data = arg->sopalin_data;
    SolverMatrix         *datacode = sopalin_data->solvmtx;
    SolverCblk           *cblk;
    Task                 *t;
    pastix_queue_t       *computeQueue;
    pastix_complex64_t   *work;
    pastix_int_t          i, ii, lwork;
    pastix_int_t          tasknbr, *tasktab, cblknum;
    int32_t               local_taskcnt = 0;
    int                   rank = ctx->rank;
    int                   dest = (ctx->rank + 1)%ctx->global_ctx->world_size;

    lwork = datacode->gemmmax;
    if ( datacode->lowrank.compress_when == PastixCompressWhenBegin ) {
        lwork = pastix_imax( lwork, 2 * datacode->blokmax );
    }
    MALLOC_INTERN( work, lwork, pastix_complex64_t );
    MALLOC_INTERN( datacode->computeQueue[rank], 1, pastix_queue_t );

    tasknbr      = datacode->ttsknbr[rank];
    tasktab      = datacode->ttsktab[rank];
    computeQueue = datacode->computeQueue[rank];
    pqueueInit( computeQueue, tasknbr );

    /* Initialize the local task queue with available cblks */
    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;

        if ( !(t->ctrbcnt) ) {
            pqueuePush1( computeQueue, t->cblknum, t->prionum );
        }
    }

    /* Make sure that all computeQueues are allocated */
    isched_barrier_wait( &(ctx->global_ctx->barrier) );

    while( arg->taskcnt > 0 )
    {
        cblknum = pqueuePop(computeQueue);

#if defined(PASTIX_WITH_MPI)
        /* Nothing to do, let's make progress on comunications */
        if( cblknum == -1 ) {
            cpucblk_zmpi_progress( PastixLCoef, datacode, rank );
            cblknum = pqueuePop(computeQueue);
        }
#endif

        /* No more local job, let's steal our neighbours */
        if( cblknum == -1 ) {
            if ( local_taskcnt ) {
                pastix_atomic_sub_32b( &(arg->taskcnt), local_taskcnt );
                local_taskcnt = 0;
            }
            cblknum = stealQueue( datacode, rank, &dest,
                                  ctx->global_ctx->world_size );
        }

        /* Still no job, let's loop again */
        if ( cblknum == -1 ) {
            continue;
        }

        cblk = datacode->cblktab + cblknum;
        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            continue;
        }
        cblk->threadid = rank;

        /* Compute */
        cpucblk_zpxtrfsp1d( datacode, cblk,
                            work, lwork );
        local_taskcnt++;
    }
    memFree_null( work );

    /* Make sure that everyone is done before freeing */
    isched_barrier_wait( &(ctx->global_ctx->barrier) );
    pqueueExit( computeQueue );
    memFree_null( computeQueue );
}

void
dynamic_zpxtrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    SolverMatrix        *datacode = sopalin_data->solvmtx;
    int32_t              taskcnt = datacode->tasknbr;
    struct args_zpxtrf_t args_zpxtrf = { sopalin_data, taskcnt };

    /* Allocate the computeQueue */
    MALLOC_INTERN( datacode->computeQueue,
                   pastix_data->isched->world_size, pastix_queue_t * );

    isched_parallel_call( pastix_data->isched, thread_zpxtrf_dynamic, &args_zpxtrf );

    memFree_null( datacode->computeQueue );

#if defined(PASTIX_WITH_MPI)
    MPI_Barrier( pastix_data->inter_node_comm );
#endif
}

static void (*zpxtrf_table[5])(pastix_data_t *, sopalin_data_t *) = {
    sequential_zpxtrf,
    static_zpxtrf,
#if defined(PASTIX_WITH_PARSEC)
    parsec_zpxtrf,
#else
    NULL,
#endif
#if defined(PASTIX_WITH_STARPU)
    starpu_zpxtrf,
#else
    NULL,
#endif
    dynamic_zpxtrf
};

void
sopalin_zpxtrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zpxtrf)(pastix_data_t *, sopalin_data_t *) = zpxtrf_table[ sched ];

    if (zpxtrf == NULL) {
        sched = PastixSchedSequential;
        zpxtrf = static_zpxtrf;
    }

    if ( (sched == PastixSchedSequential) ||
         (sched == PastixSchedStatic)     ||
         (sched == PastixSchedDynamic) )
    {
        solverRequestInit( sopalin_data->solvmtx );
        solverRecvInit( PastixLCoef, sopalin_data->solvmtx, PastixComplex64 );
    }

    zpxtrf( pastix_data, sopalin_data );

    if ( (sched == PastixSchedSequential) ||
         (sched == PastixSchedStatic)     ||
         (sched == PastixSchedDynamic) )
    {
        cpucblk_zrequest_cleanup( PastixLCoef, sched, sopalin_data->solvmtx );
        solverRequestExit( sopalin_data->solvmtx );
        solverRecvExit( sopalin_data->solvmtx );
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( pastix_data, sopalin_data->solvmtx, "pxtrf" );
#endif
}
