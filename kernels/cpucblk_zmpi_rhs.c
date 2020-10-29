/**
 *
 * @file cpucblk_zmpi_rhs.c
 *
 * Precision dependent routines to manag communications for the solve part.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.1.0
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2020-01-29
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

/**
 *******************************************************************************
 *
 * @brief Send the rhs associated to a cblk->lcolidx to the remote node.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix holding the communicator.
 *
 * @param[in] cblk
 *          The cblk which defines the part to sent.
 *
 * @param[in] b
 *          The rhs that will be sent to the cblk->ownerid
 *
 *******************************************************************************/
void
cpucblk_zsend_rhs_forward( const SolverMatrix *solvmtx,
                           SolverCblk         *cblk,
                           pastix_complex64_t *b )
{
#if defined(PASTIX_WITH_MPI)
    pastix_int_t colnbr = cblk_colnbr(cblk);
    int rc;

    assert( colnbr <= solvmtx->colmax );
    assert( cblk->cblktype & CBLK_FANIN );

#if defined (PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Fwd: Send cblk %ld to %2d at index %ld of size %ld\n",
             solvmtx->clustnum, (long)cblk->gcblknum, cblk->ownerid,
             (long)cblk->lcolidx, (long)colnbr );
#endif

    rc = MPI_Send( b + cblk->lcolidx, colnbr, PASTIX_MPI_COMPLEX64,
                   cblk->ownerid, cblk->gcblknum, solvmtx->solv_comm );
    assert( rc == MPI_SUCCESS );

    (void)rc;
#else
    (void)solvmtx;
    (void)cblk;
    (void)b;
#endif
}

/**
 *******************************************************************************
 *
 * @brief Send the rhs associated to a cblk->lcolidx to the remote node.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix holding the communicator.
 *
 * @param[in] cblk
 *          The cblk which defines the part to sent.
 *
 * @param[in] b
 *          The rhs that will be sent to the cblk->ownerid
 *
 *******************************************************************************/
void
cpucblk_zsend_rhs_backward( const SolverMatrix *solvmtx,
                            SolverCblk         *cblk,
                            pastix_complex64_t *b )
{
#if defined(PASTIX_WITH_MPI)
    pastix_int_t colnbr = cblk_colnbr(cblk);
    int rc;

    assert( colnbr <= solvmtx->colmax );
    assert( cblk->cblktype & CBLK_RECV );

#if defined (PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Bwd: Send cblk %ld to %2d at index %ld of size %ld\n",
             solvmtx->clustnum, (long)cblk->gcblknum, cblk->ownerid,
             (long)cblk->lcolidx, (long)colnbr );
#endif

    rc = MPI_Send( b + cblk->lcolidx, colnbr, PASTIX_MPI_COMPLEX64,
                   cblk->ownerid, cblk->gcblknum, solvmtx->solv_comm );

    assert( rc == MPI_SUCCESS );
    (void)rc;
#else
    (void)solvmtx;
    (void)cblk;
    (void)b;
#endif
}

/**
 *******************************************************************************
 *
 * @brief Receive the rhs associated to a cblk->lcolidx to the remote node.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix holding the communicator.
 *
 * @param[in] cblk
 *          The cblk which may define the part to sent.
 *
 * @param[inout] b
 *          The rhs that will be receive from the cblk->ownerid.
 *
 *******************************************************************************/
void
cpucblk_zrecv_rhs_backward( const SolverMatrix *solvmtx,
                            SolverCblk         *cblk,
                            pastix_complex64_t *b )
{
#if defined(PASTIX_WITH_MPI)
    MPI_Status   status;
    pastix_int_t colnbr = cblk_colnbr(cblk);
    int rc;

    assert( colnbr <= solvmtx->colmax );
    assert( cblk->cblktype & CBLK_FANIN );

#if defined (PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Bwd: Recv cblk %ld from %ld at index %ld of size %ld\n",
             solvmtx->clustnum, (long)cblk->gcblknum, (long)cblk->ownerid,
             (long)cblk->lcolidx, (long)colnbr );
#endif

    rc = MPI_Recv( b + cblk->lcolidx, colnbr, PASTIX_MPI_COMPLEX64,
                   cblk->ownerid, cblk->gcblknum, solvmtx->solv_comm, &status );
    assert( rc == MPI_SUCCESS );

#if defined (PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Bwd: Received cblk %ld from %2d\n",
             solvmtx->clustnum, (long)cblk->gcblknum, status.MPI_SOURCE );
#endif

    (void)rc;
#else
    (void)solvmtx;
    (void)cblk;
    (void)b;
#endif
}

/**
 *******************************************************************************
 *
 * @brief Receive the rhs associated to a cblk->lcolidx to the remote node.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix holding the communicator.
 *
 * @param[in] cblk
 *          The cblk which may define the part to sent.
 *
 * @param[inout] work
 *          The temporary buffer to receive the remote data
 *
 * @param[inout] b
 *          The rhs that will be updated by the reception.
 *
 * @param[in] ldb
 *          The leading dimension of the matrix b.
 *
 *******************************************************************************/
void
cpucblk_zrecv_rhs_forward( const SolverMatrix *solvmtx,
                           SolverCblk         *cblk,
                           pastix_complex64_t *work,
                           pastix_int_t        nrhs,
                           pastix_complex64_t *b,
                           pastix_int_t        ldb )
{
#if defined(PASTIX_WITH_MPI)
    MPI_Status   status;
    pastix_int_t colnbr = cblk_colnbr(cblk);
    int rc;

    assert( colnbr <= solvmtx->colmax );
    assert( cblk->cblktype & CBLK_RECV );

#if defined (PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Fwd: Recv cblk %ld from %ld at index %ld of size %ld\n",
             solvmtx->clustnum, (long)cblk->gcblknum, (long)cblk->ownerid,
             (long)cblk->lcolidx, (long)colnbr );
#endif

    rc = MPI_Recv( work, colnbr, PASTIX_MPI_COMPLEX64,
                   MPI_ANY_SOURCE, cblk->gcblknum, solvmtx->solv_comm, &status );
    assert( rc == MPI_SUCCESS );

#if defined (PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Fwd: Received cblk %ld from %2d\n",
                     solvmtx->clustnum, (long)cblk->gcblknum, status.MPI_SOURCE );
#endif

    core_zgeadd( PastixNoTrans, colnbr, nrhs,
                 1., work, ldb,
                 1., b + cblk->lcolidx, ldb );

    (void)rc;
#else
    (void)solvmtx;
    (void)cblk;
    (void)work;
    (void)nrhs;
    (void)b;
    (void)ldb;
#endif
}
