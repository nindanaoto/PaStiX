/**
 *
 * @file solver_matrix_gen.c
 *
 * PaStiX solver structure generation function.
 *
 * @copyright 1998-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Tony Delarue
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @date 2020-01-26
 *
 **/
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>

#include "common.h"
#include "symbol.h"
#include "solver.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"
#include "pastix/order.h"
#include "extendVector.h"
#include "simu.h"
#include "blendctrl.h"
#include "solver_matrix_gen_utils.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_blend
 *
 * @brief Initialize the solver matrix structure
 *
 * This function takes all the global preprocessing steps: the symbol matrix,
 * and the resul of the simulation step to generate the solver matrix that hold
 * only local information to each PaStiX process.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          On entry, the allocated pointer to a solver matrix structure.
 *          On exit, this structure holds alls the local information required to
 *          perform the numerical factorization.
 *
 * @param[in] symbmtx
 *          The global symbol matrix structure.
 *
 * @param[in] ordeptr
 *          The ordering structure.
 *
 * @param[in] simuctrl
 *          The information resulting from the simulation that will provide the
 *          data mapping, and the order of the task execution for the static
 *          scheduling.
 *
 * @param[in] ctrl
 *          The blend control structure that contains extra information
 *          computed during the analyze steps and the parameters of the analyze
 *          step.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if success.
 * @retval PASTIX_ERR_OUTOFMEMORY if one of the malloc failed.
 *
 *******************************************************************************/
int
solverMatrixGen( SolverMatrix          *solvmtx,
                 const symbol_matrix_t *symbmtx,
                 const pastix_order_t  *ordeptr,
                 const SimuCtrl        *simuctrl,
                 const BlendCtrl       *ctrl,
                 PASTIX_Comm            comm,
                 isched_t              *isched )
{
    pastix_int_t  i, j, k;
    pastix_int_t  coefnbr      = 0;
    pastix_int_t  nodenbr      = 0;
    pastix_int_t  cblknum      = 0;
    pastix_int_t  brownum      = 0;
    pastix_int_t *cblklocalnum = NULL;
    pastix_int_t *fcbklocalnum = NULL;
    pastix_int_t *pcbklocalnum = NULL;
    int          *recv_sources = NULL;
    pastix_int_t *bloklocalnum = NULL;
    pastix_int_t *tasklocalnum = NULL;
    pastix_int_t *browtmp      = NULL;
    pastix_int_t  flaglocal    = 0;
    (void)ordeptr;

    assert( symbmtx->dof == 1 );
    assert( symbmtx->baseval == 0 );

    solverInit( solvmtx );

#ifdef PASTIX_DYNSCHED
    solvmtx->btree = ctrl->btree;
#endif
    solvmtx->clustnum  = ctrl->clustnum;
    solvmtx->clustnbr  = ctrl->clustnbr;
    solvmtx->procnbr   = ctrl->total_nbcores;
    solvmtx->thrdnbr   = ctrl->local_nbthrds;
    solvmtx->bublnbr   = ctrl->local_nbctxts;
    solvmtx->solv_comm = comm;
    solvmtx->faninnbr  = 0;
    solvmtx->fanincnt  = 0;
    solvmtx->recvnbr   = 0;
    solvmtx->recvcnt   = 0;
    solvmtx->reqnum    = 0;

    /* Allocate the different local numbering arrays */
    MALLOC_INTERN( bloklocalnum, symbmtx->bloknbr,  pastix_int_t );
    MALLOC_INTERN( cblklocalnum, symbmtx->cblknbr,  pastix_int_t );
    MALLOC_INTERN( fcbklocalnum, symbmtx->cblknbr,  pastix_int_t );
    MALLOC_INTERN( pcbklocalnum, symbmtx->cblknbr,  pastix_int_t );
    MALLOC_INTERN( tasklocalnum, simuctrl->tasknbr, pastix_int_t );
    MALLOC_INTERN( recv_sources, symbmtx->cblknbr * solvmtx->clustnbr, int );

    /* Compute local indexes to compress the symbol information into solver */
    solvMatGen_fill_localnums( symbmtx, simuctrl, solvmtx,
                               cblklocalnum, bloklocalnum, tasklocalnum,
                               fcbklocalnum, pcbklocalnum, &recv_sources );

    /***************************************************************************
     * Fill in the local bloktab and cblktab
     */
    /* Allocate the cblktab and bloktab with the computed size */
    MALLOC_INTERN( solvmtx->cblktab,  solvmtx->cblknbr + 1, SolverCblk   );
    MALLOC_INTERN( solvmtx->bloktab,  solvmtx->bloknbr + 1, SolverBlok   );
    MALLOC_INTERN( solvmtx->browtab,  solvmtx->brownbr,     pastix_int_t );
    MALLOC_INTERN( browtmp,           symbmtx->browmax,     pastix_int_t );
    MALLOC_INTERN( solvmtx->gcbl2loc, symbmtx->cblknbr,     pastix_int_t );
    memset( solvmtx->gcbl2loc, 0xff,  symbmtx->cblknbr * sizeof(pastix_int_t) );
    {
        SolverCblk    *solvcblk = solvmtx->cblktab;
        SolverBlok    *solvblok = solvmtx->bloktab;
        symbol_cblk_t *symbcblk = symbmtx->cblktab;
        symbol_blok_t *symbblok = symbmtx->bloktab;
        SimuBlok      *simublok = simuctrl->bloktab;
        Cand          *candcblk = ctrl->candtab;
        pastix_int_t   blokamax = 0; /* Maximum area of a block in the global matrix */
        pastix_int_t   nbcblk2d = 0;
        pastix_int_t   nbblok2d = 0;
        pastix_int_t   gbloknm  = 0;
        pastix_int_t   bcscidx  = 0; /* Index of the local classic cblk */
        pastix_int_t   sndeidx  = 0; /* Index of the current supernode in the original elimination tree */
        pastix_int_t   recvcnt  = 0;

        solvmtx->cblkmax1d  = -1;
        solvmtx->cblkmin2d  = solvmtx->cblknbr;
        solvmtx->cblkmaxblk = 1;
        solvmtx->cblkschur  = solvmtx->cblknbr;
        solvmtx->gcblknbr   = symbmtx->cblknbr;
        solvmtx->maxrecv    = 0;
        solvmtx->colmax     = 0;
        cblknum             = 0;
        brownum             = 0;
        nodenbr             = 0;
        coefnbr             = 0;
        for ( i = 0; i < symbmtx->cblknbr; i++, symbcblk++, candcblk++ ) {
            SolverBlok  *fblokptr;
            pastix_int_t fbloknum = symbcblk[0].bloknum;
            pastix_int_t lbloknum = symbcblk[1].bloknum;
            pastix_int_t stride   = 0;
            pastix_int_t fcolnum, lcolnum, nbcols;
            pastix_int_t tasks2D  = candcblk->cblktype & CBLK_TASKS_2D;

            nbcols = solvMatGen_get_colnum( symbmtx, symbcblk, &fcolnum, &lcolnum );

            recvcnt = pcbklocalnum[i];
            for( k = 0; k < recvcnt; k++ )
            {
                fblokptr = solvblok;
                solvMatGen_init_cblk_recv( symbmtx, solvcblk, solvblok, candcblk,
                                           cblklocalnum, pcbklocalnum[i],
                                           fcolnum, lcolnum, brownum, nodenbr, i,
                                           recv_sources[solvmtx->recvcnt] );

                pastix_int_t cblksize = cblk_colnbr(solvcblk) * solvcblk->stride;
                solvmtx->maxrecv = pastix_imax( solvmtx->maxrecv, cblksize );
                pcbklocalnum[i]--;
                solvmtx->recvcnt++;
                cblknum++;
                solvcblk++;
                solvblok += (symbcblk[1].bloknum - symbcblk[0].bloknum);
            }

            flaglocal = 0;
            fblokptr  = solvblok;
            for ( j = fbloknum; j < lbloknum; j++, symbblok++, simublok++ ) {
                pastix_int_t frownum, lrownum, nbrows;

                nbrows = solvMatGen_get_rownum( symbmtx, symbblok, &frownum, &lrownum );

                /* Save max block size */
                blokamax = pastix_imax( blokamax, nbrows * nbcols );

                if ( ( simublok->ownerclust == ctrl->clustnum ) || ( fcbklocalnum[i] != -1 ) ) {
                    pastix_int_t fcblknm = ( simublok->ownerclust == ctrl->clustnum ) ? cblklocalnum[symbblok->fcblknm] : -1;
                    flaglocal = 1;

                    /* Init the blok */
                    solvMatGen_init_blok( solvblok,
                                          cblklocalnum[symbblok->lcblknm], fcblknm,
                                          frownum, lrownum, stride, nbcols,
                                          candcblk->cblktype & CBLK_LAYOUT_2D );
                    solvblok->gbloknm = gbloknm;
                    stride += nbrows;
                    solvblok++;
                }
                gbloknm++;
            }

            if ( flaglocal ) {
                pastix_int_t brownbr;
                solvmtx->colmax = pastix_imax( solvmtx->colmax, nbcols );

                /* Update the 1D/2D infos of the solvmtx through the CBLK_RECV. */
                for( k = 1; k <= recvcnt; k++ )
                {
                    solvMatGen_cblkIs2D( solvmtx, &nbcblk2d, &nbblok2d,
                                        (solvblok - fblokptr), tasks2D, cblknum - k );
                }

                /* Update the 1D/2D infos of the solvmtx through a cblk. */
                solvMatGen_cblkIs2D( solvmtx, &nbcblk2d, &nbblok2d,
                                    (solvblok - fblokptr), tasks2D, cblknum );

                /* Init the cblk */
                solvMatGen_init_cblk( solvcblk, fblokptr, candcblk, symbcblk,
                                      fcolnum, lcolnum, brownum, stride, nodenbr, i, ctrl->clustnum);

#if defined(PASTIX_WITH_MPI)
                if ( solvmtx->clustnbr > 1 )
                {
                    /* No low-rank compression in distributed for the moment */
                    if( solvcblk->cblktype & CBLK_COMPRESSED ) {
                        static int warning_compressed = 1;
                        if ( warning_compressed && (solvmtx->clustnum == 0) ) {
                            fprintf( stderr,
                                     "Warning: Low-rank compression is not yet available with multiple MPI processes\n"
                                     "         It is thus disabled and the factorization will be performed in full-rank\n" );
                            warning_compressed = 0;
                        }
                        solvcblk->cblktype &= (~CBLK_COMPRESSED);
                    }

                    /* No Schur complement in distributed for the moment */
                    if( solvcblk->cblktype & CBLK_IN_SCHUR ) {
                        static int warning_schur = 1;
                        if ( warning_schur && (solvmtx->clustnum == 0) ) {
                            fprintf( stderr,
                                     "Warning: Schur complement support is not yet available with multiple MPI processes\n"
                                     "         It is thus disabled and the factorization will be fully performed\n" );
                            warning_schur = 0;
                        }
                        solvcblk->cblktype &= (~CBLK_IN_SCHUR);
                    }
                }
#endif

                /* Store first local cblk in Schur */
                if ( (cblknum < solvmtx->cblkschur) &&
                     (solvcblk->cblktype & CBLK_IN_SCHUR) )
                {
                    solvmtx->cblkschur = cblknum;
                }

                solvmtx->gcbl2loc[i] = solvcblk - solvmtx->cblktab;

                /* If the cblk is a fanin one */
                if ( fcbklocalnum[i] != -1 ) {
                    solvcblk->cblktype |= CBLK_FANIN;
                    solvcblk->ownerid   = fcbklocalnum[i];
                    solvmtx->fanincnt++;
                }
                else {
                    solvcblk->bcscnum = bcscidx;
                    bcscidx++;
                }

                /* Compute the original supernode in the nested dissection */
                sndeidx = solvMatGen_supernode_index( solvcblk, sndeidx, ordeptr );

                /*
                 * Copy browtab information
                 * In case of 2D tasks, we reorder the browtab to first store
                 * the 1D contributions, and then the 2D updates.
                 * This might also be used for low rank compression, to first
                 * accumulate small dense contributions, and then, switch to a
                 * low rank - low rank update scheme.
                 */
                brownbr = solvMatGen_reorder_browtab( symbmtx, symbcblk, solvmtx, solvcblk,
                                                      browtmp, cblklocalnum, bloklocalnum, brownum );

                /* Diagonal bloks of CBLK_RECV are added here in the browtab */
                while( recvcnt ) {
                    fblokptr = solvcblk[-recvcnt].fblokptr;
                    solvmtx->browtab[brownum + brownbr] = fblokptr - solvmtx->bloktab;
                    fblokptr->browind = brownum + brownbr;
                    brownbr++;

                    /* Supernode index is copied in too */
                    solvcblk[-recvcnt].sndeidx = solvcblk->sndeidx;
                    recvcnt--;
                }

                brownum += brownbr;

                assert( brownum <= solvmtx->brownbr );
                assert( solvcblk->brown2d >= solvcblk->brownum );
                assert( solvcblk->brown2d <= solvcblk->brownum + brownbr );

                cblknum++;
                solvcblk++;
            }

            /* Extra statistic informations, need to be done everywhere */
            nodenbr += nbcols;
            coefnbr += stride * nbcols;
        }

        /*  Add a virtual cblk to avoid side effect in the loops on cblk bloks */
        if ( cblknum > 0 ) {
            solvMatGen_init_cblk( solvcblk, solvblok, candcblk, symbcblk,
                                  solvcblk[-1].lcolnum + 1, solvcblk[-1].lcolnum + 1,
                                  brownum, 0, nodenbr, -1, solvmtx->clustnum );
        }

        /*  Add a virtual blok to avoid side effect in the loops on cblk bloks */
        if ( solvmtx->bloknbr > 0 ) {
            solvMatGen_init_blok( solvblok, symbmtx->cblknbr + 1, symbmtx->cblknbr + 1,
                                  solvcblk[-1].lcolnum + 1, solvcblk[-1].lcolnum + 1,
                                  0, 0, 0 );
        }

        solvmtx->nodenbr = nodenbr;
        solvmtx->coefnbr = coefnbr;
        solvmtx->arftmax = blokamax;

        solvmtx->nb2dcblk = nbcblk2d;
        solvmtx->nb2dblok = nbblok2d;

        assert( solvmtx->cblkmax1d + 1 >= solvmtx->cblkmin2d );
        assert( solvmtx->brownbr  == brownum );
        assert( solvmtx->cblknbr  == cblknum );
        assert( solvmtx->faninnbr == solvmtx->fanincnt );
        assert( solvmtx->recvnbr  == solvmtx->recvcnt );
        assert( solvmtx->bloknbr  == solvblok - solvmtx->bloktab );
    }
    memFree_null( browtmp );

    /* Fill in tasktab */
    solvMatGen_fill_tasktab( solvmtx, isched, simuctrl,
                             tasklocalnum, cblklocalnum,
                             bloklocalnum, ctrl->clustnum, 0 );

    memFree_null(cblklocalnum);
    memFree_null(fcbklocalnum);
    memFree_null(pcbklocalnum);
    memFree_null(bloklocalnum);
    memFree_null(tasklocalnum);
    memFree_null(recv_sources);

    /* Compute the maximum area of the temporary buffer */
    solvMatGen_max_buffers( solvmtx );
    solvMatGen_stats_last( solvmtx );

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_blend
 *
 * @brief Initialize the solver matrix structure in sequential
 *
 * This function takes all the global preprocessing steps: the symbol matrix,
 * and the result of the simulation step to generate the solver matrix for one
 * PaStiX process.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          On entry, the allocated pointer to a solver matrix structure.
 *          On exit, this structure holds alls the local information required to
 *          perform the numerical factorization.
 *
 * @param[in] symbmtx
 *          The global symbol matrix structure.
 *
 * @param[in] ordeptr
 *          The ordering structure.
 *
 * @param[in] simuctrl
 *          The information resulting from the simulation that will provide the
 *          data mapping, and the order of the task execution for the static
 *          scheduling.
 *
 * @param[in] ctrl
 *          The blend control structure that contains extra information
 *          computed during the analyze steps and the parameters of the analyze
 *          step.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if success.
 * @retval PASTIX_ERR_OUTOFMEMORY if one of the malloc failed.
 *
 *******************************************************************************/
int
solverMatrixGenSeq( SolverMatrix          *solvmtx,
                    const symbol_matrix_t *symbmtx,
                    const pastix_order_t  *ordeptr,
                    const SimuCtrl        *simuctrl,
                    const BlendCtrl       *ctrl,
                    PASTIX_Comm            comm,
                    isched_t              *isched,
                    pastix_int_t           is_dbg )
{
    pastix_int_t  i, j;
    pastix_int_t  coefnbr = 0;
    pastix_int_t  nodenbr = 0;
    pastix_int_t  cblknum = 0;
    pastix_int_t  brownum = 0;
    pastix_int_t *browtmp = 0;
    (void)ordeptr;

    assert( symbmtx->dof == 1 );
    assert( symbmtx->baseval == 0 );

    solverInit( solvmtx );

    solvmtx->clustnum  = ctrl->clustnum;
    solvmtx->clustnbr  = ctrl->clustnbr;
    solvmtx->procnbr   = ctrl->total_nbcores;
    solvmtx->thrdnbr   = ctrl->local_nbthrds;
    solvmtx->bublnbr   = ctrl->local_nbctxts;
    solvmtx->solv_comm = comm;
    solvmtx->faninnbr  = 0;
    solvmtx->fanincnt  = 0;
    solvmtx->recvnbr   = 0;
    solvmtx->recvcnt   = 0;
    solvmtx->reqnum    = 0;

    solvmtx->tasknbr = simuctrl->tasknbr;
    solvmtx->cblknbr = symbmtx->cblknbr;
    solvmtx->bloknbr = symbmtx->bloknbr;
    solvmtx->brownbr = symbmtx->cblktab[ solvmtx->cblknbr ].brownum
                     - symbmtx->cblktab[0].brownum;

    /***************************************************************************
     * Fill in the local bloktab and cblktab
     */
    /* Allocate the cblktab and bloktab with the computed size */
    MALLOC_INTERN(solvmtx->cblktab, solvmtx->cblknbr+1, SolverCblk  );
    MALLOC_INTERN(solvmtx->bloktab, solvmtx->bloknbr+1, SolverBlok  );
    MALLOC_INTERN(solvmtx->browtab, solvmtx->brownbr,   pastix_int_t);
    MALLOC_INTERN(browtmp,          symbmtx->browmax,   pastix_int_t);
    {
        SolverCblk    *solvcblk = solvmtx->cblktab;
        SolverBlok    *solvblok = solvmtx->bloktab;
        symbol_cblk_t *symbcblk = symbmtx->cblktab;
        symbol_blok_t *symbblok = symbmtx->bloktab;
        SimuBlok      *simublok = simuctrl->bloktab;
        Cand          *candcblk = ctrl->candtab;
        pastix_int_t   blokamax = 0; /* Maximum area of a block in the global matrix */
        pastix_int_t   nbcblk2d = 0;
        pastix_int_t   nbblok2d = 0;
        pastix_int_t   sndeidx  = 0; /* Index of the current supernode in the original elimination tree */

        solvmtx->cblkmax1d  = -1;
        solvmtx->cblkmin2d  = solvmtx->cblknbr;
        solvmtx->cblkmaxblk = 1;
        solvmtx->cblkschur  = solvmtx->cblknbr;
        solvmtx->gcblknbr   = symbmtx->cblknbr;
        solvmtx->maxrecv    = 0;
        solvmtx->colmax     = 0;
        cblknum = 0;
        brownum = 0;
        nodenbr = 0;
        coefnbr = 0;
        for(i=0; i<symbmtx->cblknbr; i++, symbcblk++, candcblk++)
        {
            SolverBlok  *fblokptr = solvblok;
            pastix_int_t fbloknum = symbcblk[0].bloknum;
            pastix_int_t lbloknum = symbcblk[1].bloknum;
            pastix_int_t stride   = 0;
            pastix_int_t fcolnum, lcolnum, nbcols;
            pastix_int_t tasks2D  = candcblk->cblktype & CBLK_TASKS_2D;
            pastix_int_t brownbr;

            nbcols = solvMatGen_get_colnum( symbmtx, symbcblk, &fcolnum, &lcolnum );

            for( j=fbloknum; j<lbloknum; j++, symbblok++, simublok++ ) {
                pastix_int_t frownum, lrownum, nbrows;

                nbrows = solvMatGen_get_rownum( symbmtx, symbblok, &frownum, &lrownum );

                /* Save the max block size */
                blokamax = pastix_imax( blokamax, nbrows * nbcols );

                /* Init the blok */
                solvMatGen_init_blok( solvblok,
                                      symbblok->lcblknm, symbblok->fcblknm,
                                      frownum, lrownum, stride, nbcols,
                                      candcblk->cblktype & CBLK_LAYOUT_2D );
                stride += nbrows;
                solvblok++;
            }

            solvMatGen_cblkIs2D( solvmtx, &nbcblk2d, &nbblok2d,
                                 (solvblok - fblokptr), tasks2D, cblknum );

            /* Init the cblk */
            solvMatGen_init_cblk( solvcblk, fblokptr, candcblk, symbcblk,
                                  fcolnum, lcolnum, brownum, stride, nodenbr, i, ctrl->clustnum );
            solvcblk->bcscnum = i;
#if defined(PASTIX_WITH_MPI)
            /* No low-rank compression in distributed for the moment */
            if ( (solvmtx->clustnbr > 1) && (solvcblk->cblktype & CBLK_COMPRESSED) )
            {
                static int warning_compressed = 1;
                if ( warning_compressed && (solvmtx->clustnum == 0) ) {
                    fprintf( stderr,
                             "Warning: Low-rank compression is not available yet with multiple MPI processes\n"
                             "         It is thus disabled and the factorization will be performed in full-rank\n" );
                    warning_compressed = 0;
                }
                solvcblk->cblktype &= (~CBLK_COMPRESSED);
            }
#endif

            /* Store first local cblk in Schur */
            if ( (cblknum < solvmtx->cblkschur) &&
                 (solvcblk->cblktype & CBLK_IN_SCHUR) )
            {
                solvmtx->cblkschur = cblknum;
            }

            /* Every cblk is local */
            solvcblk->ownerid = simublok[-1].ownerclust;
            assert( nodenbr == fcolnum );

            /* Compute the original supernode in the nested dissection */
            sndeidx = solvMatGen_supernode_index( solvcblk, sndeidx, ordeptr );

            /*
             * Copy browtab information
             * In case of 2D tasks, we reorder the browtab to first store
             * the 1D contributions, and then the 2D updates.
             * This might also be used for low rank compression, to first
             * accumulate small dense contributions, and then, switch to a
             * low rank - low rank update scheme.
             */
            brownbr = solvMatGen_reorder_browtab( symbmtx, symbcblk, solvmtx, solvcblk,
                                                  browtmp, NULL, NULL, brownum );

            brownum += brownbr;

            assert( brownum <= solvmtx->brownbr );
            assert( solvcblk->brown2d >= solvcblk->brownum );
            assert( solvcblk->brown2d <= solvcblk->brownum + brownbr );

            /* Extra statistic informations */
            nodenbr += nbcols;
            coefnbr += stride * nbcols;

            cblknum++;
            solvcblk++;
        }

        /*  Add a virtual cblk to avoid side effect in the loops on cblk bloks */
        if ( cblknum > 0 ) {
            solvMatGen_init_cblk( solvcblk, solvblok, candcblk, symbcblk,
                                  solvcblk[-1].lcolnum + 1, solvcblk[-1].lcolnum + 1,
                                  symbcblk->brownum, 0, nodenbr, -1, ctrl->clustnum);
        }

        /*  Add a virtual blok to avoid side effect in the loops on cblk bloks */
        if ( solvmtx->bloknbr > 0 ) {
            solvMatGen_init_blok( solvblok, symbmtx->cblknbr + 1, symbmtx->cblknbr + 1,
                                  solvcblk[-1].lcolnum + 1, solvcblk[-1].lcolnum + 1,
                                  0, 0, 0 );
        }

        solvmtx->nodenbr = nodenbr;
        solvmtx->coefnbr = coefnbr;
        solvmtx->arftmax = blokamax;

        solvmtx->nb2dcblk = nbcblk2d;
        solvmtx->nb2dblok = nbblok2d;

        assert( solvmtx->cblkmax1d+1 >= solvmtx->cblkmin2d );
        assert( solvmtx->cblknbr == cblknum );
        assert( solvmtx->bloknbr == solvblok - solvmtx->bloktab );
    }
    memFree_null( browtmp );

    /*
     * Update browind fields
     */
    for(i=0; i<solvmtx->brownbr; i++)
    {
        pastix_int_t bloknum = solvmtx->browtab[i];
        solvmtx->bloktab[ bloknum ].browind = i;
    }

    /* Fill in tasktab */
    solvMatGen_fill_tasktab( solvmtx, isched, simuctrl,
                             NULL, NULL, NULL, ctrl->clustnum, is_dbg );

    /* Compute the maximum area of the temporary buffer */
    solvMatGen_max_buffers( solvmtx );
    solvMatGen_stats_last( solvmtx );

    return PASTIX_SUCCESS;
}
