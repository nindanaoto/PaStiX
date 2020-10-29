/**
 *
 * @file solver_check.c
 *
 * PaStiX check function fo rthe solver structure.
 *
 * @copyright 2004-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.1.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2020-01-26
 *
 **/
#include <stdio.h>
#include <assert.h>

#include "common.h"
#include "symbol.h"
#include "queue.h"
#include "solver.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"
#include "extendVector.h"
#include "blendctrl.h"
#include "simu.h"
#include "pastix_zcores.h"
#include "pastix_ccores.h"
#include "pastix_dcores.h"
#include "pastix_scores.h"

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Checks the consistency of the given solver matrix structure.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix structure to check.
 *
 *******************************************************************************
 *
 * @retval 0 if the structure is correct
 * @retval 1 if incorrect
 *
 *******************************************************************************/
int
solverCheck( const SolverMatrix *solvmtx )
{
    int          i, j, browtype, cblktype;
    SolverBlok  *blok, *fblok;
    SolverCblk  *cblk, *fcblk;
    pastix_int_t bloknum;

    assert( (solvmtx->cblknbr - solvmtx->recvnbr - solvmtx->faninnbr ) <= solvmtx->gcblknbr );

    cblk = solvmtx->cblktab;

    for ( i = 0; i < solvmtx->cblknbr; i++, cblk++ ) {
        /* Make sure the lock is initialized correctly */
        assert( cblk->lock == PASTIX_ATOMIC_UNLOCKED );

        /* Check dimensions */
        assert( cblk->fcolnum <= cblk->lcolnum );
        assert( cblk->stride >= cblk_colnbr( cblk ) );

        /* Check the range of value for the ownerid */
        assert( (cblk->ownerid  >=  0) && (cblk->ownerid  < solvmtx->clustnbr) );
        assert( (cblk->threadid >= -1) && (cblk->threadid < solvmtx->thrdnbr ) );

        /* Check that pointers have been initialized to null */
        assert( cblk->lcoeftab == NULL );
        assert( cblk->ucoeftab == NULL );

        /* Check that we have at least one block */
        assert( (cblk[1].fblokptr - cblk[0].fblokptr) >= 1 );

        fblok = cblk->fblokptr;
        fcblk = solvmtx->cblktab + fblok->fcblknm;

        if ( cblk->cblktype & CBLK_RECV ) {
            /*
             * Check the redundancy of information, the facing cblk of the first
             * block should match the index of the local cblk which refers to
             * the local non recv version of this cblk.
             */
            assert( fcblk->lcolidx  == cblk->lcolidx  );
            assert( fcblk->sndeidx  == cblk->sndeidx  );
            assert( fcblk->gcblknum == cblk->gcblknum );
            assert( fcblk->fblokptr->lcblknm == fblok->fcblknm );
            assert( i < fcblk->fblokptr->lcblknm );

            /* It has to be remote */
            assert( cblk->ownerid != solvmtx->clustnum );

            /* No index in the bcsctab */
            assert( -1 == cblk->bcscnum );

            /* It has to be included in the destination */
            assert( cblk->fcolnum >= fcblk->fcolnum );
            assert( cblk->lcolnum <= fcblk->lcolnum );

            /* It has no input contributions */
            assert( ( cblk[0].brownum == cblk[0].brown2d ) &&
                    ( cblk[0].brown2d == cblk[1].brownum ) );

            /* Check that reception blocks are included within local reception blocks */
            {
                SolverBlok *rblok;
                blok  = fblok;
                rblok = fcblk->fblokptr;
                for ( ; blok < cblk[1].fblokptr; blok++ ) {
                    while( !is_block_inside_fblock( blok, rblok ) ) {
                        rblok++;
                        assert( rblok->lcblknm == (fcblk - solvmtx->cblktab) );
                    }
                }
                assert( rblok->lcblknm == (fcblk - solvmtx->cblktab) );
            }
        }
        else {
            /* It has to be included in the destination */
            assert( cblk[0].lcolnum < cblk[1].fcolnum );

            /* Check that we have none or some contributions */
            assert( cblk[1].brownum >= cblk[0].brownum );
            assert( ( cblk[0].brownum <= cblk[0].brown2d ) &&
                    ( cblk[0].brown2d <= cblk[1].brownum ) );

            /* First diagonal block should not appear in the browtab */
            assert( fblok->browind == -1 );

            if ( cblk->cblktype & CBLK_FANIN ) {
                assert( -1 == cblk->bcscnum );

                /* It has to be remote */
                assert( cblk->ownerid != solvmtx->clustnum );

                /* Fanin targets do not contribute locally so we don't know the target */
                assert( fblok->fcblknm == -1 );

                /* Uncomment when rhs will be distributed */
                // assert( cblk->lcolidx == -1 );
            }
            else {
                /* Check that first diagonal block belongs to ourself */
                assert( fblok->fcblknm == fblok->lcblknm );

                assert( cblk->bcscnum != -1 );

                /* It has to be local */
                assert( (solvmtx->gcbl2loc == NULL) ||
                        (cblk->ownerid == solvmtx->clustnum) );

                /* Check if possible that the right hand side index are in increasing order */
                assert( (cblk[1].cblktype & CBLK_FANIN) ||
                        (!(cblk[1].cblktype & CBLK_FANIN) && (cblk[0].lcolidx < cblk[1].lcolidx)) );
            }
        }

        /* Check bloks in the current cblk */
        blok = fblok;
        for ( ; blok < cblk[1].fblokptr; blok++ ) {
            assert( blok->lcblknm == i );
            assert( blok->frownum <= blok->lrownum );

            /* Next block in the same cblk must be after */
            assert( (blok[1].lcblknm != blok[0].lcblknm) ||
                    ((blok[1].lcblknm == blok[0].lcblknm) && (blok[0].lrownum < blok[1].frownum)) );

            if ( cblk->cblktype & CBLK_FANIN ) {
                assert( blok->fcblknm == -1 );
            }
            else if ( cblk->cblktype & CBLK_RECV ) {
                if ( blok == cblk->fblokptr ) {
                    assert( blok->fcblknm != -1 );
                    assert( blok->browind != -1 );
                    assert( blok->lcblknm < blok->fcblknm );
                }
                else {
                    assert( blok->fcblknm == -1 );
                    assert( blok->browind == -1 );
                }
            }
            else {
                fcblk = solvmtx->cblktab + blok->fcblknm;
                fblok = fcblk->fblokptr;

                assert( ((blok == cblk->fblokptr) && (blok->lcblknm == blok->fcblknm)) ||
                        (blok->lcblknm < blok->fcblknm) );
                assert( blok->frownum >= fblok->frownum );
                assert( blok->lrownum <= fblok->lrownum );
            }
        }

        /* Check previous bloks in row */
        browtype = 0;
        cblktype = 0;
        for ( j = cblk[0].brownum; j < cblk[1].brownum; j++ ) {
            bloknum = solvmtx->browtab[j];
            blok    = solvmtx->bloktab + bloknum;
            fcblk   = solvmtx->cblktab + blok->lcblknm;

            assert( blok->browind == j );
            assert( blok->fcblknm == i );

            /* Blocks in Fanin sould never appear in the browtab */
            assert( !(fcblk->cblktype & CBLK_FANIN) );

            /* The contribution comes from a RECV cblk */
            if ( fcblk->cblktype & CBLK_RECV ) {
                assert( fcblk->ownerid != solvmtx->clustnum );

                /* It can come only from the diagonal block */
                assert( blok == fcblk->fblokptr );
            }
            else {
                assert( (solvmtx->gcbl2loc == NULL) || (fcblk->ownerid == solvmtx->clustnum) );
            }

            /* Check that the brow is sorted correctly (1D, 2D, RECV) */
            cblktype = fcblk->cblktype & ( CBLK_TASKS_2D | CBLK_RECV );
            assert( cblktype >= browtype );
            browtype = browtype | cblktype; /* Take the max for the next step */
        }
    }

    return 0;
}

