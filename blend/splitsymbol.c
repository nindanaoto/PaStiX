/**
 *
 * @file splitsymbol.c
 *
 * PaStiX simulation task basic functions.
 *
 * @copyright 2004-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2020-02-05
 *
 * @addtogroup blend_dev_split
 * @{
 *   This module handles the splitting of the existing symbol matrix to generate
 *   more parallelism with finer granularity tasks.
 *
 **/
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "common.h"
#include "cost.h"
#include "symbol.h"
#include "elimintree.h"
#include "extendVector.h"
#include "cand.h"
#include "queue.h"
#include "blendctrl.h"
#include "solver.h"
#include "simu.h"

#include "extracblk.h"

static inline
int pastix_blend_with_constant_split() {
    return pastix_env_is_set_to("PASTIX_BLEND_SPLIT", "CONSTANT");
}

static inline
int pastix_blend_with_smallest_upper_split() {
    return pastix_env_is_set_to("PASTIX_BLEND_SPLIT", "UPPER");
}

static inline
int pastix_blend_split_percent() {
    return
      pastix_getenv_get_value_int( "PASTIX_BLEND_SPLIT_AUTORIZED_PERCENTAGE",
                                    10);
}

/**
 *******************************************************************************
 *
 * @brief Compute the number of cut for a given cblk width and number of candidates.
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The blend control structure with the parameters.
 *
 * @param[in] candnbr
 *          The number of candidates to compute the cblk.
 *
 * @param[in] width
 *          The width of a cblk to cut.
 *
 *******************************************************************************
 *
 * @return The number of cut that should be done for this width.
 *
 *******************************************************************************/
static inline pastix_int_t
computeNbSplit( const BlendCtrl *ctrl,
                pastix_int_t     width )
{
    pastix_int_t blas_min_col;
    pastix_int_t blas_max_col;
    pastix_int_t nseq;

    /* Compute minimun broadness for splitting this cblk */
    blas_min_col = ctrl->blcolmin;
    blas_max_col = ctrl->blcolmax;

    if( width <= blas_max_col) {
        return 1;
    }
    nseq = pastix_iceil( width, blas_max_col );

    /* Make sure cblk are at least blas_min_col wide */
    if ( nseq > 1 && (width / nseq) < blas_min_col ) {
        nseq--;
    }

    assert( nseq > 0 );
    return nseq;
}

/**
 *******************************************************************************
 *
 * @brief Generate an array with the number of blocks facing each row to
 * minimize the cuts.
 *
 *******************************************************************************
 *
 * @param[in] symbmtx
 *          The symbol matrix structure.
 *
 * @param[in] frowsplit
 *          The index of the first row that might be split to avoid computing it
 *          on the full matrix.
 *
 *******************************************************************************
 *
 * @return The array with the number of block on each row from frowsplit to n.
 *
 *******************************************************************************/
static inline pastix_int_t *
computeNbBlocksPerLine( const symbol_matrix_t *symbmtx,
                        pastix_int_t frowsplit )
{
    symbol_blok_t *curblok;
    pastix_int_t  *nblocksperline;
    pastix_int_t   bloknum, line;
    pastix_int_t   size = symbmtx->nodenbr - frowsplit + 1;

    /*
     * Allocate the temporary buffer nblocksperline, nbblocksperline stores the
     * number of blocks that will be splitted if with split between line i and
     * i+1.
     */
    MALLOC_INTERN( nblocksperline, size, pastix_int_t );
    memset( nblocksperline, 0, size * sizeof(pastix_int_t) );

    curblok = symbmtx->bloktab;
    for(bloknum=0; bloknum<symbmtx->bloknbr; bloknum++, curblok++ )
    {
        if ( curblok->lrownum < frowsplit ) {
            continue;
        }

        /*
         * For each couple of rows (i,i+1) in the block, we increment
         * the number blocks in regard of the row i
         */
        for(line = pastix_imax( curblok->frownum, frowsplit);
            line < curblok->lrownum; line++ )
        {
            nblocksperline[ line-frowsplit ]++;
        }
    }
    assert( nblocksperline[ size-1 ] == 0 );

    return nblocksperline;
}

typedef pastix_int_t (*split_fct_t)( const pastix_int_t *,
                                     pastix_int_t, pastix_int_t, pastix_int_t, pastix_int_t );

/**
 *******************************************************************************
 *
 * @brief Compute the cut that minimizes the generation of off-diagonal blocks.
 *
 * This function searched first for a smaller block than the initial step to
 * generate the cut, and then for a larger one that is return only if it is
 * strictly better than the first solution found.
 *
 *******************************************************************************
 *
 * @param[in] nblocksperline
 *          Array of size max with the number of blocks in front of the lines [0;max-1]
 *
 * @param[in] step
 *          The starting cut to refine. 0 < step
 *
 * @param[in] max
 *          The size of the nblocksperline array. nmax > 0.
 *
 * @param[in] authorized_percent
 *          The authorized percentage move around the initial cut.
 *
 *******************************************************************************
 *
 * @return The optimal cut found.
 *
 *******************************************************************************/
static inline pastix_int_t
computeSmallestSplit( const pastix_int_t *nblocksperline,
                      pastix_int_t step,
                      pastix_int_t min,
                      pastix_int_t max,
                      pastix_int_t authorized_percent )
{
    pastix_int_t limit = pastix_iceil( step*authorized_percent, 100 );
    pastix_int_t i, lcolnum, nbsplit;
    pastix_int_t lmin, lmax, lavg;

    if (step >= max) {
        return max-1;
    }
    assert( step > 1 );

    lavg = step - 1;
    lmin = pastix_imax( lavg - limit - 1,  min );
    lmax = pastix_imin( lavg + limit + 1,  max );

    lcolnum = lavg;
    nbsplit = nblocksperline[ lcolnum ];

    /* Search for the minimal split */
    for(i=lavg+1; i<lmax; i++ )
    {
        if ( nblocksperline[ i ] < nbsplit )
        {
            lcolnum = i;
            nbsplit = nblocksperline[ i ];
        }
    }
    for(i=lavg-1; i>lmin; i-- )
    {
        if ( nblocksperline[ i ] < nbsplit )
        {
            lcolnum = i;
            nbsplit = nblocksperline[ i ];
        }
    }

    return lcolnum;
}

/**
 *******************************************************************************
 *
 * @brief Compute the cut that minimizes the generation of off-diagonal blocks.
 *
 * This function searched for the largest cut in the range [step -
 * authorized_percent, step + authorized_percent] that minimizes the number of
 * off-diagonal blocks generated.
 *
 *******************************************************************************
 *
 * @param[in] nblocksperline
 *          Array of size max with the number of blocks in front of the lines [0;max-1]
 *
 * @param[in] step
 *          The starting cut to refine. 0 < step
 *
 * @param[in] max
 *          The size of the nblocksperline array. nmax > 0.
 *
 * @param[in] authorized_percent
 *          The authorized percentage move around the initial cut.
 *
 *******************************************************************************
 *
 * @return The optimal cut found.
 *
 *******************************************************************************/
static inline pastix_int_t
computeSmallestSplit_max( const pastix_int_t *nblocksperline,
                          pastix_int_t step,
                          pastix_int_t min,
                          pastix_int_t max,
                          pastix_int_t authorized_percent )
{
    pastix_int_t limit = pastix_iceil( step*authorized_percent, 100 );
    pastix_int_t i, lcolnum, nbsplit;
    pastix_int_t lmin, lmax, lavg;

    if ( step >= max ) {
        return max-1;
    }
    assert( step > 1 );

    lavg = step - 1;
    lmin = pastix_imax( lavg - limit,      min );
    lmax = pastix_imin( lavg + limit + 1,  max );

    lcolnum = lmin;
    nbsplit = nblocksperline[ lcolnum ];

    /* Search for the minimal split */
    for(i=lmin; i<lmax; i++ )
    {
        if ( nblocksperline[ i ] <= nbsplit )
        {
            lcolnum = i;
            nbsplit = nblocksperline[ i ];
        }
    }

    return lcolnum;
}

/**
 *******************************************************************************
 *
 * @brief Compute a constant cut with the given parameters.
 *
 *******************************************************************************
 *
 * @param[in] nblocksperline
 *          Array of size max with the number of blocks in front of the lines [0;max-1]
 *
 * @param[in] step
 *          The starting cut to refine. 0 < step
 *
 * @param[in] max
 *          The size of the nblocksperline array. nmax > 0.
 *
 * @param[in] authorized_percent
 *          The authorized percentage move around the initial cut.
 *
 *******************************************************************************
 *
 * @return The optimal cut found.
 *
 *******************************************************************************/
static inline pastix_int_t
computeConstantSplit( const pastix_int_t *nblocksperline __attribute__((unused)),
                      pastix_int_t step,
                      pastix_int_t min __attribute__((unused)),
                      pastix_int_t max,
                      pastix_int_t authorized_percent __attribute__((unused)) )
{
    if ( step >= max ) {
        return max-1;
    }
    assert( step > 1 );
    return step-1;
}

/**
 *******************************************************************************
 *
 * @brief Split the column blocks to generate parallelism
 *
 * This algorithm can use three different strategies to cut, two of them tries
 * to minimize the number of generated off-diagnonal blocks, one with the
 * minimal size, and the seoncd with teh maximal size in the error percentage
 * around the initially computed split. The third strategy cuts regularly the
 * column blocks without paying attention to the facing off-diagonal blocks.
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The blend control structure.
 *
 * @param[in] symbmtx
 *          The symbol matrix structure.
 *
 * @param[inout] extracblk
 *          The initialized structure to store the newly created column
 *          blocks. (Blocks are generated later during merge operation, see
 *          extraCblkMerge()).
 *
 * @param[in] candtab
 *          The candidate array for the initial symbo matrix structure.
 *
 *******************************************************************************/
void
splitSmart( const BlendCtrl       *ctrl,
            const symbol_matrix_t *symbmtx,
            ExtraCblk_t           *extracblk )
{
    symbol_cblk_t *curcblk;
    symbol_blok_t *curblok;
    pastix_int_t  *nblocksperline = NULL;
    pastix_int_t   cblknum, bloknum, line;
    pastix_int_t   fsplitrow = -1;
    pastix_int_t   authorized_percent;
    split_fct_t    split_fct = computeSmallestSplit;

    if (pastix_blend_with_constant_split()) {
        split_fct = computeConstantSplit;
    }
    else if (pastix_blend_with_smallest_upper_split()) {
        split_fct = computeSmallestSplit_max;
    }

    authorized_percent = pastix_blend_split_percent();

    curcblk = symbmtx->cblktab;
    for(cblknum = 0; cblknum<symbmtx->cblknbr; cblknum++, curcblk++)
    {
        pastix_int_t fcolnum = curcblk->fcolnum;
        pastix_int_t lcolnum = curcblk->lcolnum;
        pastix_int_t step, nseq;
        pastix_int_t width, lwidth;
        pastix_int_t fcol, lcol;
        pastix_int_t nbcblk = 0;

        if ( curcblk->selevtx == SYMBCBLK_KWAY ) {
            continue;
        }

        /*
         * Compute the number of cblk to be generated by split,
         * for instance we choose to split at the maximum
         */
        width = lcolnum - fcolnum + 1;
        nseq = computeNbSplit( ctrl, width );
        if ( nseq <= 1 ) {
            continue;
        }

        /* Initialize the nbblocksperline array */
        if ( fsplitrow == -1 ) {
            assert( fcolnum != -1 );
            fsplitrow = fcolnum;
            nblocksperline = computeNbBlocksPerLine( symbmtx, fsplitrow );
            nblocksperline -= fsplitrow;
        }

        /* Create the new cblk */
        fcol = fcolnum;
        while( fcol <= lcolnum )
        {
            /* Adapt the step to the remaining width */
            nseq = computeNbSplit( ctrl, width );
            if ( nseq == 1 ) {
                extraCblkAdd( extracblk, fcol, lcolnum,
                              curcblk->selevtx );
                nbcblk++;
                fcol = lcolnum + 1;
                continue;
            }

            step = pastix_iceil( width, nseq );
            assert( step >= ctrl->blcolmin );

            /* Make sure we leave enough columns */
            lwidth = width - ctrl->blcolmin;
            assert( lwidth > ctrl->blcolmin );

            lcol = fcol + split_fct( nblocksperline + fcol,
                                     step, ctrl->blcolmin, lwidth,
                                     authorized_percent );

            assert( (lcol >= fcol) && (lcol <= lcolnum) );

            extraCblkAdd( extracblk, fcol, lcol,
                          curcblk->selevtx );
            nbcblk++;

            lwidth = lcol - fcol + 1;
            width  = width - lwidth;
            fcol   = lcol + 1;
        }

        /*
         * Mark the cblk as being splitted
         */
        extracblk->addcblk += nbcblk-1;
        extracblk->sptcblk[cblknum] = extracblk->curcblk - nbcblk + 1;
        extracblk->sptcbnb[cblknum] = nbcblk;

        /* Update the number of blocks per line */
        curblok = &(symbmtx->bloktab[curcblk->bloknum + 1]) ;
        for(bloknum = curcblk[0].bloknum + 1;
            bloknum < curcblk[1].bloknum; bloknum++, curblok++)
        {
            for(line = curblok->frownum; line < curblok->lrownum; line++ )
            {
                nblocksperline[ line ] += nbcblk-1;
            }
        }
    }

    if ( fsplitrow != -1) {
        nblocksperline += fsplitrow;
        memFree_null( nblocksperline );
    }
    else {
        assert( nblocksperline == NULL );
    }
}
/**
 *@}
 */

/**
 *******************************************************************************
 *
 * @ingroup pastix_blend
 *
 * @brief Split the column blocks of the symbol matrix to generate parallelism.
 *
 * This is the main function that cut the symbol matrix column blocks, and
 * return the new symbolMatrix. Cost matrix, elimination tree, and candidate
 * array are updated on exit of this function with respect to the newly created
 * column blocks and blocks. See splitSmart() for the cutting algorithm.
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The blend control structure.
 *          On entry, candtab must be initialized.
 *          On exit, costmtx, etree, and candtab are updated accordinglyy to the
 *          extended symbol matrix, if new cblk are generated.
 *
 * @param[in] symbmtx
 *          On entry, the symbol matrix structure to split.
 *          On exit, the new symbol matrix with the new cblk and blocks.
 *
 *******************************************************************************/
void
splitSymbol( BlendCtrl       *ctrl,
             symbol_matrix_t *symbmtx )
{
    ExtraCblk_t extracblk;

    /* Init structure to store extra cblks */
    extraCblkInit( symbmtx->cblknbr, &extracblk );

    splitSmart( ctrl, symbmtx, &extracblk );

    /* Rk: addcblk field is not erased by Exit call, so we can freely use it */
    if ( extracblk.addcblk )
    {
        /* Merge the initial matrix and the newly generated cblks */
        extraCblkMerge( &extracblk, symbmtx, &(ctrl->candtab) );
        extraCblkExit( &extracblk );

        /* Check that the generated symbol matrix is correct */
        if (ctrl->debug) {
            pastixSymbolCheck(symbmtx);
        }

        if ( ctrl->up_after_split ) {
            /* Update cost matrix to fill-in blank of newly generated blocks */
            costMatrixExit(ctrl->costmtx);
            memFree_null(ctrl->costmtx);
            ctrl->costmtx = costMatrixBuild( symbmtx,
                                             ctrl->iparm[IPARM_FLOAT],
                                             ctrl->iparm[IPARM_FACTORIZATION] );

            /* Update elimination tree */
            if (ctrl->etree != NULL) {
                eTreeExit(ctrl->etree);
            }
            ctrl->etree = eTreeBuild(symbmtx);

            /*
             * Let's update cost in the candtab for the proportionnal mapping and
             * the simulation
             */
            candUpdate( ctrl->candtab, ctrl->etree, symbmtx, ctrl->costmtx );
        }
    }

    if ( ctrl->clustnum == 0 ) {
        if (ctrl->iparm[IPARM_VERBOSE] > PastixVerboseNo) {
            pastixSymbolPrintStats( symbmtx );
        }
    }
}
