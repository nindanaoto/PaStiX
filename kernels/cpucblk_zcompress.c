/**
 *
 * @file cpucblk_zcompress.c
 *
 * Precision dependent function to compress/uncompress the coefficients
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2020-01-29
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include <lapacke.h>
#include "kernels_trace.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

/**
 *******************************************************************************
 *
 * @brief Compress a single block from full-rank to low-rank format
 *
 * The compression to low-rank format is parameterized by the input information
 * stored in the low-rank structure.
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The pointer to the low-rank structure describing the lo-rank
 *          compression parameters.
 *
 * @param[in] side
 *          Define which side of the matrix must be compressed.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *
 * @param[in] M
 *          The number of rows in the block
 *
 * @param[in] N
 *          The number of columns in the block
 *
 * @param[inout] blok
 *          The block to compress. On input, it points to a full-rank matrix. On
 *          output, if possible the matrix is compressed in block low-rank
 *          format.
 *
 *******************************************************************************
 *
 * @return The number of flops used to compress the block.
 *
 *******************************************************************************/
pastix_fixdbl_t
cpublok_zcompress( const pastix_lr_t *lowrank,
                   pastix_coefside_t side,
                   pastix_int_t M, pastix_int_t N,
                   SolverBlok *blok )
{
    pastix_fixdbl_t flops;
    pastix_lrblock_t   *lrA = blok->LRblock + side;
    pastix_complex64_t *A = lrA->u;

    assert( lrA->rk == -1   );
    assert( lrA->u  != NULL );
    assert( lrA->v  == NULL );

    kernel_trace_start_lvl2( PastixKernelLvl2_LR_init_compress );
    flops = lowrank->core_ge2lr( lowrank->use_reltol, lowrank->tolerance, -1,
                                 M, N, A, M, lrA );
    kernel_trace_stop_lvl2_rank( flops, lrA->rk );

    free( A );

    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Compress a single column block from full-rank to low-rank format
 *
 * The compression to low-rank format is parameterized by the input information
 * stored in the low-rank structure.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The pointer to the solver structure.
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] cblk
 *          The column block to compress.
 *
 *******************************************************************************
 *
 * @return The memory gain resulting from the compression to low-rank format in
 *         number of elements.
 *
 *******************************************************************************/
pastix_int_t
cpucblk_zcompress( const SolverMatrix *solvmtx,
                   pastix_coefside_t   side,
                   SolverCblk         *cblk )
{
    pastix_lrblock_t   *lrA;
    SolverBlok         *blok  = cblk[0].fblokptr + 1;
    SolverBlok         *lblok = cblk[1].fblokptr;
    pastix_int_t        ncols = cblk_colnbr( cblk );
    pastix_int_t        gain;
    pastix_int_t        gainL = 0;
    pastix_int_t        gainU = 0;
    const pastix_lr_t  *lowrank = &(solvmtx->lowrank);

    assert( cblk->cblktype & CBLK_LAYOUT_2D  );
    assert( cblk->cblktype & CBLK_COMPRESSED );

    if ( ncols < lowrank->compress_min_width ) {
        return 0;
    }

    for (; blok<lblok; blok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );
        int is_preselected = blok_is_preselected( cblk, blok, solvmtx->cblktab + blok->fcblknm );

        /*
         * If we are in the 'end' scenario, the preselected blocks must be
         * compressed.
         * If we are in the 'begin' or 'during' scenarii, the preselected blocks
         * may be compressed directly by calling cpublok_zcompress and not
         * cpucblk_zcompress.
         */
        if ( (lowrank->compress_when == PastixCompressWhenEnd) &&
             (lowrank->compress_preselect) )
        {
            is_preselected = 0;
        }

        /* Skip uncompressible blocks */
        if ( nrows < lowrank->compress_min_height ) {
            continue;
        }

        if ( is_preselected ) {
            continue;
        }

        gain = nrows * ncols;

        /* Lower part */
        if ( side != PastixUCoef ) {
            lrA = blok->LRblock;

            /* Try to compress non selected blocks */
            if ( lrA->rk == -1 ) {
                cpublok_zcompress( lowrank, PastixLCoef, nrows, ncols, blok );
            }

            if ( lrA->rk != -1 ) {
                gainL += gain - ((nrows+ncols) * lrA->rk);
            }
        }

        /* Upper part */
        if ( side != PastixLCoef ) {
            lrA = blok->LRblock + 1;

            if ( lrA->rk == -1 ) {
                cpublok_zcompress( lowrank, PastixUCoef, nrows, ncols, blok );
            }

            if ( lrA->rk != -1 ) {
                gainU += gain - ((nrows+ncols) * lrA->rk);
            }
        }
    }

    return gainL + gainU;
}

/**
 *******************************************************************************
 *
 * @brief Uncompress a single column block from low-rank format to full-rank
 * format.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] cblk
 *          The column block to uncompress.
 *
 *******************************************************************************/
void
cpucblk_zuncompress( pastix_coefside_t side,
                     SolverCblk       *cblk )
{
    SolverBlok  *blok, *lblok;
    pastix_int_t ncols = cblk_colnbr( cblk );
    int ret;

    if ( side != PastixUCoef ) {
        blok  = cblk[0].fblokptr;
        lblok = cblk[1].fblokptr;
        for (; blok<lblok; blok++)
        {
            pastix_lrblock_t lrtmp;
            pastix_int_t nrows = blok_rownbr( blok );

            memcpy( &lrtmp, blok->LRblock, sizeof(pastix_lrblock_t) );

            core_zlralloc( nrows, ncols, -1, blok->LRblock );
            ret = core_zlr2ge( PastixNoTrans, nrows, ncols,
                               &lrtmp,
                               blok->LRblock[0].u, nrows );
            assert( ret == 0 );
            core_zlrfree( &lrtmp );
        }
    }

    if ( side != PastixLCoef ) {
        blok  = cblk[0].fblokptr;
        lblok = cblk[1].fblokptr;
        for (; blok<lblok; blok++)
        {
            pastix_lrblock_t lrtmp;
            pastix_int_t nrows = blok_rownbr( blok );

            memcpy( &lrtmp, blok->LRblock+1, sizeof(pastix_lrblock_t) );

            core_zlralloc( nrows, ncols, -1, blok->LRblock+1 );
            ret = core_zlr2ge( PastixNoTrans, nrows, ncols,
                               &lrtmp,
                               blok->LRblock[1].u, nrows );
            assert( ret == 0 );
            core_zlrfree( &lrtmp );
        }
    }

    (void)ret;
}

/**
 *******************************************************************************
 *
 * @brief Return the memory gain of the low-rank form over the full-rank form
 * for a single column-block.
 *
 * This function returns the memory gain in number of elements for a single
 * column block when it is stored in low-rank format compared to a full rank
 * storage.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] solvmtx
 *          The pointer to the solver structure.
 *
 * @param[in] cblk
 *          The column block to study.
 *
 * @param[in,out] orig
 *          The structure that counts the original cost of the blocks.
 *
 * @param[in,out] gain
 *          The structure that counts gain on each type of the blocks.
 *
 *******************************************************************************/
void
cpucblk_zmemory( pastix_coefside_t   side,
                 SolverMatrix       *solvmtx,
                 SolverCblk         *cblk,
                 pastix_int_t       *orig,
                 pastix_int_t       *gain )
{
    SolverBlok *blok  = cblk[0].fblokptr + 1;
    SolverBlok *lblok = cblk[1].fblokptr;

    pastix_int_t ncols = cblk_colnbr( cblk );
    pastix_int_t size;
    pastix_int_t origblok;
    pastix_int_t gainblok, gaintmp;

    /* Compute potential gains if blocks where not compressed */
    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        cpucblk_zcompress( solvmtx, side, cblk );
    }

    for (; blok<lblok; blok++)
    {
        const SolverCblk *fcblk = solvmtx->cblktab + blok->fcblknm;
        pastix_int_t nrows = blok_rownbr( blok );
        size = nrows * ncols;
        gainblok = 0;
        origblok = size;

        /* Lower part */
        if ( (side != PastixUCoef) &&
             (blok->LRblock[0].rk >= 0) )
        {
            gaintmp = (size - ((nrows+ncols) * blok->LRblock[0].rkmax));
            assert( gaintmp >= 0 );
            gainblok += gaintmp;
        }

        /* Upper part */
        if ( (side != PastixLCoef) &&
             (blok->LRblock[1].rk >= 0) )
        {
            gaintmp = (size - ((nrows+ncols) * blok->LRblock[1].rkmax));
            assert( gaintmp >= 0 );
            gainblok += gaintmp;
        }

        if ( blok_is_preselected( cblk, blok, fcblk ) )
        {
            /* Selected block should always be inside supernode diagonal blocks */
            assert( fcblk->sndeidx == cblk->sndeidx );
            gain[LR_InSele] += gainblok;
            orig[LR_InSele] += origblok;
        }
        else{
            if ( fcblk->sndeidx == cblk->sndeidx ) {
                gain[LR_InDiag] += gainblok;
                orig[LR_InDiag] += origblok;
            }
            else {
                gain[LR_OffDiag] += gainblok;
                orig[LR_OffDiag] += origblok;
            }
        }
    }
    return;
}
