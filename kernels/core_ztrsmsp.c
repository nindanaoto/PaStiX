/**
 *
 * @file core_ztrsmsp.c
 *
 * PaStiX kernel routines
 *
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2020-02-05
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "cblas.h"
#include "blend/solver.h"
#include "kernels_trace.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static pastix_complex64_t zone =  1.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup kernel_fact_null
 *
 * @brief Apply all the trsm updates on a panel stored in 1D layout.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specify whether the A matrix appears on the left or right in the
 *          equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the A matrix is upper or lower triangular. It has to
 *          be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the A matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the A matrix is unit triangular. It has to be either
 *          PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and C pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 *******************************************************************************/
static inline void
core_ztrsmsp_1d( pastix_side_t             side,
                 pastix_uplo_t             uplo,
                 pastix_trans_t            trans,
                 pastix_diag_t             diag,
                 SolverCblk               *cblk,
                 const pastix_complex64_t *A,
                 pastix_complex64_t       *C )
{
    SolverBlok *fblok;
    pastix_int_t M, N, lda;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    lda   = cblk->stride;
    fblok = cblk->fblokptr;  /* The diagonal block */

    /* vertical dimension */
    M = lda - N;

    /* if there is an extra-diagonal bloc in column block */
    assert( fblok + 1 < cblk[1].fblokptr );
    assert( blok_rownbr( fblok) == N );
    assert(!(cblk->cblktype & CBLK_LAYOUT_2D));

    /* first extra-diagonal bloc in column block address */
    C = C + fblok[1].coefind;

    kernel_trace_start_lvl2( PastixKernelLvl2_FR_TRSM );
    cblas_ztrsm(CblasColMajor,
                (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                M, N,
                CBLAS_SADDR(zone), A, lda,
                                   C, lda);
    kernel_trace_stop_lvl2( FLOPS_ZTRSM( side, M, N ) );
}

/**
 *******************************************************************************
 *
 * @ingroup kernel_fact_null
 *
 * @brief Compute the updates associated to one off-diagonal block between two
 * cblk stored in 2D.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specify whether the A matrix appears on the left or right in the
 *          equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the A matrix is upper or lower triangular. It has to
 *          be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the A matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the A matrix is unit triangular. It has to be either
 *          PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and C pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 *******************************************************************************/
static inline void
core_ztrsmsp_2d( pastix_side_t             side,
                 pastix_uplo_t             uplo,
                 pastix_trans_t            trans,
                 pastix_diag_t             diag,
                 SolverCblk               *cblk,
                 const pastix_complex64_t *A,
                 pastix_complex64_t       *C )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda, ldc;
    pastix_complex64_t *blokC;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */
    lda   = blok_rownbr( fblok );

    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_LAYOUT_2D );

    for (blok=fblok+1; blok<lblok; blok++) {

        blokC = C + blok->coefind;
        M   = blok_rownbr(blok);
        ldc = M;

        kernel_trace_start_lvl2( PastixKernelLvl2_FR_TRSM );
        cblas_ztrsm(CblasColMajor,
                    (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                    M, N,
                    CBLAS_SADDR(zone), A, lda,
                                       blokC, ldc);
        kernel_trace_stop_lvl2( FLOPS_ZTRSM( side, M, N ) );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup kernel_fact_null
 *
 * @brief Computes the updates associated to one off-diagonal block between two
 * cblk stored in low-rank format.
 *
 *******************************************************************************
 *
 * @param[in] coef
 *          - PastixLCoef, use the lower part of the off-diagonal blocks.
 *          - PastixUCoef, use the upper part of the off-diagonal blocks
 *
 * @param[in] side
 *          Specify whether the off-diagonal blocks appear on the left or right in the
 *          equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the off-diagonal blocks are upper or lower
 *          triangular. It has to be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the off-diagonal blocks. It has
 *          to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and C pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] solvmtx
 *          The symbolic structure of pastix.
 *
 *******************************************************************************
 *
 * @return  The number of flops performed
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
core_ztrsmsp_lr( pastix_coefside_t coef, pastix_side_t side, pastix_uplo_t uplo,
                 pastix_trans_t trans, pastix_diag_t diag,
                 SolverCblk *cblk, SolverMatrix *solvmtx )
{
    SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda;
    pastix_lrblock_t *lrA, *lrC;
    pastix_complex64_t *A;

    pastix_fixdbl_t flops = 0.0;
    pastix_fixdbl_t flops_lr, flops_c;

    pastix_lr_t *lowrank = &solvmtx->lowrank;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */

    lrA   = fblok->LRblock + coef;
    A     = lrA->u;
    lda   = lrA->rkmax;

    assert( lrA->rk == -1 );
    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_COMPRESSED );
    assert( cblk->cblktype & CBLK_LAYOUT_2D  );

    for (blok=fblok+1; blok<lblok; blok++) {

        M   = blok_rownbr(blok);
        lrC = blok->LRblock + coef;
        flops_lr = 0.;
        flops_c  = 0.;

        /* Check the size of the block */
        if ( ( N >= lowrank->compress_min_width  ) &&
             ( M >= lowrank->compress_min_height ) )
        {
            int is_preselected = blok_is_preselected( cblk, blok, solvmtx->cblktab + blok->fcblknm );

            /*
             * Try to compress the block: 2 cases
             *   - We are in the compress_end version
             *   - We are in the compress_begin version, and the block was preselected
             */
            if ( (lowrank->compress_when == PastixCompressWhenEnd) &&
                 (lowrank->compress_preselect || (!is_preselected)) )
            {
                flops_lr = cpublok_zcompress( lowrank, coef, M, N, blok );
            }
            if ( (lowrank->compress_when == PastixCompressWhenBegin) &&
                 (lowrank->compress_preselect && is_preselected ) )
            {
                flops_lr = cpublok_zcompress( lowrank, coef, M, N, blok );
            }
        }

        if ( lrC->rk != 0 ) {
            if ( lrC->rk != -1 ) {
                kernel_trace_start_lvl2( PastixKernelLvl2_LR_TRSM );
                cblas_ztrsm(CblasColMajor,
                            (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                            lrC->rk, N,
                            CBLAS_SADDR(zone), A, lda,
                            lrC->v, lrC->rkmax);
                flops_c = FLOPS_ZTRSM( side, lrC->rk, N );
                kernel_trace_stop_lvl2( flops_c );
            }
            else {
                kernel_trace_start_lvl2( PastixKernelLvl2_FR_TRSM );
                cblas_ztrsm(CblasColMajor,
                            (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                            M, N,
                            CBLAS_SADDR(zone), A, lda,
                            lrC->u, lrC->rkmax);
                flops_c = FLOPS_ZTRSM( side, M, N );
                kernel_trace_stop_lvl2( flops_c );
            }
        }

        flops += flops_lr + flops_c;
    }
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Compute the updates associated to a column of off-diagonal blocks.
 *
 *******************************************************************************
 *
 * @param[in] coef
 *          Specify whether we work with the lower matrix, or the upper matrix.
 *
 * @param[in] side
 *          Specify whether the A matrix appears on the left or right in the
 *          equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the A matrix is upper or lower triangular. It has to
 *          be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the A matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the A matrix is unit triangular. It has to be either
 *          PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 * @param[in] solvmtx
 *          The symbolic structure of pastix.
 *
 *******************************************************************************/
void
cpucblk_ztrsmsp( pastix_coefside_t coef, pastix_side_t side, pastix_uplo_t uplo,
                 pastix_trans_t trans, pastix_diag_t diag,
                       SolverCblk         *cblk,
                 const pastix_complex64_t *A,
                       pastix_complex64_t *C,
                       SolverMatrix       *solvmtx )
{
    if (  cblk[0].fblokptr + 1 < cblk[1].fblokptr )
    {
        pastix_ktype_t ktype = PastixKernelLvl1Nbr;
        pastix_fixdbl_t time, flops = 0.0;
        pastix_int_t n = cblk_colnbr( cblk );
        pastix_int_t m = cblk->stride - n;

        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            ktype = PastixKernelTRSMCblkLR;
            time  = kernel_trace_start( ktype );

            flops = core_ztrsmsp_lr( coef, side, uplo, trans, diag,
                                     cblk, solvmtx );
        }
        else {
            if ( cblk->cblktype & CBLK_LAYOUT_2D ) {
                ktype = PastixKernelTRSMCblk2d;
                time  = kernel_trace_start( ktype );

                core_ztrsmsp_2d( side, uplo, trans, diag,
                                 cblk, A, C );
            }
            else {
                ktype = PastixKernelTRSMCblk1d;
                time  = kernel_trace_start( ktype );

                core_ztrsmsp_1d( side, uplo, trans, diag,
                                 cblk, A, C );
            }
            flops = FLOPS_ZTRSM( PastixRight, m, n );
        }

        kernel_trace_stop( cblk->fblokptr->inlast, ktype, m, n, 0, flops, time );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup kernel_fact_null
 *
 * @brief Compute the updates associated to one off-diagonal block between two
 * cblk stored in 2D.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specify whether the A matrix appears on the left or right in the
 *          equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the A matrix is upper or lower triangular. It has to
 *          be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the A matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the A matrix is unit triangular. It has to be either
 *          PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and C pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok_m
 *          Index of the first off-diagonal block in cblk that is solved. The
 *          TRSM is also applied to all the folowing blocks which are facing the
 *          same diagonal block
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 *******************************************************************************/
static inline void
core_ztrsmsp_2dsub( pastix_side_t side, pastix_uplo_t uplo,
                    pastix_trans_t trans, pastix_diag_t diag,
                          SolverCblk         *cblk,
                          pastix_int_t        blok_m,
                    const pastix_complex64_t *A,
                          pastix_complex64_t *C )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda, ldc, offset, cblk_m, full_m;
    pastix_complex64_t *Cptr;
    pastix_fixdbl_t flops = 0.0;
    pastix_fixdbl_t time = kernel_trace_start( PastixKernelTRSMBlok2d );

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */
    lda   = blok_rownbr( fblok );

    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_LAYOUT_2D );

    blok   = fblok + blok_m;
    offset = blok->coefind;
    cblk_m = blok->fcblknm;
    full_m = 0;

    for (; (blok < lblok) && (blok->fcblknm == cblk_m); blok++) {

        Cptr = C + blok->coefind - offset;
        M   = blok_rownbr(blok);
        ldc = M;

        cblas_ztrsm( CblasColMajor,
                     (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                     M, N,
                     CBLAS_SADDR(zone), A, lda,
                                        Cptr, ldc );

        flops += FLOPS_ZTRSM( side, M, N );
        full_m += M;
    }

    kernel_trace_stop( cblk->fblokptr->inlast, PastixKernelTRSMBlok2d,
                       full_m, N, 0, flops, time );
    return;
}

/**
 *******************************************************************************
 *
 * @ingroup kernel_fact_null
 *
 * @brief Compute the updates associated to one off-diagonal block between two
 * cblk stored in low-rank format.
 *
 *******************************************************************************
 *
 * @param[in] coef
 *          - PastixLCoef, use the lower part of the off-diagonal blocks.
 *          - PastixUCoef, use the upper part of the off-diagonal blocks
 *
 * @param[in] side
 *          Specify whether the off-diagonal blocks appear on the left or right in the
 *          equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the off-diagonal blocks are upper or lower
 *          triangular. It has to be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the off-diagonal blocks. It has
 *          to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and C pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok_m
 *          Index of the first off-diagonal block in cblk that is solved. The
 *          TRSM is also applied to all the folowing blocks which are facing the
 *          same diagonal block
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 *******************************************************************************/
static inline void
core_ztrsmsp_lrsub( pastix_coefside_t   coef,
                    pastix_side_t       side,
                    pastix_uplo_t       uplo,
                    pastix_trans_t      trans,
                    pastix_diag_t       diag,
                    SolverCblk         *cblk,
                    pastix_int_t        blok_m,
                    const pastix_lr_t  *lowrank )
{
    SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda, cblk_m, full_m, full_n;
    pastix_complex64_t *A;
    pastix_lrblock_t *lrA, *lrC;
    pastix_fixdbl_t flops = 0.0;
    pastix_fixdbl_t time = kernel_trace_start( PastixKernelTRSMBlokLR );

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */

    lrA   = fblok->LRblock + coef;
    A     = lrA->u;
    lda   = lrA->rkmax;

    assert( cblk->cblktype & CBLK_COMPRESSED );
    assert( cblk->cblktype & CBLK_LAYOUT_2D  );

    assert( blok_rownbr(fblok) == N );
    assert( lrA->rk == -1 );

    blok   = fblok + blok_m;
    cblk_m = blok->fcblknm;
    full_m = 0;
    full_n = 0;

    for (; (blok < lblok) && (blok->fcblknm == cblk_m); blok++) {

        M = blok_rownbr(blok);
        lrC = blok->LRblock + coef;

        if ( ( N >= lowrank->compress_min_width ) &&
             ( M >= lowrank->compress_min_height ) )
        {
            int is_preselected = blok_is_preselected( cblk, blok, cblk + (blok->fcblknm - blok->lcblknm) );

            /*
             * Try to compress the block: 2 cases
             *   - We are in the compress_end version
             *   - We are in the compress_begin version, and the block was preselected
             */
            if ( (lowrank->compress_when == PastixCompressWhenEnd) &&
                 (lowrank->compress_preselect || (!is_preselected)) )
            {
                flops += cpublok_zcompress( lowrank, coef, M, N, blok );
            }
            if ( (lowrank->compress_when == PastixCompressWhenBegin) &&
                 (lowrank->compress_preselect &&  is_preselected ) )
            {
                flops += cpublok_zcompress( lowrank, coef, M, N, blok );
            }
        }

        if ( lrC->rk != 0 ) {
            if ( lrC->rk != -1 ) {
                cblas_ztrsm(CblasColMajor,
                            (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                            lrC->rk, N,
                            CBLAS_SADDR(zone), A, lda,
                            lrC->v, lrC->rkmax);

                flops += FLOPS_ZTRSM( side, lrC->rk, N );
                full_n += lrC->rk;
            }
            else {
                cblas_ztrsm(CblasColMajor,
                            (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                            M, N,
                            CBLAS_SADDR(zone), A, lda,
                            lrC->u, lrC->rkmax);

                flops += FLOPS_ZTRSM( side, M, N );
                full_n += M;
            }
        }
        full_m += M;
    }

    kernel_trace_stop( cblk->fblokptr->inlast, PastixKernelTRSMBlokLR,
                       full_m, N, full_n, flops, time );
}

/**
 *******************************************************************************
 *
 * @brief Compute the updates associated to one off-diagonal block.
 *
 *******************************************************************************
 *
 * @param[in] coef
 *          - PastixLCoef, use the lower part of the off-diagonal blocks.
 *          - PastixUCoef, use the upper part of the off-diagonal blocks
 *
 * @param[in] side
 *          Specify whether the A matrix appears on the left or right in the
 *          equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the A matrix is upper or lower triangular. It has to
 *          be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the A matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the A matrix is unit triangular. It has to be either
 *          PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok_m
 *          Index of the first off-diagonal block in cblk that is solved. The
 *          TRSM is also applied to all the folowing blocks which are facing the
 *          same diagonal block
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 *******************************************************************************/
void
cpublok_ztrsmsp( pastix_coefside_t coef, pastix_side_t side, pastix_uplo_t uplo,
                 pastix_trans_t trans, pastix_diag_t diag,
                       SolverCblk         *cblk,
                       pastix_int_t        blok_m,
                 const pastix_complex64_t *A,
                       pastix_complex64_t *C,
                 const pastix_lr_t        *lowrank )
{
    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        core_ztrsmsp_lrsub( coef, side, uplo, trans, diag,
                            cblk, blok_m, lowrank );
    }
    else {
        core_ztrsmsp_2dsub( side, uplo, trans, diag,
                            cblk, blok_m, A, C );
    }
}
