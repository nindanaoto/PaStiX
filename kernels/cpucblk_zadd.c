/**
 *
 * @file cpucblk_zadd.c
 *
 * Precision dependent routines to add different cblks.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2020-01-26
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include "kernels_trace.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

/**
 *******************************************************************************
 *
 * @brief Add two column bloks in full rank format.
 *
 * The second cblk is overwritten by the sum of the two column blocks.
 *              B <- alpha * A + B
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part
 *          @arg PastixUCoef if upper part
 *
 * @param[in] alpha
 *          The scalar alpha
 *
 * @param[in] cblkA
 *          The column block of the A matrix.
 *
 * @param[inout] cblkB
 *          The column block of the B matrix
 *          On exit, cblkB coefficient arrays are overwritten by the result of
 *          alpha * A + B.
 *
 *******************************************************************************
 *
 * @return The number of flops of the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
cpucblk_zadd_frlr( pastix_coefside_t   side,
                   pastix_int_t        alpha,
                   const SolverCblk   *cblkA,
                   SolverCblk         *cblkB,
                   pastix_complex64_t *work,
                   pastix_int_t        lwork,
                   const pastix_lr_t  *lowrank )
{
    const SolverBlok   *blokA  = cblkA->fblokptr;
    const SolverBlok   *blokB  = cblkB->fblokptr;
    const SolverBlok   *lblokA = cblkA[1].fblokptr;
    const SolverBlok   *lblokB = cblkB[1].fblokptr;
    pastix_complex64_t *A;
    pastix_int_t        shift;
    pastix_fixdbl_t     flops = 0.;
    core_zlrmm_t        params;
    pastix_lrblock_t    lrA;

    assert( !(cblkA->cblktype & CBLK_COMPRESSED) );
    assert(   cblkB->cblktype & CBLK_COMPRESSED  );
    assert(   cblkA->cblktype & CBLK_LAYOUT_2D   );

    if ( side == PastixUCoef ) {
        A = cblkA->ucoeftab;
        shift = 1;
    }
    else {
        A = cblkA->lcoeftab;
        shift = 0;
    }

    assert( A != NULL );

    params.lowrank = lowrank;
    params.transA  = PastixNoTrans; /* Unused */
    params.transB  = PastixNoTrans; /* Unused */
    params.K       = -1;            /* Unused */
    params.alpha   = alpha;
    params.A       = NULL;          /* Unused */
    params.B       = NULL;          /* Unused */
    params.beta    = 1.0;
    params.work    = work;
    params.lwork   = lwork;
    params.lwused  = 0;
    params.lock    = &(cblkB->lock);

    /* Dimensions on N */
    params.N    = cblk_colnbr( cblkA );
    params.Cn   = cblk_colnbr( cblkB );
    params.offy = cblkA->fcolnum - cblkB->fcolnum;

    lrA.rk = -1;
    lrA.v  = NULL;

    for (; blokA < lblokA; blokA++) {

        /* Find facing bloknum */
        while ( !is_block_inside_fblock( blokA, blokB ) && (blokB < lblokB) ) {
            blokB++;
        }

        assert( is_block_inside_fblock( blokA, blokB ) && (blokB <= lblokB) );

        lrA.u     = A + blokA->coefind;
        lrA.rkmax = blok_rownbr( blokA );

        /* Dimensions on M */
        params.M    = blok_rownbr( blokA );
        params.Cm   = blok_rownbr( blokB );
        params.offx = blokA->frownum - blokB->frownum;
        params.C    = blokB->LRblock + shift;

        flops += core_zlradd( &params, &lrA,
                              PastixNoTrans, 0 );
    }
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Add two column bloks in full rank format.
 *
 * The second cblk is overwritten by the sum of the two column blocks.
 *              B <- alpha * A + B
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part
 *          @arg PastixUCoef if upper part
 *
 * @param[in] alpha
 *          The scalar alpha
 *
 * @param[in] cblkA
 *          The column block of the A matrix.
 *
 * @param[inout] cblkB
 *          The column block of the B matrix
 *          On exit, cblkB coefficient arrays are overwritten by the result of
 *          alpha * A + B.
 *
 *******************************************************************************
 *
 * @return The number of flops of the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
cpucblk_zadd_lrlr( pastix_coefside_t   side,
                   pastix_int_t        alpha,
                   const SolverCblk   *cblkA,
                   SolverCblk         *cblkB,
                   pastix_complex64_t *work,
                   pastix_int_t        lwork,
                   const pastix_lr_t  *lowrank )
{
    const SolverBlok   *blokA  = cblkA->fblokptr;
    const SolverBlok   *blokB  = cblkB->fblokptr;
    const SolverBlok   *lblokA = cblkA[1].fblokptr;
    const SolverBlok   *lblokB = cblkB[1].fblokptr;
    pastix_int_t        shift;
    pastix_fixdbl_t     flops = 0.;
    core_zlrmm_t        params;

    assert( (cblkA->cblktype & CBLK_COMPRESSED) );
    assert( (cblkB->cblktype & CBLK_COMPRESSED) );

    shift = (side == PastixUCoef) ? 1 : 0;

    params.lowrank = lowrank;
    params.transA  = PastixNoTrans; /* Unused */
    params.transB  = PastixNoTrans; /* Unused */
    params.K       = -1;            /* Unused */
    params.alpha   = alpha;
    params.A       = NULL;          /* Unused */
    params.B       = NULL;          /* Unused */
    params.beta    = 1.0;
    params.work    = work;
    params.lwork   = lwork;
    params.lwused  = 0;
    params.lock    = &(cblkB->lock);

    /* Dimensions on N */
    params.N    = cblk_colnbr( cblkA );
    params.Cn   = cblk_colnbr( cblkB );
    params.offy = cblkA->fcolnum - cblkB->fcolnum;

    for (; blokA < lblokA; blokA++) {

        /* Find facing bloknum */
        while ( !is_block_inside_fblock( blokA, blokB ) && (blokB < lblokB) ) {
            blokB++;
        }

        assert( is_block_inside_fblock( blokA, blokB ) && (blokB <= lblokB) );

        /* Dimensions on M */
        params.M    = blok_rownbr( blokA );
        params.Cm   = blok_rownbr( blokB );
        params.offx = blokA->frownum - blokB->frownum;
        params.C    = blokB->LRblock + shift;

        flops += core_zlradd( &params, blokA->LRblock + shift,
                              PastixNoTrans, PASTIX_LRM3_ORTHOU );
    }
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Add two column bloks in full rank format.
 *
 * The second cblk is overwritten by the sum of the two column blocks.
 *              B <- alpha * A + B
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part
 *          @arg PastixUCoef if upper part
 *
 * @param[in] alpha
 *          The scalar alpha
 *
 * @param[in] cblkA
 *          The column block of the A matrix.
 *
 * @param[inout] cblkB
 *          The column block of the B matrix
 *          On exit, cblkB coefficient arrays are overwritten by the result of
 *          alpha * A + B.
 *
 *******************************************************************************
 *
 * @return The number of flops of the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
cpucblk_zadd_frfr( pastix_coefside_t side,
                   pastix_int_t      alpha,
                   const SolverCblk *cblkA,
                   SolverCblk       *cblkB )
{
    pastix_complex64_t *A, *B;
    pastix_int_t        n = cblk_colnbr( cblkA );
    pastix_int_t        m = cblkA->stride;
    pastix_fixdbl_t     flops = m * n;

    assert( !(cblkA->cblktype & CBLK_COMPRESSED) );
    assert( !(cblkB->cblktype & CBLK_COMPRESSED) );

    if ( side == PastixUCoef ) {
        A = cblkA->ucoeftab;
        B = cblkB->ucoeftab;
    }
    else {
        A = cblkA->lcoeftab;
        B = cblkB->lcoeftab;
    }

    assert( (A != NULL) && (B != NULL) );

    /* If the cblk matches */
    if ( (n == cblk_colnbr( cblkB )) &&
         (m == cblkB->stride) ) {

        pastix_cblk_lock( cblkB );
        core_zgeadd( PastixNoTrans, m, n,
                     alpha, A, m,
                        1., B, m );
        pastix_cblk_unlock( cblkB );
    }
    else {
        pastix_complex64_t *bA, *bB;
        const SolverBlok   *blokA  = cblkA->fblokptr;
        const SolverBlok   *blokB  = cblkB->fblokptr;
        const SolverBlok   *lblokA = cblkA[1].fblokptr;
        const SolverBlok   *lblokB = cblkB[1].fblokptr;
        pastix_int_t        lda, ldb;

        /* Both cblk A and B must be stored in 2D */
        assert( cblkA->cblktype & CBLK_LAYOUT_2D );
        assert( cblkB->cblktype & CBLK_LAYOUT_2D );

        for (; blokA < lblokA; blokA++) {

            /* Find facing bloknum */
            while ( !is_block_inside_fblock( blokA, blokB ) && (blokB < lblokB) ) {
                blokB++;
            }

            assert( is_block_inside_fblock( blokA, blokB ) && (blokB <= lblokB) );

            bA = A + blokA->coefind;
            bB = B + blokB->coefind;
            lda = blok_rownbr( blokA );
            ldb = blok_rownbr( blokB );

            bB = bB + ldb * ( cblkA->fcolnum - cblkB->fcolnum ) + ( blokA->frownum - blokB->frownum );
            m = lda;

            pastix_cblk_lock( cblkB );
            core_zgeadd( PastixNoTrans, m, n,
                         alpha, bA, lda,
                            1., bB, ldb );
            pastix_cblk_unlock( cblkB );
        }
    }
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Add two column bloks in full rank format.
 *
 * The second cblk is overwritten by the sum of the two column blocks.
 *              B <- alpha * A + B
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *
 * @param[in] alpha
 *          The scalar alpha
 *
 * @param[in] cblkA
 *          The column block of the A matrix.
 *
 * @param[inout] cblkB
 *          The column block of the B matrix
 *          On exit, cblkB coefficient arrays are overwritten by the result of
 *          alpha * A + B.
 *
 *******************************************************************************/
void
cpucblk_zadd( pastix_coefside_t  side,
              double             alpha,
              const SolverCblk  *cblkA,
              SolverCblk        *cblkB,
              const pastix_lr_t *lowrank )
{
    pastix_ktype_t ktype = PastixKernelGEADDCblkFRFR;
    pastix_fixdbl_t time, flops = 0.0;
    pastix_int_t m = cblkA->stride;
    pastix_int_t n = cblk_colnbr( cblkA );

    if ( side == PastixLUCoef ) {
        n *= 2;
    }

    if ( cblkB->cblktype & CBLK_COMPRESSED ) {
        if ( cblkA->cblktype & CBLK_COMPRESSED ) {
            ktype = PastixKernelGEADDCblkLRLR;
            time  = kernel_trace_start( ktype );
            flops = cpucblk_zadd_lrlr( side, alpha, cblkA, cblkB,
                                       NULL, 0, lowrank );
        }
        else {
            ktype = PastixKernelGEADDCblkFRLR;
            time  = kernel_trace_start( ktype );
            flops = cpucblk_zadd_frlr( side, alpha, cblkA, cblkB,
                                       NULL, 0, lowrank );
        }
    }
    else {
        if ( cblkA->cblktype & CBLK_COMPRESSED ) {
            assert(0 /* We do not add a compressed cblk to a non compressed cblk */);
            time  = kernel_trace_start( ktype );
        }
        else {
            ktype = PastixKernelGEADDCblkFRFR;
            time  = kernel_trace_start( ktype );
            flops = cpucblk_zadd_frfr( side, alpha, cblkA, cblkB );
        }
    }

    kernel_trace_stop( cblkB->fblokptr->inlast, ktype, m, n, 0, flops, time );
}

