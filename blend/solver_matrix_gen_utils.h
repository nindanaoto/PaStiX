/**
 *
 * @file solver_matrix_gen_utils.h
 *
 * PaStiX solver structure generation functions to factorize
 * solver_matric_gen.c .
 *
 * @copyright 1998-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.1.0
 * @author Tony Delarue
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @date 2020-01-29
 *
 * @addtogroup blend_dev_solver
 * @{
 *
 **/
#ifndef _solver_matrix_gen_utils_h_
#define _solver_matrix_gen_utils_h_

/**
 *******************************************************************************
 *
 * @brief Get the expanded column indexes of a symbol_cblk.
 *
 *******************************************************************************
 *
 * @param[in] symbmtx
 *          Pointer to the symbol matrix.
 *
 * @param[in] symbcblk
 *          The pointer to the current symbol_cblk.
 *
 * @param[inout] fcolnum
 *          First column index of the current cblk.
 *
 * @param[inout] symbcblk
 *          Last column index of the current cblk.
 *
 * @return The number of columns of the expanded cblk.
 *
 *******************************************************************************/
static inline pastix_int_t
solvMatGen_get_colnum( const symbol_matrix_t *symbmtx,
                       symbol_cblk_t         *symbcblk,
                       pastix_int_t          *fcolnum,
                       pastix_int_t          *lcolnum )
{
    if ( symbmtx->dof < 0 ) {
        *fcolnum = symbmtx->dofs[symbcblk->fcolnum];
        *lcolnum = symbmtx->dofs[symbcblk->lcolnum + 1] - 1;
    }
    else {
        *fcolnum = symbmtx->dof *   symbcblk->fcolnum;
        *lcolnum = symbmtx->dof * ( symbcblk->lcolnum + 1 ) - 1;
    }
    return (*lcolnum) - (*fcolnum) + 1;
}

/**
 *******************************************************************************
 *
 * @brief Get the expanded row index of a symbol_blok.
 *
 *******************************************************************************
 *
 * @param[in] symbmtx
 *          Pointer to the symbol matrix.
 *
 * @param[in] symbblok
 *          The pointer to the current symbol_blok.
 *
 * @param[inout] frownum
 *          First row index of the current blok.
 *
 * @param[inout] lrownum
 *          Last row index of the current blok.
 *
 * @return The number of rows of the expanded blok.
 *
 *******************************************************************************/
static inline pastix_int_t
solvMatGen_get_rownum( const symbol_matrix_t *symbmtx,
                       symbol_blok_t         *symbblok,
                       pastix_int_t          *frownum,
                       pastix_int_t          *lrownum )
{
    if ( symbmtx->dof < 0 ) {
        *frownum = symbmtx->dofs[symbblok->frownum];
        *lrownum = symbmtx->dofs[symbblok->lrownum + 1] - 1;
    }
    else {
        *frownum = symbmtx->dof *   symbblok->frownum;
        *lrownum = symbmtx->dof * ( symbblok->lrownum + 1 ) - 1;
    }
    return (*lrownum) - (*frownum) + 1;
}

/**
 *******************************************************************************
 *
 * @brief Init a local solver block, virtual or not.
 *
 *******************************************************************************
 *
 * @param[inout] solvblok
 *          The pointer to the block.
 *
 * @param[in] fcblknm
 *          First column block index.
 *
 *  @param[in] lcblknm
 *          Last column block index.
 *
 *  @param[in] frownum
 *          First row of the block.
 *
 * @param[in] lrownum
 *          Last row of the block.
 *
 * @param[in] stride
 *          Stride of the block.
 *
 * @param[in] nbcols
 *          Amount of columns in the cblk. Not set in the cblk at this time.
 *
 * @param[in] layout2D
 *          Parameter which indicates if it is a 2D block
 *
 ********************************************************************************/
static inline void
solvMatGen_init_blok( SolverBlok  *solvblok,
                      pastix_int_t lcblknm,
                      pastix_int_t fcblknm,
                      pastix_int_t frownum,
                      pastix_int_t lrownum,
                      pastix_int_t stride,
                      pastix_int_t nbcols,
                      pastix_int_t layout2D )
{
    assert( fcblknm >= -1 );
    assert( lcblknm >= 0 );
    assert( (fcblknm == -1) || (lcblknm <= fcblknm) );
    assert( frownum >= 0 );
    assert( lrownum >= frownum );
    assert( stride  >= 0 );
    assert( nbcols  >= 0 );

    solvblok->handler[0] = NULL;
    solvblok->handler[1] = NULL;
    solvblok->fcblknm    = fcblknm;
    solvblok->lcblknm    = lcblknm;
    solvblok->frownum    = frownum;
    solvblok->lrownum    = lrownum;
    solvblok->coefind    = layout2D ? stride * nbcols : stride;
    solvblok->browind    = -1;
    solvblok->gpuid      = GPUID_UNDEFINED;
    solvblok->inlast     = 0;
    solvblok->LRblock    = NULL;
}

/**
 *******************************************************************************
 *
 * @brief Init a local solver cblk, virtual or not.
 *
 *******************************************************************************
 *
 * @param[inout] solvcblk
 *          The pointer to the cblk.
 *
 * @param[inout] fblokptr
 *          The pointer to the first block.
 *
 * @param[in] candcblk
 *          Allow us to know the type of the cblk.
 *
 * @param[in] symbcblk
 *          Corresponding global cblk in the symbol matrix.
 *
 * @param[in] fcolnum
 *          First column of the cblk.
 *
 * @param[in] lcolnum
 *          Last column of the cblk.
 *
 * @param[in] brownum
 *          Index in th browtab.
 *
 * @param[in] stride
 *          Stride of the cblk.
 *
 * @param[in] nodenbr
 *          Current global column index.
 *
 * @param[in] cblknum
 *          Global cblkindex. -1 if virtual.
 *
 *******************************************************************************/
static inline void
solvMatGen_init_cblk( SolverCblk    *solvcblk,
                      SolverBlok    *fblokptr,
                      Cand          *candcblk,
                      symbol_cblk_t *symbcblk,
                      pastix_int_t   fcolnum,
                      pastix_int_t   lcolnum,
                      pastix_int_t   brownum,
                      pastix_int_t   stride,
                      pastix_int_t   nodenbr,
                      pastix_int_t   cblknum,
                      int            ownerid )
{
    assert( fblokptr != NULL );
    assert( fcolnum >= 0 );
    assert( lcolnum >= fcolnum );
    assert( stride  >= 0 );
    assert( nodenbr >= 0 );
    assert( brownum >= 0 );

    /* Init the cblk */
    solvcblk->lock       = PASTIX_ATOMIC_UNLOCKED;
    solvcblk->ctrbcnt    = -1;
    solvcblk->cblktype   = (cblknum == -1) ? 0 : candcblk->cblktype;
    solvcblk->gpuid      = GPUID_UNDEFINED;
    solvcblk->fcolnum    = fcolnum;
    solvcblk->lcolnum    = lcolnum;
    solvcblk->fblokptr   = fblokptr;
    solvcblk->stride     = stride;
    solvcblk->lcolidx    = nodenbr;
    solvcblk->brownum    = brownum;
    solvcblk->gcblknum   = cblknum;
    solvcblk->bcscnum    = -1;
    solvcblk->selevtx    = (symbcblk->selevtx == SYMBCBLK_PROJ) ? 1 : 0;
    solvcblk->ownerid    = ownerid;
    solvcblk->lcoeftab   = NULL;
    solvcblk->ucoeftab   = NULL;
    solvcblk->handler[0] = NULL;
    solvcblk->handler[1] = NULL;
    solvcblk->threadid   = -1;
}

/**
 *******************************************************************************
 *
 * @brief Make the link between the supernodes and the cblks.
 *
 *******************************************************************************
 *
 * @param[inout] solvcblk
 *          The pointer to the cblk.
 *
 * @param[in] sndeidx
 *          The index of the last visited supernode.
 *
 * @param[in] ordeptr
 *          The ordering structure.
 *
 * @return  The current sndeidx that can be used to reduce the cost of further calls.
 *
 *******************************************************************************/
static inline pastix_int_t
solvMatGen_supernode_index( SolverCblk           *solvcblk,
                            pastix_int_t          sndeidx,
                            const pastix_order_t *ordeptr )
{
    while ( (sndeidx < ordeptr->sndenbr ) &&
            (ordeptr->sndetab[sndeidx+1] <= solvcblk->lcolnum) )
    {
        sndeidx++;
    }
    assert( (ordeptr->sndetab[sndeidx]   <= solvcblk->fcolnum) &&
            (ordeptr->sndetab[sndeidx+1] >  solvcblk->lcolnum) );
    solvcblk->sndeidx = sndeidx;

    /* Register the cblk as being part of the last supernode */
    if ( solvcblk->sndeidx+1 == ordeptr->sndenbr ) {
        solvcblk->cblktype = solvcblk->cblktype | CBLK_IN_LAST;
    }

    return sndeidx;
}

/**
 *******************************************************************************
 *
 * @brief Update the 1D/2D infos of the solvmtx through a cblk.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          Pointer to the solver matrix.
 *
 * @param[inout] nbcblk2d
 *          Amount of 2D cblk.
 *
 * @param[inout] nbblok2d
 *          Amount of 2D blocks.
 *
 * @param[in] nbbloks
 *          Amount blocks in the current cblk.
 *
 * @param[in] tasks2D
 *          Boolean which indicate if the task is 2D.
 *
 * @param[in] cblknum
 *          Current cblk index.
 *
 *******************************************************************************/
static inline void
solvMatGen_cblkIs2D( SolverMatrix *solvmtx,
                     pastix_int_t *nbcblk2d,
                     pastix_int_t *nbblok2d,
                     pastix_int_t  nbbloks,
                     pastix_int_t  tasks2D,
                     pastix_int_t  cblknum )
{
    /*
     * 2D tasks: Compute the number of cblk split in 2D tasks, and
     * the smallest id
     */
    if ( tasks2D ) {
        solvmtx->cblkmin2d = pastix_imin( solvmtx->cblkmin2d, cblknum );
        *nbcblk2d += 1;
        *nbblok2d += nbbloks;
    }
    else {
        solvmtx->cblkmax1d = pastix_imax( solvmtx->cblkmax1d, cblknum );
    }

    /*
     * Compute the maximum number of block per cblk for data
     * structure in PaRSEC/StarPU
     */
    if ( cblknum >= solvmtx->cblkmin2d ) {
        solvmtx->cblkmaxblk = pastix_imax( solvmtx->cblkmaxblk, nbbloks );
    }
}

void solvMatGen_fill_localnums( const symbol_matrix_t *symbmtx,
                                const SimuCtrl  *simuctrl,
                                SolverMatrix    *solvmtx,
                                pastix_int_t    *cblklocalnum,
                                pastix_int_t    *bloklocalnum,
                                pastix_int_t    *tasklocalnum,
                                pastix_int_t    *fcbklocalnum,
                                pastix_int_t    *pcbklocalnum,
                                int            **recv_sources_ptr );

void solvMatGen_init_cblk_recv( const symbol_matrix_t *symbmtx,
                                      SolverCblk      *solvcblk,
                                      SolverBlok      *solvblok,
                                      Cand            *candcblk,
                                      pastix_int_t    *cblklocalnum,
                                      pastix_int_t     recvidx,
                                      pastix_int_t     fcolnum,
                                      pastix_int_t     lcolnum,
                                      pastix_int_t     brownum,
                                      pastix_int_t     nodenbr,
                                      pastix_int_t     cblknum,
                                      int              ownerid );

pastix_int_t
solvMatGen_reorder_browtab( const symbol_matrix_t *symbmtx,
                            symbol_cblk_t         *symbcblk,
                            SolverMatrix          *solvmtx,
                            SolverCblk            *solvcblk,
                            pastix_int_t          *browtmp,
                            pastix_int_t          *cblklocalnum,
                            pastix_int_t          *bloklocalnum,
                            pastix_int_t           brownum );

void solvMatGen_fill_tasktab( SolverMatrix   *solvmtx,
                              isched_t       *isched,
                              const SimuCtrl *simuctrl,
                              pastix_int_t   *tasklocalnum,
                              pastix_int_t   *cblklocalnum,
                              pastix_int_t   *bloklocalnum,
                              pastix_int_t    clustnum,
                              int             is_dbg );

void solvMatGen_stats_last( SolverMatrix *solvmtx );
void solvMatGen_max_buffers( SolverMatrix *solvmtx );

#endif /* _solver_matrix_gen_utils_h_ */

/**
 *@}
 */
