/**
 *
 * @file solver_io.c
 *
 * PaStiX solver structure I/O routines.
 *
 * @copyright 2004-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2019-12-09
 *
 **/
#include "common.h"
#include "symbol.h"
#include "queue.h"
#include "solver.h"

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Load a solver matrix structure from a file.
 *
 *******************************************************************************
 *
 * @param[inout] solvptr
 *          The allocated pointer to a solver structure initialized. No need to
 *          initialized it with solverInit().
 *
 * @param[in] stream
 *          The stream where to read the informations.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_FILE if a problem occurs during the read.
 *
 *******************************************************************************/
int
solverLoad( SolverMatrix *solvptr,
            FILE         *stream  )
{
    pastix_int_t   i,j;
    pastix_int_t   clustnbr, clustnum;
    SolverCblk    *cblkptr;
    SolverCblk    *cblktnd;
    SolverBlok    *blokptr;
    SolverBlok    *bloktnd;
    Task          *taskptr;
    Task          *tasknd;

    pastix_int_t   versval;
    pastix_int_t   baseval;
    pastix_int_t   nodenbr;
    pastix_int_t   cblknbr;
    pastix_int_t   cblknum;
    pastix_int_t   bloknbr;
    pastix_int_t   bloknum;

    solverInit(solvptr);

    /** Load the symbol matrix **/
    if ((intLoad (stream, &versval) +               /* Read header */
         intLoad (stream, &cblknbr) +
         intLoad (stream, &bloknbr) +
         intLoad (stream, &nodenbr) +
         intLoad (stream, &baseval) != 5) ||
        (versval < 0)                     ||        /* Version should be 0 or 1 */
        (versval > 2)                     ||
        (bloknbr < cblknbr)               ||
        (nodenbr < cblknbr)) {
        errorPrint ("solverLoad: bad input (1)");
        return PASTIX_ERR_FILE;
    }

    if ( versval > 1 ) {
        errorPrint("solverLoad: Version 0 and 1 of the solver files are not supported anymore");
        return PASTIX_ERR_FILE;
    }

    MALLOC_INTERN(solvptr->cblktab, cblknbr + 1, SolverCblk);
    MALLOC_INTERN(solvptr->bloktab, bloknbr,     SolverBlok);
    if (solvptr->cblktab == NULL || solvptr->bloktab == NULL) {
        errorPrint ("solverLoad: out of memory");
        solverExit (solvptr);
        solverInit (solvptr);
        return PASTIX_ERR_FILE;
    }
    solvptr->baseval = baseval;
    solvptr->cblknbr = cblknbr;
    solvptr->bloknbr = bloknbr;
    solvptr->nodenbr = nodenbr;

    for (cblknum = 0; cblknum < cblknbr; cblknum ++) {
        if ((intLoad (stream, &solvptr->cblktab[cblknum].fcolnum) + /* Read column blocks */
             intLoad (stream, &solvptr->cblktab[cblknum].lcolnum) +
             intLoad (stream, &bloknum) != 3) ||
            (solvptr->cblktab[cblknum].fcolnum > solvptr->cblktab[cblknum].lcolnum)) {
            errorPrint ("solverLoad: bad input (2)");
            /* solverExit (solvptr); */
            /* solverInit (solvptr); */
            return PASTIX_ERR_FILE;
        }
        solvptr->cblktab[cblknum].fblokptr = solvptr->bloktab+bloknum;
    }
    solvptr->cblktab[cblknbr].fcolnum =             /* Set last column block */
        solvptr->cblktab[cblknbr].lcolnum = nodenbr + baseval;
    solvptr->cblktab[cblknbr].fblokptr = solvptr->bloktab + bloknbr;

    for (bloknum = 0; bloknum < bloknbr; bloknum ++) {
        if ((intLoad (stream, &solvptr->bloktab[bloknum].frownum) + /* Read column blocks */
             intLoad (stream, &solvptr->bloktab[bloknum].lrownum) +
             intLoad (stream, &solvptr->bloktab[bloknum].fcblknm) != 3) ||
            (solvptr->bloktab[bloknum].frownum > solvptr->bloktab[bloknum].lrownum)) {
            errorPrint ("solverLoad: bad input (3)");
            solverExit (solvptr);
            solverInit (solvptr);
            return PASTIX_ERR_FILE;
        }

        {
            pastix_int_t levfval;
            if ((versval == 0) &&
                ((intLoad (stream, &levfval) != 1) ||
                 (levfval < 0))) {
                errorPrint ("solverLoad: bad input (4)");
                solverExit (solvptr);
                solverInit (solvptr);
                return PASTIX_ERR_FILE;
            }
        }
    }


    if( intLoad (stream, &solvptr->coefnbr) +
        intLoad (stream, &solvptr->gemmmax) +
        intLoad (stream, &solvptr->nbftmax) +
        intLoad (stream, &solvptr->arftmax) +
        intLoad (stream, &clustnum)         +
        intLoad (stream, &clustnbr)         +
        intLoad (stream, &solvptr->tasknbr) +
        intLoad (stream, &solvptr->procnbr) +
        intLoad (stream, &solvptr->thrdnbr)
        != 11 )
    {
        errorPrint ("solverLoad: bad input (1)");
        return PASTIX_ERR_FILE;
    }

    solvptr->clustnbr = (pastix_int_t)clustnbr;
    solvptr->clustnum = (pastix_int_t)clustnum;

    if (((solvptr->cblktab = (SolverCblk *)   memAlloc((solvptr->cblknbr + 1) * sizeof(SolverCblk)    )) == NULL) ||
        ((solvptr->bloktab = (SolverBlok *)   memAlloc( solvptr->bloknbr      * sizeof(SolverBlok)    )) == NULL) ||
        ((solvptr->tasktab = (Task *)         memAlloc((solvptr->tasknbr+1)   * sizeof(Task)          )) == NULL) ||
        ((solvptr->ttsknbr = (pastix_int_t *) memAlloc((solvptr->thrdnbr)     * sizeof(pastix_int_t)  )) == NULL) ||
        ((solvptr->ttsktab = (pastix_int_t **)memAlloc((solvptr->thrdnbr)     * sizeof(pastix_int_t *))) == NULL) )
    {
        errorPrint ("solverLoad: out of memory (1)");
        if (solvptr->cblktab != NULL) {
            memFree_null (solvptr->cblktab);
        }
        if (solvptr->bloktab != NULL) {
            memFree_null (solvptr->bloktab);
        }
        if (solvptr->tasktab != NULL) {
            memFree_null (solvptr->tasktab);
        }
        return PASTIX_ERR_FILE;
    }

    for (cblkptr = solvptr->cblktab,                /* Read column block data */
             cblktnd = cblkptr + solvptr->cblknbr;
         cblkptr < cblktnd; cblkptr ++)
    {
        if (intLoad (stream, &cblkptr->stride) != 1)
        {
            errorPrint ("solverlLoad: bad input (2)");
            solverExit (solvptr);
            return     PASTIX_ERR_FILE;
        }
    }

    for (blokptr = solvptr->bloktab,                /* Read block data */
             bloktnd = blokptr + solvptr->bloknbr;
         blokptr < bloktnd; blokptr ++)
    {
        if (intLoad (stream, &blokptr->coefind) != 1)
        {
            errorPrint ("solverLoad: bad input (3)");
            solverExit (solvptr);
            return     PASTIX_ERR_FILE;
        }

    }

    for (taskptr = solvptr->tasktab,                /** Read Task data **/
             tasknd = taskptr + solvptr->tasknbr +1;
         (taskptr < tasknd); taskptr ++)
    {
        pastix_int_t temp;

        intLoad(stream, &(taskptr->taskid));
        intLoad(stream, &(taskptr->prionum));
        intLoad(stream, &(taskptr->cblknum));
        intLoad(stream, &(taskptr->bloknum));
        {
            /* volatile pb alpha */
            intLoad(stream, &temp);
            taskptr->ctrbcnt = temp;
        }
    }

    for(i=0;i<solvptr->thrdnbr;i++)                 /** Read task by thread data **/
    {
        intLoad(stream, &(solvptr->ttsknbr[i]));
        MALLOC_INTERN(solvptr->ttsktab[i], solvptr->ttsknbr[i], pastix_int_t);
        if (solvptr->ttsktab[i] == NULL)
        {
            errorPrint ("solverLoad: out of memory (1)");
            return 1;
        }
        for (j=0;j<solvptr->ttsknbr[i];j++)
        {
            intLoad(stream, &(solvptr->ttsktab[i][j]));
        }
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Save a solver matrix structure into a file.
 *
 *******************************************************************************
 *
 * @param[inout] solvptr
 *          The solver matrix structure to dump to disk.
 *
 * @param[in] stream
 *          The stream where to write the ordering.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if the ordeptr structure is incorrect,
 * @retval PASTIX_ERR_FILE if a problem occurs during the write.
 *
 *******************************************************************************/
int
solverSave( const SolverMatrix *solvptr,
            FILE               *stream  )
{
    pastix_int_t   i, j, o;
    SolverCblk    *cblkptr;
    SolverCblk    *cblktnd;
    SolverBlok    *blokptr;
    SolverBlok    *bloktnd;
    Task          *taskptr;

    /* Save the solver matrix */
    {
        const SolverCblk *cblktnd;
        const SolverCblk *cblkptr;
        const SolverBlok *bloktnd;
        const SolverBlok *blokptr;

        o = (fprintf (stream, "2\n%ld\t%ld\t%ld\t%ld\n", /* Write file header */
                      (long) solvptr->cblknbr,
                      (long) solvptr->bloknbr,
                      (long) solvptr->nodenbr,
                      (long) solvptr->baseval) == EOF);
        for (cblkptr = solvptr->cblktab, cblktnd = cblkptr + solvptr->cblknbr;
             (cblkptr < cblktnd) && (o == 0); cblkptr ++) {
            o = (fprintf (stream, "%ld\t%ld\t%ld\n",
                          (long) cblkptr->fcolnum,
                          (long) cblkptr->lcolnum,
                          (long) (cblkptr->fblokptr - solvptr->bloktab)) == EOF);
        }
        for (blokptr = solvptr->bloktab, bloktnd = blokptr + solvptr->bloknbr;
             (blokptr < bloktnd) && (o == 0); blokptr ++) {
            o = (fprintf (stream, "%ld\t%ld\t%ld\n",/* "%ld\t%ld\t%ld\t%ld\n", */
                          (long) blokptr->frownum,
                          (long) blokptr->lrownum,
                          (long) blokptr->fcblknm/* , */
                          /* (long) blokptr->levfval */) == EOF);
        }
    }

    /* Write file header */
    o = (fprintf (stream, "\n%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",
                  (long) solvptr->coefnbr,
                  (long) solvptr->gemmmax,
                  (long) solvptr->nbftmax,
                  (long) solvptr->arftmax,
                  (long) solvptr->clustnum,
                  (long) solvptr->clustnbr,
                  (long) solvptr->tasknbr,
                  (long) solvptr->procnbr,
                  (long) solvptr->thrdnbr
                  ) == EOF);

    /* write cblk data */
    for (cblkptr = solvptr->cblktab, cblktnd = cblkptr + solvptr->cblknbr;
         (cblkptr < cblktnd) && (o == 0); cblkptr ++)
    {
        o = (fprintf (stream, "%ld\n",
                      (long) cblkptr->stride) == EOF);
    }

    /* write blok data */
    for (blokptr = solvptr->bloktab,
             bloktnd = blokptr + solvptr->bloknbr;
         (blokptr < bloktnd) && (o == 0); blokptr ++)
    {
        o = (fprintf (stream, "%ld\n",(long) blokptr->coefind) == EOF);
    }

    fprintf(stream, "\n");
    fprintf(stream, "\n");

    /* Write Task data */
    {
        Task *taskend = solvptr->tasktab + solvptr->tasknbr;
        for (taskptr = solvptr->tasktab;
             (taskptr < taskend) && (o==0); taskptr ++)
        {
            fprintf(stream, "%ld\t%ld\t%ld\t%ld\t%ld\n",
                    (long)taskptr->taskid, (long)taskptr->prionum, (long)taskptr->cblknum, (long)taskptr->bloknum,
                    (long)taskptr->ctrbcnt);
        }
    }

    /* Write ttsktab */
    for (i=0; i<solvptr->thrdnbr; i++)
    {
        fprintf(stream, "%ld\n", (long)solvptr->ttsknbr[i]);
        for (j=0; j<solvptr->ttsknbr[i]; j++)
        {
            fprintf(stream, "%ld\n", (long)solvptr->ttsktab[i][j]);
        }
    }

    return o ? PASTIX_ERR_FILE : PASTIX_SUCCESS;
}
