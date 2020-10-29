/**
 *
 * @file coeftab.h
 *
 * PaStiX coefficient array routines header.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2020-01-06
 *
 * @addtogroup coeftab
 * @{
 *   This group collects all the functions that operate on the full matrix and
 *   which are not factorization/solve routines.
 *
 **/
#ifndef _coeftab_h_
#define _coeftab_h_

#include "sopalin/coeftab_z.h"
#include "sopalin/coeftab_c.h"
#include "sopalin/coeftab_d.h"
#include "sopalin/coeftab_s.h"

void coeftabInit( pastix_data_t     *pastix_data,
                  pastix_coefside_t  side );
void coeftabExit( SolverMatrix      *solvmtx );

pastix_int_t coeftabCompress( pastix_data_t *pastix_data );

#if defined(PASTIX_WITH_MPI)
void coeftab_scatter( SolverMatrix *solvmtx,
                      PASTIX_Comm   comm,
                      pastix_int_t  root,
                      pastix_coeftype_t typesze );
void coeftab_gather ( SolverMatrix *solvmtx,
                      PASTIX_Comm   comm,
                      pastix_int_t  root,
                      pastix_coeftype_t typesze );
void coeftab_nullify( SolverMatrix *solvmtx );
#endif

/**
 * @brief Type of the memory gain functions
 */
typedef void (*coeftab_fct_memory_t)( SolverMatrix * );

/**
 * @brief List of functions to compute the memory gain in low-rank per precision.
 */
coeftab_fct_memory_t coeftabMemory[4];

/**
 * @}
 */
#endif /* _coeftab_h_ */
