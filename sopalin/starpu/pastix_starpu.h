/**
 *
 * @file pastix_starpu.h
 *
 * StarPU support for the numerical factorization and solve of PaStiX.
 *
 * @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2020-01-06
 *
 * @addtogroup pastix_starpu
 * @{
 *   This module describes the functionnality provided by the runtime system
 *   StarPU for the numerical factorization and solve.
 *
 **/
#ifndef _pastix_starpu_h_
#define _pastix_starpu_h_

#include "common.h"
#include "solver.h"

#if defined(PASTIX_WITH_MPI)
#include <starpu_mpi.h>
#else
#include <starpu.h>
#endif

#include <starpu_profiling.h>

#if defined(PASTIX_WITH_CUDA) && !defined(PASTIX_STARPU_SIMULATION)
#include <starpu_scheduler.h>
#include <starpu_cuda.h>

#include <cublas.h>
#include <starpu_cublas.h>
#if defined(PASTIX_WITH_CUBLAS_V2)
#include <cublas_v2.h>
#include <starpu_cublas_v2.h>
#endif
#endif

typedef struct starpu_conf starpu_conf_t;

#if defined(PASTIX_STARPU_SYNC)
#define TASK_SYNCHRONOUS , STARPU_TASK_SYNCHRONOUS, 1
#else
#define TASK_SYNCHRONOUS
#endif

#if defined(PASTIX_WITH_MPI)
#define starpu_insert_task starpu_mpi_insert_task
#define pastix_codelet(_codelet_) sopalin_data->solvmtx->solv_comm, _codelet_ TASK_SYNCHRONOUS
#else
#define pastix_codelet(_codelet_) _codelet_ TASK_SYNCHRONOUS
#endif

/**
 * @brief Additional StarPU handlers for a column-block when using 2D kernels.
 *
 * Handle requirements for contiguous allocation of the block handlers when
 * using StarPU data partitioning.
 */
typedef struct starpu_cblk_s {
    pastix_int_t          handlenbr; /**< Number of 2D block handlers in the column-block */
    starpu_data_handle_t *handletab; /**< Array of 2D block handlers for the column-block */
} starpu_cblk_t;

/**
 * @brief StarPU descriptor stucture for the sparse matrix.
 */
typedef struct starpu_sparse_matrix_desc_s {
    int64_t         mpitag;         /**< MPI id of StarPU */
    int             typesze;        /**< Arithmetic size                                                                              */
    int             mtxtype;        /**< Matrix structure: PastixGeneral, PastixSymmetric or PastixHermitian.                         */
    SolverMatrix   *solvmtx;        /**< Solver matrix structure that describes the problem and stores the original data              */
    starpu_cblk_t  *cblktab_handle; /**< Array of 2D column-block handlers (NULL when using 1D kernels only)                          */
    void          **gpu_blocktab;     /**< Pointer to GPU arrays that contains frownum,lrownum of each block for Fermi (NULL otherwise) */
} starpu_sparse_matrix_desc_t;

/**
 * @brief StarPU descriptor for the vectors linked to a given sparse matrix.
 */
typedef struct starpu_dense_matrix_desc_s {
    int64_t               mpitag;    /**< MPI id of StarPU */
    int                   ncol;      /**< Number of columns of the matrix                                                 */
    int                   typesze;   /**< Arithmetic size                                                                 */
    SolverMatrix         *solvmtx;   /**< Solver matrix structure that describes the problem and stores the original data */
    starpu_data_handle_t *handletab; /**< Array of handlers for the blocks */
    void                 *dataptr;   /**< Store the main data pointer to check that the descriptor matches the reference  */
} starpu_dense_matrix_desc_t;

void starpu_sparse_matrix_init    ( SolverMatrix *solvmtx,
                                    int typesize, int mtxtype,
                                    int nodes, int myrank );
void starpu_sparse_matrix_destroy ( starpu_sparse_matrix_desc_t *desc );
void starpu_sparse_matrix_getoncpu( starpu_sparse_matrix_desc_t *desc );

void starpu_dense_matrix_init     ( SolverMatrix *solvmtx,
                                    pastix_int_t ncol, char *A, pastix_int_t lda,
                                    int typesze, int nodes, int myrank );
void starpu_dense_matrix_destroy  ( starpu_dense_matrix_desc_t *desc );
void starpu_dense_matrix_getoncpu ( starpu_dense_matrix_desc_t *desc );

void pastix_starpu_init( pastix_data_t *pastix,
                         int *argc, char **argv[],
                         const int *bindtab );
void pastix_starpu_finalize( pastix_data_t *pastix );

int64_t pastix_starpu_get_tag( );

void
pastix_starpu_partition_submit( pastix_coefside_t side,
                                SolverCblk       *cblk,
                                starpu_cblk_t    *starpu_cblk );
void
pastix_starpu_unpartition_submit( const starpu_sparse_matrix_desc_t *spmtx,
                                  int rank, pastix_coefside_t side,
                                  SolverCblk    *cblk,
                                  starpu_cblk_t *starpu_cblk );

#endif /* _pastix_starpu_h_ */

/**
 *@}
 */
