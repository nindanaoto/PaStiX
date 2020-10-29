/**
 *
 * @file bvec_zcompute.c
 *
 *  Functions computing operations on the BCSC.
 *
 * @copyright 2004-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Théophile terraz
 * @date 2020-01-26
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <math.h>
#include "lapacke.h"
#include "bcsc.h"
#include "bcsc_z.h"
#include "frobeniusupdate.h"
#include "cblas.h"
#include "solver.h"

#if defined(PASTIX_WITH_MPI)
void
bvec_zmpi_frb_merge( double       *dist,
                     double       *loc,
                     int          *len,
                     MPI_Datatype *dtype )
{
    assert( *len == 2 );
    frobenius_merge( dist[0], dist[1], loc, loc+1 );
    (void)len;
    (void)dtype;
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Compute the norm 2 of a vector. (Sequential version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          Provide information to the parallel version, to know the global
 *          context (Number of thread, barrier, ...).
 *
 * @param[in] n
 *          The size of the vector x.
 *
 * @param[in] x
 *          The vector x of size n.
 *
 *******************************************************************************
 *
 * @retval the norm 2 of x.
 *
 *******************************************************************************/
double
bvec_znrm2_seq( pastix_data_t            *pastix_data,
                pastix_int_t              n,
                const pastix_complex64_t *x )
{
    SolverMatrix  *solvmtx = pastix_data->solvmatr;
    SolverCblk    *scblk;
    pastix_bcsc_t *bcsc    = pastix_data->bcsc;
    bcsc_cblk_t   *bcblk   = bcsc->cscftab;
    pastix_int_t   cblknbr, colnbr;
    double         data[] = { 0., 1. }; /* Scale, Sum */
    double         norm;
    const double  *valptr;
    pastix_int_t   i, j;

    cblknbr = bcsc->cscfnbr;
    for( i = 0; i < cblknbr; i++, bcblk++ ) {
        scblk  = solvmtx->cblktab + bcblk->cblknum;
        colnbr = cblk_colnbr( scblk );
        valptr = (const double*)(x + scblk->lcolidx);

        for( j=0; j < colnbr; j++, valptr++ )
        {
            /* Real part */
            frobenius_update( 1, data, data + 1, valptr );
#if defined(PRECISION_z) || defined(PRECISION_c)
            /* Imaginary part */
            valptr++;
            frobenius_update( 1, data, data + 1, valptr );
#endif
        }
    }

#if defined(PASTIX_WITH_MPI)
    {
        MPI_Op merge;

        MPI_Op_create( (MPI_User_function *)bvec_zmpi_frb_merge, 1, &merge );
        MPI_Allreduce( MPI_IN_PLACE, data, 2, MPI_DOUBLE, merge, solvmtx->solv_comm );
        MPI_Op_free( &merge );
    }
#endif

    norm = data[0] * sqrt( data[1] );

    (void)n;
    return norm;
}

struct z_argument_nrm2_s
{
    pastix_int_t              n;
    const pastix_complex64_t *x;
    pastix_atomic_lock_t      lock;
    double                    scale;
    double                    sumsq;
};

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Compute the norm 2 of a vector. (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          The context of the current thread
 *
 * @param[inout] args
 *          The parameter that providenumber of elements of x,
 *          and the vector which norm2 is to be computed, and the norm value
 *
 *******************************************************************************/
static inline void
pthread_bvec_znrm2( isched_thread_t *ctx,
                   void            *args )
{
    struct z_argument_nrm2_s *arg = (struct z_argument_nrm2_s*)args;
    pastix_int_t              n = arg->n;
    const pastix_complex64_t *x = arg->x;
    double                    scale = 0.;
    double                    sumsq = 1.;
    double                   *valptr = (double*)x;
    pastix_int_t              i, rank, size;
    pastix_int_t              begin, end;

    size = ctx->global_ctx->world_size;
    rank = ctx->rank;

    begin = (n / size) * rank;
    if (rank == (size - 1)) {
        end = n; /* One iteration more */
    } else {
        end = (n / size) * (rank + 1);
    }

    valptr += begin;
#if defined(PRECISION_z) || defined(PRECISION_c)
    valptr += begin;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    for( i = begin; i < end; i++, valptr++ )
    {
        frobenius_update( 1, &scale, &sumsq, valptr );
#if defined(PRECISION_z) || defined(PRECISION_c)
        valptr ++;
        frobenius_update( 1, &scale, &sumsq, valptr );
#endif
    }

    /* If we computed something */
    if ( scale != 0. ) {
        pastix_atomic_lock( &(arg->lock) );
        frobenius_merge( scale, sumsq, &(arg->scale), &(arg->sumsq) );
        pastix_atomic_unlock( &(arg->lock) );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Compute the norm 2 of a vector. (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          Provide information to the parallel version, to know the global
 *          context (Number of thread, barrier, ...).
 *
 * @param[in] x
 *          The vector which norm2 is to be computed
 *
 * @param[in] n
 *          The number of elements of x
 *
 *******************************************************************************
 *
 * @return The norm 2 of the vector
 *
 *******************************************************************************/
double
bvec_znrm2_smp( pastix_data_t            *pastix_data,
                pastix_int_t              n,
                const pastix_complex64_t *x )
{
    struct z_argument_nrm2_s arg = { n, x, PASTIX_ATOMIC_UNLOCKED, 0., 1. };
    isched_parallel_call( pastix_data->isched, pthread_bvec_znrm2, &arg );

    return arg.scale * sqrt( arg.sumsq );
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Scale a vector by the scalar alpha. (Sequential version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          Provide information to the parallel version, to know the global
 *          context (Number of thread, barrier, ...).
 *
 * @param[in] n
 *          The size of the vector x.
 *
 * @param[in] alpha
 *          The scalar to scale the vector x.
 *
 * @param[inout] x
 *          The vector x to scale.
 *
 *******************************************************************************/
void
bvec_zscal_seq( pastix_data_t      *pastix_data,
                pastix_int_t        n,
                pastix_complex64_t  alpha,
                pastix_complex64_t *x )
{
#if defined(PASTIX_WITH_MPI) && 0
    SolverMatrix  *solvmtx = pastix_data->solvmatr;
    SolverCblk    *scblk   = solvmtx->cblktab;
    pastix_bcsc_t *bcsc    = pastix_data->bcsc;
    bcsc_cblk_t   *bcblk   = bcsc->cscftab;
    pastix_int_t   i, cblknbr;

    cblknbr = bcsc->cscfnbr;
    for( i = 0; i < cblknbr; i++, bcblk++ ) {
        scblk  = solvmtx->cblktab + bcblk->cblknum;
        n = cblk_colnbr( scblk );

        cblas_zscal( n, CBLAS_SADDR(alpha), x + scblk->lcolidx, 1 );
    }
#else
    (void)pastix_data;
    cblas_zscal( n, CBLAS_SADDR(alpha), x, 1 );
#endif

}

struct z_argument_scal_s
{
  pastix_int_t        n;
  pastix_complex64_t  alpha;
  pastix_complex64_t *x;
};

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Scale a vector (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          The information about number of thread and rank of the actual thread
 *
 * @param[inout] args
 *          The argument providing number of elements of the vector, the scaling
 *          parameter and, and the vector to be scaled
 *
 *******************************************************************************/
static inline void
pthread_bvec_zscal( isched_thread_t *ctx,
                    void            *args )
{
    struct z_argument_scal_s *arg   = (struct z_argument_scal_s*)args;
    pastix_complex64_t       *x     = arg->x;
    pastix_int_t              n     = arg->n;
    pastix_complex64_t        alpha = arg->alpha;
    pastix_int_t              size  = ctx->global_ctx->world_size;
    pastix_int_t              rank;
    pastix_int_t              begin, end;

    if( x == NULL ) {
        return;
    }

    rank = (pastix_int_t)ctx->rank;
    begin = (n/size) * rank;
    if (rank == (size - 1)) {
        end = n;
    } else {
        end = (n/size) * (rank + 1);
    }

    if ( (end - begin) > 0 ) {
        cblas_zscal( end - begin, CBLAS_SADDR(alpha), x + begin, 1 );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Scale a vector (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
 *
 * @param[in] n
 *          The number of elements of the vector
 *
 * @param[in] alpha
 *          The scaling parameter
 *
 * @param[inout] x
 *          The vector to be scaled
 *
 *******************************************************************************/
void
bvec_zscal_smp( pastix_data_t      *pastix_data,
                pastix_int_t        n,
                pastix_complex64_t  alpha,
                pastix_complex64_t *x )
{
    struct z_argument_scal_s arg = {n, alpha, x};
    isched_parallel_call( pastix_data->isched, pthread_bvec_zscal, &arg );
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Compute y <- alpha * x + y. (Sequential version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
-*
 * @param[in] n
 *          The size of the vectors.
 *
 * @param[in] alpha
 *          A scalar.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[inout] y
 *          The vector y.
 *
 *******************************************************************************/
void
bvec_zaxpy_seq( pastix_data_t            *pastix_data,
                pastix_int_t              n,
                pastix_complex64_t        alpha,
                const pastix_complex64_t *x,
                pastix_complex64_t       *y)
{
#if defined(PASTIX_WITH_MPI) && 0
    SolverMatrix  *solvmtx = pastix_data->solvmatr;
    SolverCblk    *scblk   = solvmtx->cblktab;
    pastix_bcsc_t *bcsc    = pastix_data->bcsc;
    bcsc_cblk_t   *bcblk   = bcsc->cscftab;
    pastix_int_t   i, cblknbr;

    cblknbr = bcsc->cscfnbr;
    for( i = 0; i < cblknbr; i++, bcblk++ ){
        scblk   = solvmtx->cblktab + bcblk->cblknum;
        n = cblk_colnbr( scblk );

        cblas_zaxpy( n, CBLAS_SADDR(alpha),
                     x + scblk->lcolidx, 1,
                     y + scblk->lcolidx, 1 );
    }
#else
    (void)pastix_data;
    cblas_zaxpy( n, CBLAS_SADDR(alpha), x, 1, y, 1 );
#endif
}

struct z_argument_axpy_s
{
    pastix_int_t              n;
    pastix_complex64_t        alpha;
    const pastix_complex64_t *x;
    pastix_complex64_t       *y;
};

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Compute y <- alpha * x + y (Parallel version).
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          The context about the current thread
 *
 * @param[inout] args
 *          The parameter providing the size of the vectors, a scalar,
 *          the vectors x and y.
 *
 *******************************************************************************/
static inline void
pthread_bvec_zaxpy( isched_thread_t *ctx,
                    void            *args)
{
    struct z_argument_axpy_s  *arg = (struct z_argument_axpy_s*)args;
    pastix_int_t               n = arg->n;
    pastix_complex64_t         alpha = arg->alpha;
    const pastix_complex64_t  *x = arg->x;
    pastix_complex64_t        *y = arg->y;
    pastix_int_t               rank, size;
    pastix_int_t               begin, end;

    if( (y == NULL) || (x == NULL) ) {
        return;
    }

    if( alpha == (pastix_complex64_t)0.0 ) {
        return;
     }

    size = (pastix_int_t)ctx->global_ctx->world_size;
    rank = (pastix_int_t)ctx->rank;

    begin = (n/size) * rank;
    if (rank  == (size - 1)) {
        end = n;
    } else {
        end = (n/size) * (rank + 1);
    }

    if ( (end - begin) > 0 ) {
        cblas_zaxpy( end - begin, CBLAS_SADDR(alpha), x + begin, 1, y + begin, 1 );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Perform y = alpha * x + y (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
 *
 * @param[in] n
 *          The number of elements of vectors x and y
 *
 * @param[in] alpha
 *          The scalar to scale x
 *
 * @param[in] x
 *          The vector to be scaled
 *
 * @param[inout] y
 *          The resulting solution
 *
 *******************************************************************************/
void
bvec_zaxpy_smp( pastix_data_t            *pastix_data,
                pastix_int_t              n,
                pastix_complex64_t        alpha,
                const pastix_complex64_t *x,
                pastix_complex64_t       *y )
{
    struct z_argument_axpy_s args = {n, alpha, x, y};
    isched_parallel_call( pastix_data->isched, pthread_bvec_zaxpy, &args );
}

struct z_argument_dotc_s
{
    pastix_int_t              n;
    const pastix_complex64_t *x;
    const pastix_complex64_t *y;
    pastix_atomic_lock_t      lock;
    pastix_complex64_t        sum;
};

#if defined(PRECISION_z) || defined(PRECISION_c)
/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Compute the scalar product x.conj(y). (Sequential version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
 *
 * @param[in] n
 *          The size of the vectors.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] y
 *          The vector y.
 *
 *******************************************************************************
 *
 * @retval the scalar product of x and conj(y).
 *
 *******************************************************************************/
pastix_complex64_t
bvec_zdotc_seq( pastix_data_t            *pastix_data,
                pastix_int_t              n,
                const pastix_complex64_t *x,
                const pastix_complex64_t *y )
{
    SolverMatrix             *solvmtx = pastix_data->solvmatr;
    SolverCblk               *scblk   = solvmtx->cblktab;
    pastix_bcsc_t            *bcsc    = pastix_data->bcsc;
    bcsc_cblk_t              *bcblk   = bcsc->cscftab;
    pastix_int_t              i, j, cblknbr;
    const pastix_complex64_t *xptr;
    const pastix_complex64_t *yptr;
    pastix_complex64_t        r = 0.0;

    cblknbr = bcsc->cscfnbr;
    for( i = 0; i < cblknbr; i++, bcblk++ ) {
        scblk  = solvmtx->cblktab + bcblk->cblknum;
        n = cblk_colnbr( scblk );

        xptr = x + scblk->lcolidx;
        yptr = y + scblk->lcolidx;
        for( j=0; j<n; j++, xptr++, yptr++ ) {
            r += (*xptr) * conj(*yptr);
        }
    }

#if defined(PASTIX_WITH_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &r, 1, PASTIX_MPI_COMPLEX64,
                   MPI_SUM, solvmtx->solv_comm );
#endif

    return r;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Compute the scalar product x.conj(y). (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          The context of the current thread
 *
 * @param[in] args
 *          The argument provding the vectors x and y and their size.
 *
 *******************************************************************************/
static inline void
pthread_bvec_zdotc( isched_thread_t *ctx,
                    void            *args )
{
    struct z_argument_dotc_s *arg  = (struct z_argument_dotc_s*)args;
    pastix_int_t n = arg->n;
    int i;
    const pastix_complex64_t *xptr = arg->x;
    const pastix_complex64_t *yptr = arg->y;
    pastix_int_t begin, end, rank, size;
    pastix_complex64_t r = 0.0;

    rank = (pastix_int_t)ctx->rank;
    size = (pastix_int_t)ctx->global_ctx->world_size;

    begin = (n/size) * rank;
    if (rank != size - 1) {
        end = (n/size) * (rank + 1);
    } else { /*The last one computes the calcul for the rest of the sum*/
        end = n;
    }

    xptr += begin;
    yptr += begin;

    for ( i = begin; i < end; i++, xptr++, yptr++ )
    {
        r += (*xptr) * conj(*yptr);
    }

    if ( cabs(r) > 0. ) {
        pastix_atomic_lock( &(arg->lock) );
        arg->sum += r;
        pastix_atomic_unlock( &(arg->lock) );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Compute a scalar product between complex vectors: x.conj(y)
 * (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
 *
 * @param[in] n
 *          The number of elements of vectors x, y and r
 *
 * @param[in] y
 *          The first vector of the scalar product
 *
 * @param[in] n
 *          The second vector of the scalar product
 *
 * @param[out] r
 *          The result of the scalar product
 *
 *******************************************************************************/
pastix_complex64_t
bvec_zdotc_smp( pastix_data_t            *pastix_data,
                pastix_int_t              n,
                const pastix_complex64_t *x,
                const pastix_complex64_t *y )
{
    struct z_argument_dotc_s arg = {n, x, y, PASTIX_ATOMIC_UNLOCKED, 0.0};
    isched_parallel_call( pastix_data->isched, pthread_bvec_zdotc, &arg );

    return arg.sum;
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Compute the scalar product x.y. (Sequential version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] y
 *          The vector y.
 *
 * @param[in] n
 *          The size of the vectors.
 *
 *******************************************************************************
 *
 * @retval the scalar product of x and y.
 *
 *******************************************************************************/
pastix_complex64_t
bvec_zdotu_seq( pastix_data_t            *pastix_data,
                pastix_int_t              n,
                const pastix_complex64_t *x,
                const pastix_complex64_t *y )
{
    SolverMatrix             *solvmtx = pastix_data->solvmatr;
    SolverCblk               *scblk   = solvmtx->cblktab;
    pastix_bcsc_t            *bcsc    = pastix_data->bcsc;
    bcsc_cblk_t              *bcblk   = bcsc->cscftab;
    pastix_int_t              i, j, cblknbr;
    const pastix_complex64_t *xptr;
    const pastix_complex64_t *yptr;
    pastix_complex64_t        r = 0.0;

    cblknbr = bcsc->cscfnbr;
    for( i = 0; i < cblknbr; i++, bcblk++ ) {
        scblk  = solvmtx->cblktab + bcblk->cblknum;
        n = cblk_colnbr( scblk );

        xptr = x + scblk->lcolidx;
        yptr = y + scblk->lcolidx;
        for( j=0; j<n; j++, xptr++, yptr++ ) {
            r += (*xptr) * (*yptr);
        }
    }

#if defined(PASTIX_WITH_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &r, 1, PASTIX_MPI_COMPLEX64,
                   MPI_SUM, solvmtx->solv_comm );
#endif

    return r;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Compute the scalar product x.y. (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          The context of the current thread.
 *
 * @param[in] args
 *          The argument providing the vectors x and y and their size n.
 *
 *******************************************************************************/
static inline void
pthread_bvec_zdotu( isched_thread_t *ctx,
                    void            *args )
{
    struct z_argument_dotc_s *arg = (struct z_argument_dotc_s*)args;
    int i;
    pastix_int_t              n = arg->n;
    const pastix_complex64_t *x = arg->x;
    const pastix_complex64_t *y = arg->y;
    const pastix_complex64_t *xptr;
    const pastix_complex64_t *yptr;
    pastix_complex64_t        r = 0.0;
    pastix_int_t              size, rank;
    pastix_int_t              begin, end;

    rank = (pastix_int_t)ctx->rank;
    size = (pastix_int_t)ctx->global_ctx->world_size;

    begin = (n / size) * rank;
    if ( rank == (size - 1)) {
        end = n;
    } else {
        end = (n / size) * (rank + 1);
    }

    xptr = x + begin;
    yptr = y + begin;

    for (i=begin; i<end; i++, xptr++, yptr++)
    {
        r += (*xptr) * (*yptr);
    }

    if ( cabs(r) > 0. ) {
        pastix_atomic_lock( &(arg->lock) );
        arg->sum += r;
        pastix_atomic_unlock( &(arg->lock) );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Compute a regular scalar product x.y (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
 *
 * @param[in] n
 *          The number of elements of vectors x, y and r
 *
 * @param[in] x
 *          The first vector of the scalar product
 *
 * @param[in] y
 *          The second vector of the scalar product
 *
 * @param[out] r
 *          The result of the scalar product
 *
 *******************************************************************************
 *
 * @return The allocated vector
 *
 *******************************************************************************/
pastix_complex64_t
bvec_zdotu_smp( pastix_data_t            *pastix_data,
                pastix_int_t              n,
                const pastix_complex64_t *x,
                const pastix_complex64_t *y )
{
    struct z_argument_dotc_s arg = {n, x, y, PASTIX_ATOMIC_UNLOCKED, 0.0};
    isched_parallel_call( pastix_data->isched, pthread_bvec_zdotu, &arg );

    return arg.sum;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Apply a row permutation to a matrix A (LAPACK xlatmr)
 *
 *******************************************************************************
 *
 * @param[in] thread_safe
 *          Boolean to switch between the thread-safe implementation that
 *          exploits an additional workspace, or the non thread-safe version
 *          that has no memory overhead.
 *
 * @param[in] m
 *          The number of rows in the matrix A, and the number of elements in
 *          perm.
 *
 * @param[in] n
 *          The number of columns in the matrix A.
 *
 * @param[inout] A
 *          A matrix of size lda-by-n.
 *          On exit, rowas are permuted and A contains P A.
 *
 * @param[in] lda
 *          The leading dimension of A.
 *
 * @param[inout] perm
 *          The permutation array. Must be 0 based.
 *          If thread_safe is true, then perm array is used only as input, and a
 *          temporary array is allocated to follow the cycles. If thread_safe is
 *          false, then perm array is modified during the swap and restored at
 *          the end of the call.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
bvec_zlapmr( int thread_safe,
             pastix_dir_t        dir,
             pastix_int_t        m,
             pastix_int_t        n,
             pastix_complex64_t *A,
             pastix_int_t        lda,
             pastix_int_t       *perm )
{
    pastix_complex64_t tmp;
    pastix_int_t i, j, k, jj;
    pastix_int_t *perm_cpy;

    if ( thread_safe ) {
        perm_cpy = malloc( m * sizeof(pastix_int_t) );
        memcpy( perm_cpy, perm, m * sizeof(pastix_int_t) );
    }
    else {
        perm_cpy = perm;
    }

    if ( dir == PastixDirBackward ) {
        for(k=0; k<m; k++) {
            i = k;
            j = perm_cpy[i];

            /* Cycle already seen */
            if ( j < 0 ) {
                continue;
            }

            /* Mark the i^th element as being seen */
            perm_cpy[i] = -j-1;

            while( j != k ) {

                for(jj=0; jj<n; jj++) {
                    tmp             = A[j + jj * lda];
                    A[j + jj * lda] = A[k + jj * lda];
                    A[k + jj * lda] = tmp;
                }

                i = j;
                j = perm_cpy[i];
                perm_cpy[i] = -j-1;

                assert( (j != i) && (j >= 0) );
            }
        }
    }
    else {
        for(k=0; k<m; k++) {
            i = k;
            j = perm_cpy[i];
            perm_cpy[i] = -j-1;

            /* Cycle already seen */
            if ( j < 0 ) {
                continue;
            }

            i = perm_cpy[j];

            /* Mark the i^th element as being seen */
            while( i >= 0 ) {

                for(jj=0; jj<n; jj++) {
                    tmp             = A[j + jj * lda];
                    A[j + jj * lda] = A[i + jj * lda];
                    A[i + jj * lda] = tmp;
                }

                perm_cpy[j] = -i-1;
                j = i;
                i = perm_cpy[j];

                assert( j != i );
            }
        }
    }

    if ( thread_safe ) {
        free( perm_cpy );
    }
    else {
        /* Restore perm array */
        for(k=0; k<m; k++) {
            assert(perm[k] < 0);
            perm[k] = - perm[k] - 1;
        }
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Copy a vector y = x (Sequential version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
 *
 * @param[in] n
 *          The number of elements of vectors x and y
 *
 * @param[in] x
 *          The vector to be copied
 *
 * @param[inout] y
 *          The vector copy of x
 *
 *******************************************************************************/
void
bvec_zcopy_seq( pastix_data_t            *pastix_data,
                pastix_int_t              n,
                const pastix_complex64_t *x,
                pastix_complex64_t       *y )
{
#if defined(PASTIX_WITH_MPI) && 0
    SolverMatrix  *solvmtx = pastix_data->solvmatr;
    SolverCblk    *scblk   = solvmtx->cblktab;
    pastix_bcsc_t *bcsc    = pastix_data->bcsc;
    bcsc_cblk_t   *bcblk   = bcsc->cscftab;
    pastix_int_t   i, cblknbr;

    cblknbr = bcsc->cscfnbr;
    for( i = 0; i < cblknbr; i++, bcblk++ ) {
        scblk  = solvmtx->cblktab + bcblk->cblknum;
        n = cblk_colnbr( scblk );

        memcpy( y + scblk->lcolidx, x + scblk->lcolidx, n * sizeof(pastix_complex64_t) );
    }
#else
    (void)pastix_data;
    memcpy( y, x, n * sizeof(pastix_complex64_t) );
#endif
}

struct argument_copy_s {
    pastix_int_t              n;
    const pastix_complex64_t *x;
    pastix_complex64_t       *y;
};

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Copy a vector y = x (parallel version)
 *
 *        Perform the coopy of x into y.
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          The context of the current thread
 *
 * @param[inout] args
 *          The argument containing the vector copy and the one to copy and their
 *          size.
 *
 *******************************************************************************/
static inline void
pthread_bvec_zcopy( isched_thread_t *ctx,
                    void            *args )
{
    struct argument_copy_s   *arg = (struct argument_copy_s*)args;
    pastix_int_t              n = arg->n;
    pastix_int_t              size, rank;
    pastix_int_t              begin, end;

    size = (pastix_int_t)ctx->global_ctx->world_size;
    rank = (pastix_int_t)ctx->rank;

    begin = (n/size) * rank;

    if (rank == (size - 1)) {
        end = n;
    } else {
        end = (n/size) * (rank + 1);
    }

    if ( (end - begin) > 0 ) {
        memcpy( arg->y + begin, arg->x + begin, (end - begin) * sizeof(pastix_complex64_t) );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Copy a vector y = x (parallel version)
 *
 *        Initialise argument for the parallel function
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
 *
 * @param[in] n
 *          The number of elements of vectors x and y
 *
 * @param[in] x
 *          The vector to be copied
 *
 * @param[inout] y
 *          The vector copy of x
 *
 *******************************************************************************/
void
bvec_zcopy_smp( pastix_data_t            *pastix_data,
                pastix_int_t              n,
                const pastix_complex64_t *x,
                pastix_complex64_t       *y )
{
    struct argument_copy_s args = {n, x, y};
    isched_parallel_call( pastix_data->isched, pthread_bvec_zcopy, &args );
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Solve A x = b with A the sparse matrix
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 * @param[inout] d
 *          On entry, the right hand side
 *          On exit, the solution of tha problem A x = b
 *
 *******************************************************************************/
void bcsc_zspsv( pastix_data_t      *pastix_data,
                 pastix_complex64_t *b )
{
    pastix_int_t n = pastix_data->bcsc->gN;
    pastix_data->iparm[IPARM_VERBOSE]--;
    pastix_subtask_solve( pastix_data, 1, b, n );
    pastix_data->iparm[IPARM_VERBOSE]++;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Compute \f[ y = \alpha A x + \beta y \f] (Sequential version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
 *
 * @param[in] m
 *          The number of rows of the matrix A, and the size of y.
 *
 * @param[in] n
 *          The number of columns of the matrix A, and the size of x.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The dense matrix A of size lda-by-n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1,m)
 *
 * @param[in] x
 *          The vector x of size n.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[inout] y
 *          On entry, the initial vector y of size m.
 *          On exit, the updated vector.
 *
 *******************************************************************************/
void
bvec_zgemv_seq( pastix_data_t            *pastix_data,
                pastix_int_t              m,
                pastix_int_t              n,
                pastix_complex64_t        alpha,
                const pastix_complex64_t *A,
                pastix_int_t              lda,
                const pastix_complex64_t *x,
                pastix_complex64_t        beta,
                pastix_complex64_t       *y )
{
#if defined(PASTIX_WITH_MPI) && 0
    SolverMatrix  *solvmtx = pastix_data->solvmatr;
    SolverCblk    *scblk   = solvmtx->cblktab;
    pastix_bcsc_t *bcsc    = pastix_data->bcsc;
    bcsc_cblk_t   *bcblk   = bcsc->cscftab;
    pastix_int_t   i, cblknbr;

    cblknbr = bcsc->cscfnbr;
    for( i = 0; i < cblknbr; i++, bcblk++ ) {
        scblk  = solvmtx->cblktab + bcblk->cblknum;
        m = cblk_colnbr( scblk );

        cblas_zgemv( CblasColMajor, CblasNoTrans, m, n,
                     CBLAS_SADDR(alpha), A + scblk->lcolidx, lda, x, 1,
                     CBLAS_SADDR(beta), y + scblk->lcolidx, 1 );
    }
#else
    (void)pastix_data;
    cblas_zgemv( CblasColMajor, CblasNoTrans, m, n,
                 CBLAS_SADDR(alpha), A, lda, x, 1,
                 CBLAS_SADDR(beta), y, 1 );
#endif

}

struct z_gemv_s
{
    pastix_int_t              m;
    pastix_int_t              n;
    pastix_complex64_t        alpha;
    const pastix_complex64_t *A;
    pastix_int_t              lda;
    const pastix_complex64_t *x;
    pastix_complex64_t        beta;
    pastix_complex64_t       *y;
};

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Compute \f[ y = \alpha A x + \beta y \f] (Parallel version)
 *
 *        This is the function called by bvec_zgemv_smp to perform gemv
 *        (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          Information about number of thread and rank of current thread.
 *
 * @param[inout] args
 *          The number of rows (m) and columns (n) of the matrix A, vecteors y
 *          (size m) and x (size n), scalars alpha and beta. The dense matrix A
 *          of size lda-by-n, and its leading dimension lda >= max(1,m), and
 *          the vector x of size n. The result is stored in y.
 *
 *******************************************************************************/
static inline void
pthread_bvec_zgemv( isched_thread_t *ctx,
                    void            *args )
{
    struct z_gemv_s          *arg = (struct z_gemv_s*)args;
    pastix_int_t              m = arg->m, sub_m;
    pastix_int_t              n = arg->n;
    pastix_complex64_t        alpha = arg->alpha;
    const pastix_complex64_t *A = arg->A, *Aptr;
    pastix_int_t              lda = arg->lda;
    const pastix_complex64_t *x = arg->x, *xptr;
    pastix_complex64_t        beta = arg->beta;
    pastix_complex64_t       *y = arg->y, *yptr;
    pastix_int_t              size, rank;

    size = (pastix_int_t)ctx->global_ctx->world_size;
    rank = (pastix_int_t)ctx->rank;

    sub_m = (m / size);

    /* Subdivision of A */
    Aptr = A + sub_m * rank;

    xptr = x;
    yptr = y + sub_m * rank;

    if (rank == (size - 1)) {
        sub_m += m % size; /* Last thread has to do more tasks */
    }

    if ( sub_m > 0 ) {
        cblas_zgemv( CblasColMajor, CblasNoTrans, sub_m, n,
                     CBLAS_SADDR(alpha), Aptr, lda, xptr, 1,
                     CBLAS_SADDR(beta), yptr, 1 );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Compute \f[ y = \alpha A x + \beta y \f] (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
 *
 * @param[in] m
 *          The number of rows of the matrix A, and the size of y.
 *
 * @param[in] n
 *          The number of columns of the matrix A, and the size of x.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The dense matrix A of size lda-by-n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1,m)
 *
 * @param[in] x
 *          The vector x of size n.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[inout] y
 *          On entry, the initial vector y of size m.
 *          On exit, the updated vector.
 *
 *******************************************************************************/
void
bvec_zgemv_smp( pastix_data_t            *pastix_data,
                pastix_int_t              m,
                pastix_int_t              n,
                pastix_complex64_t        alpha,
                const pastix_complex64_t *A,
                pastix_int_t              lda,
                const pastix_complex64_t *x,
                pastix_complex64_t        beta,
                pastix_complex64_t       *y )
{
    struct z_gemv_s arg = {m, n, alpha, A, lda, x, beta, y};

    isched_parallel_call( pastix_data->isched, pthread_bvec_zgemv, &arg );
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Set to 0 remote coefficients
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
 *
 * @param[inout] y
 *          On entry, the initial vector y of size m.
 *          On exit, the y vector with remote section set to 0.
 *
 *******************************************************************************/
void
bvec_znullify_remote( const pastix_data_t *pastix_data,
                      pastix_complex64_t  *y )
{
#if defined( PASTIX_WITH_MPI )
    const SolverMatrix *solvmtx = pastix_data->solvmatr;
    const SolverCblk   *cblk = solvmtx->cblktab;
    pastix_int_t        cblknbr;
    pastix_int_t        i, lastindex = 0;

    cblknbr = solvmtx->cblknbr;
    for ( i = 0; i < cblknbr; i++, cblk++ ) {
        if ( cblk->cblktype & (CBLK_FANIN|CBLK_RECV) ) {
            continue;
        }

        if ( cblk->fcolnum != lastindex ) {
            /* Set to 0 all remote data bewtween previous local cblk, and current cblk */
            memset(
                y + lastindex, 0, ( cblk->fcolnum - lastindex ) * sizeof( pastix_complex64_t ) );
        }
        lastindex = cblk->lcolnum + 1;
    }

    if ( lastindex < solvmtx->nodenbr ) {
        /* Set to 0 all remote data bewtween previous local cblk, and current cblk */
        memset( y + lastindex, 0, ( solvmtx->nodenbr - lastindex ) * sizeof( pastix_complex64_t ) );
    }
#endif
    (void)pastix_data;
    (void)y;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Apply an all reduce of the vector on all nodes
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The information about sequential and parallel version (Number of
 *          thread, ...).
 *
 * @param[inout] y
 *          On entry, the initial vector y of size m.
 *          On exit, the y vector with remote section set to 0.
 *
 *******************************************************************************/
void
bvec_zallreduce( const pastix_data_t *pastix_data,
                 pastix_complex64_t  *y )
{
#if defined( PASTIX_WITH_MPI )
    /* Reduce the partial sums on all nodes */
    MPI_Allreduce( MPI_IN_PLACE,
                   y, pastix_data->bcsc->gN,
                   PASTIX_MPI_COMPLEX64, MPI_SUM,
                   pastix_data->inter_node_comm );
#endif
    (void)pastix_data;
    (void)y;
}
