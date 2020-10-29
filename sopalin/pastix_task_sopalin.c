/**
 *
 * @file pastix_task_sopalin.c
 *
 *  PaStiX factorization routines
 *
 * @copyright 2004-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2020-01-21
 *
 **/
#include "common.h"
#include <lapacke.h>
#include "isched.h"
#include "spm.h"
#include "bcsc.h"
#include "blend/solver.h"
#include "coeftab.h"
#include "sopalin_data.h"
#include "kernels/pastix_lowrank.h"
#include "kernels/pastix_zlrcores.h"
#include "kernels/pastix_clrcores.h"
#include "kernels/pastix_dlrcores.h"
#include "kernels/pastix_slrcores.h"
#include "kernels/kernels_trace.h"

#if defined(PASTIX_WITH_PARSEC)
#include "sopalin/parsec/pastix_parsec.h"
#endif

#if defined(PASTIX_WITH_STARPU)
#include "sopalin/starpu/pastix_starpu.h"
#endif

static void (*sopalinFacto[5][4])(pastix_data_t *, sopalin_data_t*) =
{
    { sopalin_spotrf, sopalin_dpotrf, sopalin_cpotrf, sopalin_zpotrf },
    { sopalin_ssytrf, sopalin_dsytrf, sopalin_csytrf, sopalin_zsytrf },
    { sopalin_sgetrf, sopalin_dgetrf, sopalin_cgetrf, sopalin_zgetrf },
    { sopalin_spotrf, sopalin_dpotrf, sopalin_cpxtrf, sopalin_zpxtrf },
    { sopalin_ssytrf, sopalin_dsytrf, sopalin_chetrf, sopalin_zhetrf }
};

/**
 *******************************************************************************
 *
 * @ingroup pastix_numfact
 *
 * @brief Fill the internal block CSC structure.
 *
 * This internal block CSC can be used by the refinement step with or without
 * the preconditioner obtained by the numerical factorization.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the internal block CSC is filled with entries from
 *          the spm matrix.
 *
 * @param[inout] spm
 *          The sparse matrix descriptor that describes problem instance.
 *          On exit, if IPARM_FREE_CSCUSER is set, the spm is freed.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastix_subtask_spm2bcsc( pastix_data_t *pastix_data,
                         spmatrix_t    *spm )
{
    double time;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_spm2bcsc: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (spm == NULL) {
        errorPrint("pastix_subtask_spm2bcsc: wrong spm parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_ANALYSE) ) {
        errorPrint("pastix_subtask_spm2bcsc: All steps from pastix_task_init() to pastix_task_blend() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    /*
     * Compute the norm of A, to scale the epsilon parameter for pivoting
     */
    {
        pastix_data->dparm[ DPARM_A_NORM ] = spmNorm( SpmFrobeniusNorm, spm );
        if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNo ) {
            pastix_print( pastix_data->inter_node_procnum, 0,
                          "    ||A||_2  =                            %e\n",
                          pastix_data->dparm[ DPARM_A_NORM ] );
        }
    }

    /*
     * Temporary switch of the solver pointer for the runtime
     */
    if ( !(pastix_data->steps & STEP_CSC2BCSC) )
    {
        SolverMatrix *tmp;
        switch( pastix_data->iparm[IPARM_SCHEDULER] ){
        case PastixSchedParsec:
        case PastixSchedStarPU:
            tmp = pastix_data->solvmatr;
            pastix_data->solvmatr = pastix_data->solvglob;
            pastix_data->solvglob = tmp;
            break;
        default:
            break;
        }
    }

    /*
     * Fill in the internal blocked CSC. We consider that if this step is called
     * the spm values have changed so we need to update the blocked csc.
     */
    if ( pastix_data->bcsc != NULL )
    {
        bcscExit( pastix_data->bcsc );
        memFree_null( pastix_data->bcsc );
    }

    MALLOC_INTERN( pastix_data->bcsc, 1, pastix_bcsc_t );

    time = bcscInit( spm,
                     pastix_data->ordemesh,
                     pastix_data->solvmatr,
                     (pastix_data->iparm[IPARM_FACTORIZATION] == PastixFactLU), /*&& (! pastix_data->iparm[IPARM_ONLY_REFINE]) )*/
                     pastix_data->bcsc );

    if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
        pastix_print( pastix_data->inter_node_procnum, 0, OUT_BCSC_TIME, time );
    }

    if ( pastix_data->iparm[IPARM_FREE_CSCUSER] ) {
        spmExit( spm );
    }

    /*
     * Invalidate following step, and add current step to the ones performed
     */
    pastix_data->steps &= ~STEP_BCSC2CTAB;
    pastix_data->steps |= STEP_CSC2BCSC;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_numfact
 *
 * @brief Fill the internal solver matrix structure.
 *
 * This step is linked with the pastix_subtask_sopalin() since this structure is
 * only used during the numerical factorization.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, [IPARM_FACTORIZATION.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the internal solver structure is filled with entries from
 *          the internal block CSC matrix.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastix_subtask_bcsc2ctab( pastix_data_t *pastix_data )
{
    pastix_bcsc_t *bcsc;
    pastix_lr_t   *lr;
    Clock timer;
    int mtxtype;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_bcsc2ctab: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_CSC2BCSC) ) {
        errorPrint("pastix_subtask_bcsc2ctab: All steps from pastix_task_init() to pastix_stask_blend() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (pastix_data->bcsc == NULL) {
        errorPrint("pastix_subtask_bcsc2ctab: wrong pastix_data->bcsc parameter");
        return PASTIX_ERR_BADPARAMETER;
    }

    clockSyncStart( timer, pastix_data->inter_node_comm );

    /* Initialize low-rank parameters */
    lr = &(pastix_data->solvmatr->lowrank);
    lr->compress_when       = pastix_data->iparm[IPARM_COMPRESS_WHEN];
    lr->compress_method     = pastix_data->iparm[IPARM_COMPRESS_METHOD];
    lr->compress_min_width  = pastix_data->iparm[IPARM_COMPRESS_MIN_WIDTH];
    lr->compress_min_height = pastix_data->iparm[IPARM_COMPRESS_MIN_HEIGHT];
    lr->compress_preselect  = pastix_data->iparm[IPARM_COMPRESS_PRESELECT];
    lr->use_reltol          = pastix_data->iparm[IPARM_COMPRESS_RELTOL];
    lr->tolerance           = pastix_data->dparm[DPARM_COMPRESS_TOLERANCE];

    pastix_lr_minratio      = pastix_data->dparm[DPARM_COMPRESS_MIN_RATIO];
    pastix_lr_ortho         = pastix_data->iparm[IPARM_COMPRESS_ORTHO];

    bcsc = pastix_data->bcsc;
    lr->core_ge2lr = ge2lrMethods[ pastix_data->iparm[IPARM_COMPRESS_METHOD] ][bcsc->flttype-2];
    lr->core_rradd = rraddMethods[ pastix_data->iparm[IPARM_COMPRESS_METHOD] ][bcsc->flttype-2];

    if ( pastix_data->iparm[IPARM_COMPRESS_WHEN] == PastixCompressWhenBegin ) {
        if ( lr->compress_preselect == -1 ) {
            lr->compress_preselect = 1;
        }
        core_get_rklimit = core_get_rklimit_begin;
    }
    else {
        if ( lr->compress_preselect == -1 ) {
            if ( lr->tolerance <= 1.e-10 ) {
                lr->compress_preselect = 0;
            }
            else {
                lr->compress_preselect = 1;
            }
        }
        core_get_rklimit = core_get_rklimit_end;
    }

    pastix_data->solvmatr->factotype = pastix_data->iparm[IPARM_FACTORIZATION];

    /*
     * Fill in the internal coeftab structure. We consider that if this step is
     * called the bcsc values have changed, or a factorization have already been
     * performed, so we need to update the coeftab arrays.
     */
    if (pastix_data->bcsc != NULL)
    {
        coeftabExit( pastix_data->solvmatr );
    }

    coeftabInit( pastix_data,
                 pastix_data->iparm[IPARM_FACTORIZATION] == PastixFactLU ? PastixLUCoef : PastixLCoef );

    switch( pastix_data->iparm[IPARM_FACTORIZATION] ) {
    case PastixFactLLH:
    case PastixFactLDLH:
        mtxtype = PastixHermitian;
        break;

    case PastixFactLLT:
    case PastixFactLDLT:
        mtxtype = PastixHermitian;
        break;

    case PastixFactLU:
    default:
        mtxtype = PastixGeneral;
    }

#if defined(PASTIX_WITH_PARSEC)
    if ( pastix_data->iparm[IPARM_SCHEDULER] == PastixSchedParsec )
    {
        /* Start PaRSEC if not already started */
        if (pastix_data->parsec == NULL) {
            int argc = 0;
            pastix_parsec_init( pastix_data, &argc, NULL, NULL );
        }
        /* Create the matrix descriptor */
        parsec_sparse_matrix_init( pastix_data->solvmatr,
                                   pastix_size_of( bcsc->flttype ), mtxtype,
                                   pastix_data->inter_node_procnbr,
                                   pastix_data->inter_node_procnum );
    }
#endif

#if defined(PASTIX_WITH_STARPU)
    if ( pastix_data->iparm[IPARM_SCHEDULER] == PastixSchedStarPU )
    {
        /* Create the matrix descriptor */
        starpu_sparse_matrix_init( pastix_data->solvmatr,
                                   pastix_size_of( bcsc->flttype ), mtxtype,
                                   pastix_data->inter_node_procnbr,
                                   pastix_data->inter_node_procnum );
    }
#endif

    clockSyncStop( timer, pastix_data->inter_node_comm );
    if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
        pastix_print( pastix_data->inter_node_procnum, 0, OUT_COEFTAB_TIME,
                      clockVal(timer) );
    }

    /* Invalidate following step, and add current step to the ones performed */
    pastix_data->steps &= ~STEP_NUMFACT;
    pastix_data->steps |= STEP_BCSC2CTAB;

    (void)mtxtype;
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_numfact
 *
 * @brief Factorize the given problem using Cholesky or LU decomposition.
 *
 * The user can call the pastix_task_solve() to obtain the solution.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION, IPARM_COMPRESS_WHEN.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the solver matrix structure stores the factorization of the
 *          given problem.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastix_subtask_sopalin( pastix_data_t *pastix_data )
{
    sopalin_data_t  sopalin_data;
    SolverBackup_t *sbackup;
    pastix_bcsc_t  *bcsc;
    PASTIX_Comm     pastix_comm;
    pastix_int_t   *iparm;
    double         *dparm;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_sopalin: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_ANALYSE) ) {
        errorPrint("pastix_subtask_sopalin: All steps from pastix_task_init() to pastix_task_analyze() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_CSC2BCSC) ) {
        errorPrint("pastix_subtask_sopalin: All steps from pastix_task_init() to pastix_task_analyze() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_BCSC2CTAB) ) {
        errorPrint("pastix_subtask_sopalin: All steps from pastix_task_init() to pastix_task_analyze() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    bcsc = pastix_data->bcsc;
    if (bcsc == NULL) {
        errorPrint("pastix_subtask_sopalin: wrong pastix_data_bcsc parameter");
        return PASTIX_ERR_BADPARAMETER;
    }

    pastix_comm = pastix_data->inter_node_comm;
    iparm = pastix_data->iparm;
    dparm = pastix_data->dparm;

    /* Prepare the sopalin_data structure */
    {
        double threshold;

        sopalin_data.solvmtx = pastix_data->solvmatr;

        /* TODO: might change the behavior: if the user wants a ratio of the norm, it could compute it himself */
        if ( dparm[ DPARM_EPSILON_MAGN_CTRL ] < 0. ) {
            threshold = - dparm[ DPARM_EPSILON_MAGN_CTRL ];
        }
        else if ( dparm[ DPARM_EPSILON_MAGN_CTRL ] == 0. ) {
            /*
             * Use the rule presented in "Making Sparse Gaussian Elimination
             * Scalable by Static Pivoting", S. X. Li, J. Demmel
             *
             * sqrt(eps) * ||A||_1 is a small half precision pertrbation
             * intriduced in the pb when the pivot is too small
             */
            double eps;
            if ( (bcsc->flttype == PastixFloat) || (bcsc->flttype == PastixComplex32) ) {
                eps = LAPACKE_slamch_work( 'e' );
            }
            else {
                eps = LAPACKE_dlamch_work( 'e' );
            }
            threshold =  sqrt(eps) * dparm[DPARM_A_NORM];
        }
        else {
            threshold = dparm[ DPARM_EPSILON_MAGN_CTRL ] * dparm[DPARM_A_NORM];
        }

        sopalin_data.solvmtx->diagthreshold = threshold;
        sopalin_data.solvmtx->nbpivots      = 0;

        sopalin_data.cpu_coefs = &(pastix_data->cpu_models->coefficients[bcsc->flttype-2]);
        sopalin_data.gpu_coefs = &(pastix_data->gpu_models->coefficients[bcsc->flttype-2]);
    }

    sbackup = solverBackupInit( pastix_data->solvmatr );
    pastix_data->solvmatr->restore = 2;
    {
        void (*factofct)( pastix_data_t *, sopalin_data_t *);
        double timer, timer_local, flops, flops_g = 0.;

        factofct = sopalinFacto[ pastix_data->iparm[IPARM_FACTORIZATION] ][bcsc->flttype-2];
        assert(factofct);

        kernelsTraceStart( pastix_data );
        clockSyncStart( timer, pastix_comm );
        clockStart(timer_local);

        factofct( pastix_data, &sopalin_data );

        clockStop(timer_local);
        clockSyncStop( timer, pastix_comm );
        kernelsTraceStop( pastix_data );

        /* Output time and flops */
        pastix_data->dparm[DPARM_FACT_TIME] = clockVal(timer);
        flops = pastix_data->dparm[DPARM_FACT_THFLOPS] / pastix_data->dparm[DPARM_FACT_TIME];
        pastix_data->dparm[DPARM_FACT_FLOPS] = ((flops / 1024.) / 1024.) / 1024.;

        /* Compute the number of flops performed (may be different from THFLOPS for low-rank) */
        {
            double flops_l = pastix_data->dparm[DPARM_FACT_RLFLOPS];
            MPI_Reduce( &flops_l, &flops_g, 1, MPI_DOUBLE, MPI_SUM, 0, pastix_data->inter_node_comm );
        }

        /* Reduce the number of pivots */
        {
            int nbpivot_l = sopalin_data.solvmtx->nbpivots;
            int nbpivot_g;
            MPI_Allreduce( &nbpivot_l, &nbpivot_g, 1, MPI_INT, MPI_SUM, pastix_data->inter_node_comm );
            pastix_data->iparm[IPARM_STATIC_PIVOTING] = nbpivot_g;
        }

        if ( iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
            pastix_print( pastix_data->inter_node_procnum, 0, OUT_SOPALIN_TIME,
                          clockVal(timer),
                          pastix_print_value( flops ),
                          pastix_print_unit(  flops ),
                          pastix_print_value( flops_g ),
                          pastix_print_unit(  flops_g ),
                          (long)pastix_data->iparm[IPARM_STATIC_PIVOTING] );
        }

#if defined(PASTIX_WITH_PARSEC) && defined(PASTIX_DEBUG_PARSEC)
        {
            int i;
            fprintf(stderr, "-- Check status of PaRSEC nbtasks counters\n" );
            for (i=0; i<32; i++) {
                if ( parsec_nbtasks_on_gpu[i] != 0 ) {
                    fprintf(stderr, "Device %d => %d\n",
                            i, parsec_nbtasks_on_gpu[i] );
                }
                parsec_nbtasks_on_gpu[i] = 0;
            }
        }
#endif /* defined(PASTIX_WITH_PARSEC) && defined(PASTIX_DEBUG_PARSEC) */
    }
    solverBackupRestore( pastix_data->solvmatr, sbackup );
    solverBackupExit( sbackup );

#if defined(PASTIX_NUMFACT_DUMP_SOLVER)
    {
        FILE *stream = NULL;

        pastix_gendirectories( pastix_data );
        stream = pastix_fopenw( pastix_data->dir_local, "solver.eps", "w" );
        if ( stream ) {
            solverDraw( pastix_data->solvmatr,
                        stream, iparm[IPARM_VERBOSE],
                        pastix_data->dir_local );
            fclose(stream);
        }
    }
#endif

    if ( (pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot) &&
         (pastix_data->iparm[IPARM_COMPRESS_WHEN] != PastixCompressNever) )
    {
        /* Compute the memory gain */
        coeftabMemory[bcsc->flttype-2]( pastix_data->solvmatr );
    }

    /* Invalidate following steps, and add factorization step to the ones performed */
    pastix_data->steps &= ~( STEP_BCSC2CTAB |
                             STEP_SOLVE |
                             STEP_REFINE );
    pastix_data->steps |= STEP_NUMFACT;

    (void)pastix_comm;
    return EXIT_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_users
 *
 * @brief Perform all the numerical factorization steps:
 * fill the internal block CSC and the solver matrix structures,
 * then apply the factorization step.
 *
 * The user can call the pastix_task_solve() to obtain the solution.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the internal block CSC is filled with entries from the
 *          spm matrix and the solver matrix structure stores the factorization
 *          of the given problem.
 *
 * @param[inout] spm
 *          The sparse matrix descriptor that describes problem instance.
 *          On exit, the spm structure is freed if IPARM_FREE_CSCUSER is true.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastix_task_numfact( pastix_data_t *pastix_data,
                     spmatrix_t    *spm )
{
    pastix_int_t *iparm;
    pastix_int_t  procnum;
    pastix_int_t  steps;
    int rc;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_task_sopalin: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (spm == NULL) {
        errorPrint("pastix_task_sopalin: wrong spm parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_ANALYSE) ) {
        errorPrint("pastix_task_sopalin: All steps from pastix_task_init() to pastix_task_blend() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    iparm   = pastix_data->iparm;
    procnum = pastix_data->inter_node_procnum;

    if (iparm[IPARM_VERBOSE] > PastixVerboseNot) {
        pastix_print( procnum, 0, OUT_STEP_SOPALIN,
                      pastixFactotypeStr( iparm[IPARM_FACTORIZATION] ) );
    }

    /* Invalidate upcoming steps */
    steps = ~(STEP_CSC2BCSC  |
              STEP_BCSC2CTAB |
              STEP_NUMFACT   |
              STEP_SOLVE     |
              STEP_REFINE    );
    pastix_data->steps &= steps;

    rc = pastix_subtask_spm2bcsc( pastix_data, spm );
    if (rc != PASTIX_SUCCESS) {
        return rc;
    }

    rc = pastix_subtask_bcsc2ctab( pastix_data );
    if (rc != PASTIX_SUCCESS) {
        return rc;
    }

    rc = pastix_subtask_sopalin( pastix_data );
    if (rc != PASTIX_SUCCESS) {
        return rc;
    }

    return EXIT_SUCCESS;
}
