/**
 *
 * @file starpu.c
 *
 * @copyright 2017-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.1.0
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2020-01-06
 *
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#include "common.h"
#if !defined(PASTIX_WITH_STARPU)
#error "This file should not be compiled if Starpu is not enabled"
#endif
#include <stdio.h>
#include "pastix_starpu.h"

/**
 *  Set the tag sizes
 */
#if defined(PASTIX_WITH_MPI)

/* Take 24 bits for the tile id, and 7 bits for descriptor id.
 These values can be changed through the call CHAMELEON_user_tag_size(int tag_width, int tag_sep) */
#define TAG_WIDTH_MIN 20
static int starpu_tag_width = 64;
static int starpu_tag_sep   = 32;
static volatile int32_t starpu_tag_counter = 0;
static int _tag_mpi_initialized_ = 0;

/**
 *******************************************************************************
 *
 * @brief Init a new StarPU tag.
 *
 * This function initializes the StarPU tags thanks to global values
 *
 *******************************************************************************
 *
 * @param[inout] pastix
 *          The main pastix_data structure.
 *
 ******************************************************************************/
static int
pastix_starpu_tag_init( pastix_data_t *pastix )
{
    if (!_tag_mpi_initialized_) {
        int ok = 0;
        uintptr_t tag_ub;

        void *tag_ub_p = NULL;

        starpu_tag_width = 64;
        starpu_tag_sep   = 24;

#if defined(HAVE_STARPU_MPI_COMM_GET_ATTR) || 1
        starpu_mpi_comm_get_attr( pastix->inter_node_comm, STARPU_MPI_TAG_UB, &tag_ub_p, &ok );
#else
        MPI_Comm_get_attr( pastix->inter_node_comm, MPI_TAG_UB, &tag_ub_p, &ok );
#endif
        tag_ub = (uintptr_t)tag_ub_p;

        if ( !ok ) {
            pastix_print_error("pastix_starpu_tag_init: MPI_TAG_UB not known by StarPU\n");
        }

        while ( ((uintptr_t)((1UL<<starpu_tag_width) - 1) > tag_ub ) &&
                (starpu_tag_width >= TAG_WIDTH_MIN) )
        {
            starpu_tag_width--;
            starpu_tag_sep = starpu_tag_width / 2;
        }

        if ( starpu_tag_width < TAG_WIDTH_MIN ) {
            pastix_print_error("pastix_starpu_tag_init: MPI_TAG_UB may be too small to identify all the pieces of data\n");
            return PASTIX_ERR_INTERNAL;
        }

        _tag_mpi_initialized_ = 1;
        return PASTIX_SUCCESS;
    }
    else {
        return PASTIX_ERR_INTERNAL;
    }
}

/**
 *******************************************************************************
 *
 * @brief Get the StarPU unique current tag.
 *
 * This function returns a new tag for the StarPU distributed version thanks to
 * global variables.
 *
 ******************************************************************************/
int64_t
pastix_starpu_get_tag( )
{
    return ((int64_t)pastix_atomic_inc_32b( &starpu_tag_counter )) << starpu_tag_sep;
}
#else /* defined(PASTIX_WITH_MPI) */
int64_t
pastix_starpu_get_tag(){
    return 1;
}
#endif

/**
 *******************************************************************************
 *
 * @brief Startup the StarPU runtime system.
 *
 * This function initialize and startup the StarPU runtime system with PaStix
 * configuration variables
 *
 *******************************************************************************
 *
 * @param[inout] pastix
 *          The main pastix_data structure.
 *
 * @param[inout] argc
 *          The number of arguments of the main program.
 *
 * @param[inout] argv
 *          The list of argument given to the main program.
 *
 * @param[in] bindtab
 *          The binding array of size the number of threads if a specific
 *          binding is required, NULL otherwise.
 *
 ******************************************************************************/
void
pastix_starpu_init( pastix_data_t *pastix,
                    int *argc, char **argv[],
                    const int *bindtab )
{
    struct starpu_conf *conf;
    pastix_int_t *iparm = pastix->iparm;
    int rc;

    if ( pastix->starpu != NULL )
        return;

    pastix->starpu = malloc(sizeof(struct starpu_conf));
    starpu_conf_init( pastix->starpu );

    /* Force no GPUs if CUDA has not been enabled in PaStiX */
#if !defined(PASTIX_WITH_CUDA)
    iparm[IPARM_GPU_NBR] = 0;
#endif

    conf = pastix->starpu;
    conf->ncpus = iparm[IPARM_THREAD_NBR];
    conf->ncuda = iparm[IPARM_GPU_NBR];
    conf->nopencl = 0;

    if (conf->ncuda > 0) {
#if defined(PASTIX_GENERATE_MODEL)
        pastix_print( pastix->procnum, 0,
                      "WARNING: PaStiX compiled with -DPASTIX_GENERATE_MODEL forces:\n"
                      "    - a single event per stream\n"
                      "    - a single stream per GPU\n"
                      "    - restore the automatic detection of the number of threads\n" );

        conf->ncpus = -1;

        pastix_setenv( "STARPU_NWORKER_PER_CUDA", "1", 1 );
        pastix_setenv( "STARPU_CUDA_PIPELINE", "0", 1 );
#endif

        conf->sched_policy_name = "dmdas";
    }
    else {
        /*
         * Set scheduling to "ws"/"lws" if no cuda devices used because it
         * behaves better on homogneneous architectures. If the user wants
         * to use another scheduling strategy, he can set STARPU_SCHED
         * env. var. to whatever he wants
         */
#if (STARPU_MAJOR_VERSION > 1) || ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION >= 2))
        conf->sched_policy_name = "lws";
#else
        conf->sched_policy_name = "ws";
#endif
    }

    if ( bindtab != NULL ) {
        int i;

        assert( iparm[IPARM_THREAD_NBR] < STARPU_NMAXWORKERS );
        conf->use_explicit_workers_bindid = 1;

        for(i=0; i < pastix_imin( iparm[IPARM_THREAD_NBR], STARPU_NMAXWORKERS ); i++) {
            conf->workers_bindid[i] = bindtab[i];
        }
    }
    rc = starpu_init( conf );

    starpu_malloc_on_node_set_default_flags( STARPU_MAIN_RAM,
                                             STARPU_MALLOC_PINNED
                                             | STARPU_MALLOC_COUNT
#if defined(PASTIX_STARPU_SIMULATION)
                                             | STARPU_MALLOC_SIMULATION_FOLDED
#endif
                                             );

#if defined(PASTIX_WITH_MPI)
#if defined(PASTIX_DEBUG_MPI) && !defined(PASTIX_STARPU_SIMULATION)
    {
        int flag = 0;
        MPI_Initialized( &flag );
        assert( flag );
    }
#endif
    starpu_mpi_init_comm( argc, argv, 0, pastix->inter_node_comm );
    pastix_starpu_tag_init( pastix );
#endif

#if defined(PASTIX_WITH_CUDA) && !defined(PASTIX_STARPU_SIMULATION)
    starpu_cublas_init();
#endif

    /* Suspend threads until we need them */
    starpu_pause();

    assert( pastix->starpu != NULL );

    (void)argc; (void)argv;
    (void)rc;
}

/**
 *******************************************************************************
 *
 * @brief Finalize the StarPU runtime system.
 *
 * This function stop the StarPU runtime system.
 *
 *******************************************************************************
 *
 * @param[inout] pastix
 *          The main pastix_data structure.
 *
 ******************************************************************************/
void
pastix_starpu_finalize( pastix_data_t *pastix )
{
    if (pastix->starpu != NULL) {
        starpu_resume();

#if defined(PASTIX_WITH_MPI)
        starpu_mpi_shutdown();
#endif
#if defined(PASTIX_WITH_CUDA) && !defined(PASTIX_STARPU_SIMULATION)
        starpu_cublas_shutdown();
#endif
        starpu_shutdown();

        free( pastix->starpu );
        pastix->starpu = NULL;
    }
}

/**
 * @}
 */
