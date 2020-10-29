/**
 *
 * @file isched_hwloc.c
 *
 * @copyright 2008-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2010-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 *
 * PaStiX thread binding routines.
 *
 * @version 6.0.3
 * @author Mathieu Faverge
 * @date 2019-11-12
 *
 */
#include "common.h"
#include "isched_hwloc.h"

#if defined(HAVE_HWLOC)

static hwloc_topology_t topology;
static int first_init = 0;
static int initialized = 0;
static volatile pastix_atomic_lock_t topo_lock = PASTIX_ATOMIC_UNLOCKED;

#if defined(HAVE_HWLOC_PARENT_MEMBER)
#define HWLOC_GET_PARENT(OBJ)  (OBJ)->parent
#else
#define HWLOC_GET_PARENT(OBJ)  (OBJ)->father
#endif  /* defined(HAVE_HWLOC_PARENT_MEMBER) */

#if !defined(HAVE_HWLOC_BITMAP)
#define hwloc_bitmap_t        hwloc_cpuset_t
#define hwloc_bitmap_alloc    hwloc_cpuset_alloc
#define hwloc_bitmap_free     hwloc_cpuset_free
#define hwloc_bitmap_dup      hwloc_cpuset_dup
#define hwloc_bitmap_singlify hwloc_cpuset_singlify
#define hwloc_bitmap_free     hwloc_cpuset_free
#endif

int isched_hwloc_init(void)
{
    int rc = 0;
    pastix_atomic_lock( &topo_lock );
    if ( first_init == 0 ) {
        hwloc_bitmap_t cpuset = hwloc_bitmap_alloc();

        unsigned version = hwloc_get_api_version();
        if ((version >> 16) != (HWLOC_API_VERSION >> 16)) {
            fprintf(stderr,
                    "isched_hwloc_init: PaStiX is compiled for hwloc API 0x%x but running on incompatible library API 0x%x.\n",
                    HWLOC_API_VERSION, version );
            exit(EXIT_FAILURE);
        }

        rc = hwloc_topology_init( &topology );
        if ( rc != 0 ) {
            fprintf(stderr,
                    "isched_hwloc_init: Failed to initialize HwLoc topology. Binding will not be available\n");
            first_init++;
            pastix_atomic_unlock( &topo_lock );
            return -1;
        }

        rc = hwloc_topology_load( topology );
        if ( rc != 0 ) {
            fprintf(stderr,
                    "isched_hwloc_init: Failed to load the HwLoc topology. Binding will not be available\n");
            first_init++;
            pastix_atomic_unlock( &topo_lock );
            return -1;
        }

        rc = hwloc_get_cpubind( topology, cpuset, HWLOC_CPUBIND_PROCESS );
        if ( rc == 0 ) {
#if HWLOC_API_VERSION >= 0x20000
            rc = hwloc_topology_restrict( topology, cpuset, HWLOC_RESTRICT_FLAG_REMOVE_CPULESS );
#else
            rc = hwloc_topology_restrict( topology, cpuset, 0 );
#endif
            if ( rc != 0 ) {
                fprintf(stderr,
                        "isched_hwloc_init: Failed to restrict the topology to the correct cpuset\n"
                        "                   This may generate incorrect bindings\n");
            }
        }
        hwloc_bitmap_free(cpuset);
    }

    initialized = 1;
    first_init++;
    pastix_atomic_unlock( &topo_lock );
    return rc;
}

int isched_hwloc_destroy(void)
{
    pastix_atomic_lock( &topo_lock );
    first_init--;
    if ( (first_init == 0) && initialized ) {
        hwloc_topology_destroy(topology);
    }
    pastix_atomic_unlock( &topo_lock );
    return 0;
}

unsigned int isched_hwloc_nb_cores_per_obj( hwloc_obj_type_t type, int index )
{
    hwloc_obj_t obj = hwloc_get_obj_by_type(topology, type, index);
    assert( obj != NULL );
    return hwloc_get_nbobjs_inside_cpuset_by_type(topology, obj->cpuset, HWLOC_OBJ_CORE);
}

int isched_hwloc_world_size()
{
    return isched_hwloc_nb_cores_per_obj( HWLOC_OBJ_MACHINE, 0 );
}

int isched_hwloc_bind_on_core_index(int cpu_index)
{
    hwloc_obj_t    core;     /* Hwloc object    */
    hwloc_bitmap_t cpuset;   /* Hwloc cpuset    */

    /* Get the core of index cpu_index */
    core = hwloc_get_obj_by_type(topology, HWLOC_OBJ_CORE, cpu_index);
    if (!core) {
        fprintf(stderr,
                "isched_hwloc_bind_on_core_index: unable to get the core of index %i (nb physical cores = %i )\n",
                cpu_index, isched_hwloc_world_size());
        return -1;
    }

    /* Get a copy of its cpuset that we may modify.  */
    cpuset = hwloc_bitmap_dup(core->cpuset);
    hwloc_bitmap_singlify(cpuset);

    /* And try to bind ourself there.  */
    if (hwloc_set_cpubind(topology, cpuset, HWLOC_CPUBIND_THREAD)) {
        char *str = NULL;
        hwloc_bitmap_asprintf(&str, core->cpuset);
        fprintf(stderr, "isched_hwloc: couldn't bind to cpuset %s\n", str);
        free(str);

        /* Free our cpuset copy */
        hwloc_bitmap_free(cpuset);
        return -1;
    }

    /* Get the number at Proc level*/
    cpu_index = core->os_index;

    /* Free our cpuset copy */
    hwloc_bitmap_free(cpuset);
    return cpu_index;
}

int isched_hwloc_unbind()
{
#if defined(HAVE_HWLOC_BITMAP)
    hwloc_obj_t obj;

    /* HwLoc has not been initialized */
    if ( first_init <= 0 ) {
        return -1;
    }

    /* Get last one.  */
    obj = hwloc_get_obj_by_type( topology, HWLOC_OBJ_MACHINE, 0 );
    if (!obj) {
        fprintf(stderr, "isched_hwloc_unbind: Could not get object\n");
        return PASTIX_ERR_UNKNOWN;
    }

    if (hwloc_set_cpubind(topology, obj->cpuset, HWLOC_CPUBIND_THREAD)) {
        char *str = NULL;
        hwloc_bitmap_asprintf(&str, obj->cpuset);
        fprintf(stderr, "isched_hwloc_unbind: Couldn't unbind with cpuset %s\n", str);
        free(str);
        return -1;
    }
#endif
    return PASTIX_SUCCESS;
}

#endif /* defined(HAVE_HWLOC) */
