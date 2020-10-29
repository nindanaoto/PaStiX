/**
 *
 * @file kernels.c
 *
 * Non precision dependent routines and variables associated to the kernels.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.1.0
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2020-01-29
 *
 **/
#include "common.h"

pthread_mutex_t    pastix_comm_lock = PTHREAD_MUTEX_INITIALIZER;
volatile pthread_t pastix_comm_tid  = (pthread_t)-1;


