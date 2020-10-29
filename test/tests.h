/**
 *
 * @file tests.h
 *
 * Tests functions header.
 *
 * @copyright 2018-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Mathieu Faverge
 * @date 2019-11-12
 *
 **/
#ifndef _tests_h_
#define _tests_h_

#include <stdio.h>
#include "kernels/pastix_lowrank.h"

struct test_param_s {
    int n[3];         /**< Matrix size (min, max, step)            */
    int mode[3];      /**< Matrix generation mode (min, max, step) */
    int prank[3];     /**< Matrix rank percentage (min, max, step) */
    int method[3];    /**< Compression method (min, max, step)     */
    int nb_runs;      /**< Number of run per case                  */
    int use_reltol;   /**< Enable/Disable relative tolerance       */
    double tol_gen;   /**< Tolerance for the matrix generation     */
    double tol_cmp;   /**< Tolerance for the matrix compression    */
    double threshold; /**< Tolerance for the matrix compression    */
    FILE  *output;
};
typedef struct test_param_s test_param_t;

typedef struct test_matrix_s {
    int              m;    /* Number of rows of the test matrix    */
    int              n;    /* Number of columns of the test matrix */
    int              ld;   /* Leading dimension of the test matrix */
    int              rk;   /* Required rank of the test matrix     */
    double           norm; /* Norm of the full rank matrix         */
    void            *fr;   /* Full rank matrix                     */
    pastix_lrblock_t lr;   /* Low rank matrix                      */
} test_matrix_t;

void testGetOptions( int argc, char **argv, test_param_t *params, double eps );

#endif /* _tests_h_ */
