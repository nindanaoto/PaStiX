/**
 *
 * @file lowrank.c
 *
 * PaStiX low-rank common structures to store pointer to the multiple functions.
 *
 * @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2019-12-19
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "pastix_zlrcores.h"
#include "pastix_clrcores.h"
#include "pastix_dlrcores.h"
#include "pastix_slrcores.h"

double pastix_lr_minratio = 1.0;
pastix_int_t pastix_lr_ortho = 0;

const char *compmeth_shnames[PastixCompressMethodNbr] = {
    "SVD",
    "PQRCP",
    "RQRCP",
    "TQRCP",
    "RQRRT",
};

const char *compmeth_lgnames[PastixCompressMethodNbr] = {
    "Singular Values Decomposition",
    "Partial QR with Column Pivoting",
    "Randomized QR with Column Pivoting",
    "Truncated QR with Column Pivoting",
    "Randomized QR with QR rotation",
};

const fct_ge2lr_t ge2lrMethods[PastixCompressMethodNbr][4] =
{
    { core_sge2lr_svd,   core_dge2lr_svd,   core_cge2lr_svd,   core_zge2lr_svd   },
    { core_sge2lr_pqrcp, core_dge2lr_pqrcp, core_cge2lr_pqrcp, core_zge2lr_pqrcp },
    { core_sge2lr_rqrcp, core_dge2lr_rqrcp, core_cge2lr_rqrcp, core_zge2lr_rqrcp },
    { core_sge2lr_tqrcp, core_dge2lr_tqrcp, core_cge2lr_tqrcp, core_zge2lr_tqrcp },
    { core_sge2lr_rqrrt, core_dge2lr_rqrrt, core_cge2lr_rqrrt, core_zge2lr_rqrrt },
};

const fct_rradd_t rraddMethods[PastixCompressMethodNbr][4] =
{
    { core_srradd_svd,   core_drradd_svd,   core_crradd_svd,   core_zrradd_svd   },
    { core_srradd_pqrcp, core_drradd_pqrcp, core_crradd_pqrcp, core_zrradd_pqrcp },
    { core_srradd_rqrcp, core_drradd_rqrcp, core_crradd_rqrcp, core_zrradd_rqrcp },
    { core_srradd_tqrcp, core_drradd_tqrcp, core_crradd_tqrcp, core_zrradd_tqrcp },
    { core_srradd_pqrcp, core_drradd_pqrcp, core_crradd_pqrcp, core_zrradd_pqrcp }
};
