/**
 *
 * @file bvec.h
 *
 * @copyright 2004-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.1.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2019-11-12
 *
 * @addtogroup bcsc
 * @{
 *
 **/
#ifndef _bvec_h_
#define _bvec_h_

#include <stdlib.h>

void *bvec_malloc( size_t size );

void  bvec_free( void *x );

#endif /* _bvec_h_ */
/**
 * @}
 */
