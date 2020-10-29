/**
 *
 * @file tests_auxiliary.c
 *
 * Tests option reader.
 *
 * @copyright 2018-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Esragul Korkmaz
 * @author Mathieu Faverge
 * @date 2019-11-12
 *
 **/
#include "common.h"
#include <unistd.h>
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */
#include <string.h>
#include <lapacke.h>
#include "tests.h"

static int
get_range( char *range, int params[3] ) {
    int rc;

    rc = sscanf( range, "%d:%d:%d", params, params+1, params+2 );
    if ( rc != 3 ) {
        rc = sscanf( range, "%d:%d", params, params+1 );
        rc++;
    }
    if ( rc != 3 ) {
        rc = sscanf( range, "%d", params );
        rc+=2;
    }
    if ( rc != 3 ) {
        return -1;
    }
    return 0;
}

/**
 * @brief Print default usage for low-rank test binaries
 */
static inline void
test_usage(void)
{
    fprintf(stderr,
            "Low-rank test parameters:\n"
            " -n --n            : Defines the matrix dimension\n"
            "    --n_range      : Defines the matrix dimension range nmin:nmax:nstep\n"
            " -m --mode         : Defines the matrix mode generation\n"
            "    --m_range      : Defines the matrix mode generation range mmin:mmax:mstep\n"
            " -p --percent      : Defines the rank percentage of the matrix dimension\n"
            "    --p_range      : Defines the rank percentage range pmin:pmax:pstep\n"
            " -x --method       : Defines the compression method\n"
            "    --x_range      : Defines the compression method range xmin:xmax:xstep\n"
            " -i --niter        : Defines the number of iteation per test case\n"
            " -r --reltol       : Defines is the system uses relative or absolute tolerance\n"
            " -t --tol_gen      : Defines the precision with which the matrix is generated\n"
            "    --tol_cmp      : Defines the compression tolerance (Default: tol_gen)\n"
            " -s --threshold    : Defines the precision threshold with which the matrix is generated\n"
            "\n"
            " -o --output       : Defines the output filename (default:stdout)\n"
            "\n"
            " -h --help         : this message\n"
            "\n"
            );
}

/**
 * @brief Define the options and their requirement used by the low-rank tests
 */
#define GETOPT_STRING "i:n:N:m:M:p:P:x:X:r:t:T:o:s:h"

#if defined(HAVE_GETOPT_LONG)
/**
 * @brief Define the long options when getopt_long is available for the low-rank test parameters
 */
static struct option long_options[] =
{
    {"n",         required_argument, 0, 'n'},
    {"n_range",   required_argument, 0, 'N'},
    {"m",         required_argument, 0, 'm'},
    {"m_range",   required_argument, 0, 'M'},
    {"p",         required_argument, 0, 'p'},
    {"p_range",   required_argument, 0, 'P'},
    {"r",         required_argument, 0, 'r'},
    {"reltol",    required_argument, 0, 'r'},
    {"x",         required_argument, 0, 'x'},
    {"x_range",   required_argument, 0, 'X'},
    {"tol_gen",   required_argument, 0, 't'},
    {"tol_cmp",   required_argument, 0, 'T'},
    {"threshold", required_argument, 0, 's'},
    {"output",    required_argument, 0, 'o'},
    {"help",      no_argument,       0, 'h'},
    {"i",         required_argument, 0, 'i'},
    {"niter",     required_argument, 0, 'i'},
    {0, 0, 0, 0}
};
#endif  /* defined(HAVE_GETOPT_LONG) */

void
testGetOptions( int argc, char **argv,
                test_param_t *params,
                double        eps )
{
    int c;

    /*
     * Let's initialize the default parameters
     */
    params->n[0] = 300;
    params->n[1] = 500;
    params->n[2] = 100;

    params->mode[0] = 0;
    params->mode[1] = 2;
    params->mode[2] = 1;

    params->prank[0] = 5;
    params->prank[1] = 20;
    params->prank[2] = 5;

    params->method[0] = PastixCompressMethodSVD;
    params->method[1] = PastixCompressMethodNbr-1;
    params->method[2] = 1;

    params->nb_runs    = 1;
    params->use_reltol = 0;

    params->tol_gen = sqrt( eps );
    params->tol_cmp = params->tol_gen;
    params->threshold = params->tol_gen * params->tol_gen;

    params->output = stdout;

    do
    {
#if defined(HAVE_GETOPT_LONG)
        c = getopt_long( argc, argv, GETOPT_STRING,
                         long_options, NULL );
#else
        c = getopt( argc, argv, GETOPT_STRING );
#endif  /* defined(HAVE_GETOPT_LONG) */

        switch(c)
        {
        case 'i':
            params->nb_runs = atoi( optarg );
            break;

        case 'n':
            params->n[0] = atoi( optarg );
            params->n[1] = params->n[0];
            params->n[2] = 1;
            break;

        case 'N':
            get_range( optarg, params->n );
            break;

        case 'm':
            params->mode[0] = atoi( optarg );
            params->mode[1] = params->mode[0];
            params->mode[2] = 1;
            break;

        case 'M':
            get_range( optarg, params->mode );
            break;

        case 'p':
            params->prank[0] = atoi( optarg );
            params->prank[1] = params->prank[0];
            params->prank[2] = 1;
            break;

        case 'P':
            get_range( optarg, params->prank );
            break;

        case 'r':
            params->use_reltol = atoi( optarg ) ? 1 : 0;
            break;

        case 'x':
            params->method[0] = atoi( optarg );
            params->method[1] = params->method[0];
            params->method[2] = 1;
            break;

        case 'X':
            get_range( optarg, params->method );
            break;

        case 't':
            params->tol_gen = atof( optarg );
            params->tol_cmp = params->tol_gen;
            params->threshold = params->tol_gen * params->tol_gen;
            break;

        case 'T':
            params->tol_cmp = atof( optarg );
            break;

        case 's':
            params->threshold = atof( optarg );
            break;

        case 'o':
            params->output = fopen( optarg, "w" );
            if ( !params->output ) {
                perror( "Failed to open output file" );
                exit(EXIT_FAILURE);
            }
            break;

        case 'h':
            test_usage(); exit(EXIT_FAILURE);

        case ':':
            fprintf(stderr, "\nOption %c is missing an argument\n\n", c );
            goto unknown_option;

        case '?': /* getopt_long already printed an error message. */
            test_usage(); exit(EXIT_FAILURE);
        default:
            break;
        }
    } while(-1 != c);

    return;

  unknown_option:
    test_usage(); exit(EXIT_FAILURE);
}
