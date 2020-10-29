#
# @package Tools
#
# @copyright (c) 2009-2018 The University of Tennessee and The
#                          University of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                          Univ. Bordeaux. All rights reserved.
#
# @project MORSE
#
# MORSE is a software package provided by:
#      Inria Bordeaux - Sud-Ouest,
#      Univ. of Tennessee,
#      King Abdullah Univesity of Science and Technology
#      Univ. of California Berkeley,
#      Univ. of Colorado Denver.
#
# @author Mark Gates
# @author Mathieu Faverge
#
# @date 2019-4-2
#
# This python script holds the main dictionnary for the substitutions
# of the precision generation tool codegen.py.
#
# Substitutions used in codegen.py
#
# Substitutions are applied in the order listed. This is important in cases
# where multiple substitutions could match, or when one substitution matches
# the result of a previous substitution. For example, these rules are correct
# in this order:
#
#    ('real',   'double precision',  'real',   'double precision' ),  # before double
#    ('float',  'double',            'float',  'double'           ),
#
# but if switched would translate 'double precision' -> 'float precision',
# which is wrong.
#
from __future__ import print_function

import re
import sys, traceback
from os import path
from subs_blas import blas_mixed, blas, lapack

# ===========================================================================
# utility functions

# ----------------------------------------
def _upper( table ):
    '''
    Maps double-nested list of strings to upper case.
    [ ['Foo', 'bar'], ['baz', 'ZAB'] ]
    becomes
    [ ['FOO', 'BAR'], ['BAZ', 'ZAB'] ]
    '''
    ucase = [ [ x.upper() for x in row ] for row in table ]
    return ucase
# end

# ----------------------------------------
def _lower( table ):
    '''
    Maps double-nested list of strings to lower case.
    [ ['Foo', 'BAR'], ['BAZ', 'zab'] ]
    becomes
    [ ['foo', 'bar'], ['baz', 'zab'] ]
    '''
    lcase = [ [ x.lower() for x in row ] for row in table ]
    return lcase
# end

# ----------------------------------------
def _title( table ):
    '''
    Maps double-nested list of strings to Title case. Useful for cuBLAS.
    [ ['FOO', 'bar'], ['Baz', 'Zab'] ]
    becomes
    [ ['Foo', 'Bar'], ['Baz', 'Zab'] ]
    '''
    tcase = [ [ x.title() for x in row ] for row in table ]
    return tcase
# end

# Dictionary is keyed on substitution type (mixed, normal, etc.)
_subs = {
    # ------------------------------------------------------------
    # replacements applied to mixed precision files.
    'mixed' : [
        # double/single,          double/single-complex
        #'12345678901234567890', '12345678901234567890')

        # ----- Special line indicating column types
        ['ds',                   'zc'                  ],

        # ----- Text
        (r'\ba symmetric',      r'\ban hermitian'      ),
        ('symmetric',            'hermitian'           ),
        ('Symmetric',            'Hermitian'           ),
        ('SYMMETRIC',            'HERMITIAN'           ),
        ('orthogonal',           'unitary'             ),

        # ----- CBLAS
        ('',                     'CBLAS_SADDR'         ),
        ('saxpy',                'caxpy'               ),

        # ----- LAPACKE
        ('slange',               'clange'              ),
        ('slarnv',               'clarnv'              ),

        # ----- Complex numbers
        # See note in "normal" section below about regexps
        (r'',                   r'\bconj\b'            ),

        # ----- Constants
        ('Trans',                'ConjTrans'           ),

        # ----- BLAS & LAPACK
    ]
    + _title( blas_mixed )  # e.g., Dgemm, as in cuBLAS, before lowercase (e.g., for Zdrot)
    + _lower( blas_mixed )  # e.g., dgemm
    + _upper( blas_mixed )  # e.g., DGEMM
    + [

        # ----- Fortran Types
        ('real\(kind=c_double\)', 'complex\(kind=c_double_complex\)' ),
        ('real\(kind=c_float\)',  'real\(kind=c_float_complex\)'     ),

        # ----- Data types
        ('double',               'double2'             ),
        ('float',                'float2'              ),
        ('double',               'cuDoubleComplex'     ),
        ('float',                'cuFloatComplex'      ),
        ('DOUBLE PRECISION',    r'COMPLEX\*16'         ),
        ('DOUBLE PRECISION',     'COMPLEX_16'          ),
        ('SINGLE PRECISION',     'COMPLEX'             ),
        ('RealFloat',            'ComplexFloat'        ),
        ('RealDouble',           'ComplexDouble'       ),

        # ----- header files
        (r'_ds\.h\b',           r'_zc\.h\b'            ),
        (r'_ds_h_\b',           r'_zc_h_\b'            ),
        (r'_d\.h\b',            r'_z\.h\b'             ),
        (r'_s\.h\b',            r'_c\.h\b'             ),

        # ----- Prefixes
        (r'\bds_',             r'\bzc_'                ),
        (r'\bd_',              r'\bz_'                 ),
        (r'\bs_',              r'\bc_'                 ),
        (r'\bDS_',             r'\bZC_'                ),
        (r'\bD_',              r'\bZ_'                 ),
        (r'\bS_',              r'\bC_'                 ),

        # See note in "normal" section below
        #('LAPACKE_d',            'LAPACKE_z'           ),
        #('LAPACKE_s',            'LAPACKE_c',          ),
        #('plasma_d',             'plasma_z'            ),
        #('plasma_s',             'plasma_c'            ),

        # ----- Fortran examples
        ('real\(',               'complex\(',          ),
        ('\(transpose\(',        'conjg\(transpose\('  ),
    ],  # end mixed

    # ------------------------------------------------------------
    # replacements applied to most files.
    'normal' : [
        # pattern                single                  double                  single-complex          double-complex
        #'12345678901234567890', '12345678901234567890', '12345678901234567890', '12345678901234567890', '12345678901234567890')

        # ----- header (identifies precision; not a substitution)
        ['p',                    's',                    'd',                    'c',                    'z'                   ],

        # ----- Text
        (r'\ba symmetric',      r'\ba symmetric',      r'\ba symmetric',        r'\ban hermitian',      r'\ba hermitian'       ),
        ('symmetric',            'symmetric',            'symmetric',            'hermitian',            'hermitian'           ),
        ('Symmetric',            'Symmetric',            'Symmetric',            'Hermitian',            'Hermitian'           ),
        ('SYMMETRIC',            'SYMMETRIC',            'SYMMETRIC',            'HERMITIAN',            'HERMITIAN'           ),
        ('orthogonal',           'orthogonal',           'orthogonal',           'unitary',              'unitary'             ),
        ('\^T',                  '\^T',                  '\^T',                  '\^H',                  '\^H'                 ),

        # ----- CBLAS
        ('',                     '',                     '',                     'CBLAS_SADDR',          'CBLAS_SADDR'         ),

        # ----- Complex numbers
        # \b regexp here avoids conjugate -> conjfugate, and fabs -> fabsf -> fabsff.
        # Note r for raw string literals, otherwise \b is a bell character.
        # The \b is deleted from replacement strings.
        # conj() and fabs() are overloaded in MAGMA, so don't need substitution.
        (r'',                   r'',                    r'',                    r'\bconjf\b',           r'\bconj\b'            ),
        (r'',                   r'\bfabsf\b',           r'\bfabs\b',            r'\bfabsf\b',           r'\bfabs\b'            ),
        (r'',                   r'\bfabsf\b',           r'\bfabs\b',            r'\bcabsf\b',           r'\bcabs\b'            ),
        (r'',                   r'\bsqrtf\b',           r'\bsqrt\b',            r'\bcsqrtf\b',          r'\bcsqrt\b'            ),
        (r'',                   r'\bsqrtf\b',           r'\bsqrt\b',            r'\bsqrtf\b',           r'\bsqrt\b'            ),
        (r'',                   r'',                    r'',                    r'\bcrealf\b',          r'\bcreal\b'           ),
        (r'',                   r'',                    r'',                    r'\bcimagf\b',          r'\bcimag\b'           ),

        # ----- Constants
        ('Trans',                'Trans',                'Trans',                'ConjTrans',            'ConjTrans'           ),

        # ----- BLAS & LAPACK
    ]
    + _title( blas )    # e.g., Dgemm, as in cuBLAS, before lowercase (e.g., for Zdrot)
    + _lower( blas )    # e.g., dgemm
    + _upper( blas )    # e.g., DGEMM
    + _lower( lapack )  # e.g., dgetrf
    + _title( lapack )  # e.g., Dgetrf
    + _upper( lapack )  # e.g., DGETRF
    + [
        # ----- Fortran Types
        ('int\(',                'real\(',               'real\(',               'complex\(',            'complex\('           ),
        ('int',                  'real',                 'double precision',     'real',                r'\bdouble precision\b'),  # before double
        ('_int',                 '_float',               '_double',              '_complex',             '_double_complex'     ),
        (r'\bint',              r'\bfloat',             r'\bdouble',            r'\bcomplex',           r'\bdouble_complex'    ),
        (r'\bint',              r'\bfloat',             r'\bdouble',            r'\bfloat',             r'\bdouble'            ),
        ('int',                  'real',                 'double precision',     'complex',             r'\bcomplex\*16'       ),
        ('INT',                  'REAL',                 'DOUBLE_PRECISION',     'COMPLEX',             r'\bCOMPLEX_16'        ),
        ('INT',                  'REAL',                 'DOUBLE PRECISION',     'COMPLEX',             r'\bDOUBLE COMPLEX'    ),
        ('INT',                  'REAL',                 'DOUBLE PRECISION',     'REAL',                r'\bDOUBLE PRECISION'  ),
        ('MPI_INT',              'MPI_FLOAT',            'MPI_DOUBLE',           'MPI_COMPLEX',          'MPI_DOUBLE_COMPLEX'  ),
        ('MPI_INT',              'MPI_FLOAT',            'MPI_DOUBLE',           'MPI_FLOAT',            'MPI_DOUBLE'          ),

        # ----- Data types
        # C++
        ('int',                  'float',                'double',               'float _Complex',      r'\bdouble _Complex'   ),
        # CUDA
        ('int',                  'float',                'double',               'make_cuFloatComplex',  'make_cuDoubleComplex'),
        (r'\bint',              r'\bfloat',             r'\bdouble',            r'\bcuFloatComplex',    r'\bcuDoubleComplex'   ),
        (r'\bint',              r'\breal',              r'\breal',              r'\bcomplex',           r'\bcomplex'           ),
        (r'\bint',              r'\bfloat',             r'\bdouble',            r'\bfloat2',            r'\bdouble2'           ),

        # ----- header files
        (r'_p\.h\b',            r'_s\.h\b',            r'_d\.h\b',             r'_c\.h\b',             r'_z\.h\b'             ),
        ('_p_',                  '_s_',                 '_d_',                  '_c_',                  '_z_'                 ),
        ('_P_',                  '_S_',                 '_D_',                  '_C_',                  '_Z_'                 ),
        (r'\bp_',               r'\bs_',               r'\bd_',                r'\bc_',                r'\bz_'                ),
    ],  # end normal
} #end _subs

class Substitution( object ):
    def __init__( self, subsfiles=[] ):
        # Fill in subs_search with same structure as subs, but containing None values.
        # Later in substitute(), we'll cache compiled regexps in subs_search.
        # We could pre-compile them here, but that would compile many unneeded ones.
        #
        for file in subsfiles:
            fullpath = path.realpath(file)
            (filepath, filename) = path.split( fullpath )

            if filename != "local_subs.py":
                print( "Error in dictionnary: Must be named local_subs.py and not", filename, file=sys.stderr )
                exit(1)

            remove = False
            if filepath not in sys.path:
                sys.path.append( filepath )
                remove = True
            try:
                from local_subs import subs
                for key in subs.keys():
                    _subs[key] = _subs[key] + subs[key]
            except Exception as err:
                print( "Error: in importing:", file, file=sys.stderr )
                traceback.print_exc()
                exit(1)
            if remove:
                sys.path.remove( filepath )

        # Fill in subs_replace with pre-processed version of subs, removing regexp escapes.
        try:
            subs = _subs
            subs_search  = {}
            subs_replace = {}
            for key in subs.keys():
                nrow = len( subs[key]    )
                ncol = len( subs[key][0] )
                subs_search [key] = [ [ None for j in range(ncol) ] for i in range(nrow) ]
                subs_replace[key] = [ [ None for j in range(ncol) ] for i in range(nrow) ]
                for (i, row) in enumerate( subs[key] ):
                    for (j, sub) in enumerate( row ):
                        sub = sub.replace( r'\b',  r''  )
                        sub = sub.replace( r'\*',  r'*' )
                        sub = sub.replace( r'\(',  r'(' )
                        sub = sub.replace( r'\)',  r')' )
                        sub = sub.replace( r'\.',  r'.' )
                        sub = sub.replace( r'\^',  r'^' )
                        subs_replace[key][i][j] = sub
                    # end
                # end
            # end
        except Exception as err:
            print( "Error: in subs:", file=sys.stderr )
            if 'key' in locals() and 'i' in locals():
                print( "row %d of substitution table '%s': %s" %
                       (i, key, row), file=sys.stderr )
            traceback.print_exc()
            exit(1)

        self.subs = subs
        self.subs_search  = subs_search
        self.subs_replace = subs_replace
