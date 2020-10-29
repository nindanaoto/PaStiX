subs = {
    # ------------------------------------------------------------
    # replacements applied to mixed precision files.
    'normal': [
        # pattern                single                  double                  single-complex          double-complex
        #'12345678901234567890', '12345678901234567890', '12345678901234567890', '12345678901234567890', '12345678901234567890')
        ('int',                  'float',                'double',               'spm_complex32_t',     r'\bspm_complex64_t'   ),
        ('SpmPattern',           'SpmFloat',             'SpmDouble',            'SpmComplex32',        r'\bSpmComplex64'      ),
        ('SpmPattern',           'SpmFloat',             'SpmDouble',            'SpmFloat',            r'\bSpmDouble'         ),
        ('MPI_INT',              'MPI_FLOAT',            'MPI_DOUBLE',           'MPI_COMPLEX32',        'MPI_COMPLEX64'       ),

        # ----- Variables
        (r'\b',                 r'szero\b',             r'dzero\b',             r'czero\b',             r'zzero\b'             ),
        (r'\b',                 r'sone\b',              r'done\b',              r'cone\b',              r'zone\b'              ),

        # ----- SPM Prefixes
        ('spm_p',                'spm_s',                'spm_d',                'spm_c',                'spm_z'               ),
    ]
}
