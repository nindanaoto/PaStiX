#!/usr/bin/env python
"""@package Tools

 @copyright (c) 2009-2018 The University of Tennessee and The
                          University of Tennessee Research Foundation.
                          All rights reserved.
 @copyright (c) 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                          Univ. Bordeaux. All rights reserved.

 @project MORSE

MORSE is a software package provided by:
     Inria Bordeaux - Sud-Ouest,
     Univ. of Tennessee,
     King Abdullah Univesity of Science and Technology
     Univ. of California Berkeley,
     Univ. of Colorado Denver.

This python script is responsible for precision generation
replacements as well as replacements of any kind in other
files. Different types of replacements can be defined such that no two
sets can conflict.  This script is a issued from the updated version
in plasma-openmp developed by Mark Gates, and inspired by the one
developed by Alvaro Wes in the original Plasma.  This version is
modified as the original one to fit the need of the cmake
RulesPrecision module of morse package used in packages such as
Chameleon, PaStiX, or DPLASMA.

@author Mark Gates
@author Mathieu Faverge

@date 2019-4-2

"""
from __future__ import print_function
import sys, traceback
import re
from datetime import datetime
from os import path;
from optparse import OptionParser,OptionGroup
from subs import Substitution

__version__=2.0
_verbose=False

# ------------------------------------------------------------
class SourceFile( object ):
    '''SourceFile encapsulates a single file.
    If the file contains @precisions, it is a generator file, otherwise the script does nothing.
    It handles determining output files, cmake rules, and generating other precisions.'''

    # --------------------
    # matches "@precisions normal z -> s d c"
    #                                          ($1_)  ($2_)      ($3_________)
    precisions_re = re.compile( r"@precisions +(\w+) +(\w+) +-> +(\w+( +\w+)*)" )
    generated_re  = re.compile( r"@generated" )

    # --------------------
    def __init__( self, subs, filename, filepath=None, precisions=None ):
        '''Creates a single file.
        Determines whether it can do precision generation and if so,
        its substitution table, source and destination precisions.
        '''
        self._subs = subs
        self._filename = filename
        if filepath:
          self._filepath = filepath
        else:
          self._filepath = ''

        self._fullname = path.realpath( path.join(self._filepath, filename) )
        fd = open( self._fullname, 'r' )
        self._text = fd.read()
        fd.close()
        m = self.precisions_re.search( self._text )
        if m:
            self._is_generator = True
            self._table = m.group(1)          # e.g.:  normal or mixed
            self._src   = m.group(2).lower()  # e.g.:  z
            self._dsts  = set(m.group(3).lower().split())  # e.g.:  s, d, c
            self._dsts.add( self._src )
            if precisions:
              self._dsts.intersection_update( precisions )
        else:
            self._is_generator = False
            self._table = None
            self._src   = None
            self._dsts  = set()
    # end

    # --------------------
    def is_generator( self ):
        '''True if this file can generate other precisions (has @precisions).'''
        return self._is_generator

    # --------------------
    def get_filename( self, precision, output_path=None ):
        '''Returns the output filename for a given precision.'''
        if not self.is_generator():
            return None

        outname = self._substitute( self._filename, precision )

        if output_path:
            basename= list(path.split(outname))[1]
            outname = path.join( output_path, basename );

        # Get the relative path to the output
        outname = path.relpath(outname);

        # If full pathnames are identical, no file is generated
        if path.realpath(outname) != self._fullname:
            return outname
        return None

    # --------------------
    def get_filenames( self, precisions=None, output_path=None ):
        '''Returns (files, precs) for the given precisions.
        files is list of generated filenames,
        precs is list of corresponding precisions.
        If precision is None, returns for all of file's precisions.
        If file is not a generator, returns empty lists.'''
        files = []
        ps    = []
        if self.is_generator():
            for prec in self._get_precisions( precisions ):
                outname = self.get_filename( prec, output_path )
                if outname:
                    files.append( outname )
                    ps   .append( prec )
        return (files, ps)

    # --------------------
    def get_cmake_deps( self, precisions=None, output_path=None ):
        string=""
        (files, precs) = self.get_filenames( precisions, output_path )
        for (outname, prec) in zip( files, precs ):
            string += self._filename + ',' + prec + ',' + outname + ';'
        return string

    # --------------------
    def generate_file( self, precision, output_path=None ):
        '''Generate the file for the given precision.
        If file is not a generator, does nothing.'''
        if not self.is_generator():
            return

        outname = self.get_filename( precision, output_path )
        if outname:
            if _verbose:
                print( "generating", outname, file=sys.stderr )
            output = self._substitute( self._text, precision )
            output = self._substitute_precision( output, precision )
            fd = open( outname, 'w' )
            fd.write( output )
            fd.close()

    # --------------------
    def generate_files( self, precisions=None, output_path=None ):
        '''Generates files for the given precisions.
        If precision is None, generates for all of file's precisions.
        If file is not a generator, does nothing.'''
        for prec in self._get_precisions( precisions ):
            self.generate_file( prec, output_path )

    # --------------------
    def _get_precisions( self, precisions=None ):
        '''Given a precision or list of precisions,
        returns list of those that apply to this file.
        If precision is None, returns all of file's precisions.
        '''
        ps = self._dsts
        if precisions:
            ps.intersection_update( precisions );
        return ps
    # end

    # --------------------
    def _substitute( self, text, precision ):
        '''Apply substitutions to text for given precision.'''
        try:
            # Get substitution table based on self._table
            subs_o = self._subs.subs[         self._table ]  # original
            subs_s = self._subs.subs_search[  self._table ]  # compiled as search regexp
            subs_r = self._subs.subs_replace[ self._table ]  # with regexp removed for replacement

            # Get which column is from and to.
            header = subs_o[0]
            jfrom = header.index( self._src )
            jto   = header.index( precision )
        except Exception as err:
            print( "Error: bad table or precision in '%s', @precisions %s %s -> %s:" %
                   (self._filename, self._table, self._src, self._dsts), file=sys.stderr )
            traceback.print_exc()
            exit(1)

        # Apply substitutions
        try:
            line = 0
            for (orig, search, replace) in zip( subs_o[1:], subs_s[1:], subs_r[1:] ):
                line += 1
                if search[jfrom] is None:
                    search[jfrom] = re.compile( orig[jfrom] )
                text = re.sub( search[jfrom], replace[jto], text )
            # end
        except Exception as err:
            print( "Error: in row %d of substitution table '%s': %s" %
                   (line, self._table, subs_o[line]), file=sys.stderr )
            traceback.print_exc()
            exit(1)

        return text
    # end

    # --------------------
    def _substitute_precision( self, text, precision ):
        '''Subtitute the @precision line'''

        # Replace @precision with @generated, file, rule, and timestamp
        gen = "@generated from %s, %s %s -> %s, %s" % (
            self._filename, self._table, self._src, precision, datetime.now().ctime())
        text = re.sub( self.precisions_re, gen, text )
        return text
    # end
# end SourceFile

def main():
  """Create the options parser for detecting options on the command line."""
  parser = OptionParser(usage="Usage: %prog [options]",version='%prog '+str(__version__));
  group = OptionGroup(parser,"Main Options","Only one of these two options is valid");
  group.add_option("-c", "--cmake",    help='Print the list of dependenencies for the RulesPrecision.cmake package', action='store_true', dest='cmake',    default=False);
  group.add_option("-g", "--genfiles", help='Generate the converted files.',                                         action='store_true', dest='genfiles', default=False);
  parser.add_option_group(group);
  group = OptionGroup(parser,"Settings","These options specify how the work should be done.");
  group.add_option("-s", "--srcdir",      help='The input source directory.',                                 action='store', dest='srcdir',                     default=None);
  group.add_option("-P", "--prefix",      help='The output directory if different from the input directory.', action='store', dest='prefix',                     default=None);
  group.add_option("-f", "--file",        help='Specify a file(s) on which to operate.',                      action='store', dest='fileslst',    type='string', default="");
  group.add_option("-p", "--prec",        help='Specify a precision(s) on which to operate.',                 action='store', dest='precslst',    type='string', default="");
  group.add_option("-e", "--filetypes",   help='Specify file extensions on which to operate when walking.',   action='store', dest='fileexts',    type='string', default="");
  group.add_option("-D", "--dictionnary", help='Specify the dictionnary to use in names conversion.',         action='store', dest='dictionnary', type='string', default="");
  parser.add_option_group(group);

  (options, args) = parser.parse_args();

  subs = Substitution( options.dictionnary.split() )

  if ((( not options.cmake ) and ( not options.genfiles )) or
      ((     options.cmake ) and (     options.genfiles ))):
    parser.usage()
    exit(1)

  cmakestr = ""
  for file in options.fileslst.split():
    """For each valid conversion file found."""
    try:
      """Try creating and executing a converter."""
      srcfile = SourceFile( subs, file, filepath=options.srcdir );
      if options.genfiles:
        srcfile.generate_files( precisions=options.precslst.split(),
                                output_path=options.prefix )
      elif options.cmake:
        cmakestr += srcfile.get_cmake_deps( precisions=options.precslst.split(),
                                            output_path=options.prefix )

    except Exception as err:
      print("Unexpected error (%s):"%(file), sys.exc_info()[0], file=sys.stderr)
      traceback.print_exc()
      exit(1)

  if options.cmake:
    print(cmakestr);

if __name__ == "__main__":
    main();
