#
#  @file update_release.sh
#
#  @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.1.0
#  @author Mathieu Faverge
#  @date 2019-11-12
#
#!/usr/bin/env sh

#
# Steps to update header information before doing the release
#
#   1) First update the date of the files with the following lines
#
for f in `git ls-files`
do
    date=`git log -1 --format=%cd --date=short CMakeLists.txt`
    echo $date $f
    sed -i "s/date [-0-9]*\$/date $date/" $f
done

#
#   2) Update the release number
#
git grep -E "version 6\.0_.[01]" | awk -F ":" '{ print $1 }' | sort -u | xargs sed -i 's/version 6\.0\.[01]/version 6.1.0/'

#
# Or this one to update only changed files since last release
#
for i in $( git diff v6.0.2 --name-only ); do if [ -f $i ]; then sed -i 's/@version [0-9].[0-9].[0-9]/@version 6.0.3/' $i; fi; done

#
#   3) Update manually the version number in CMakeLists.txt
#

#
#   4) If necessary, update the copyright information
#
git grep -E "copyright [0-9]{4}-[0-9]{4}" | awk -F ":" '{ print $1 }' | sort -u | xargs sed -i 's/copyright \([0-9]*\)-[0-9]* Bordeaux/copyright \1-2020 Bordeaux/'

#
#   5) Check that the fortran/python wrappers have been updated (see gen_wrappers.py)
#   6) Check header files with check_headers.sh
#   7) Update homebrew formula (only after release)
#
