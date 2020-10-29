#!/usr/bin/env sh
###
#
#  @file filelist.sh
#  @copyright 2013-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @brief Generate the filelist for the static analysis
#
#  @version 6.0.3
#  @author Mathieu Faverge
#  @date 2019-11-12
#
###

if [ $# -gt 0 ]
then
    BUILDDIR=$1
fi
BUILDDIR=${BUILDDIR-=build}

echo $PWD
rm -f filelist.txt

git ls-files | grep "\.[ch]"   >  filelist.txt
git ls-files | grep "\.py"     >> filelist.txt
find $BUILDDIR -name '*\.[ch]' >> filelist.txt
echo "${BUILDDIR}/include/pastix/config.h" >> filelist.txt
echo "wrappers/python/examples/pypastix/enum.py" >> filelist.txt

# Remove files in kernel/gpus that are C++ and not our own files.
sed -i "/kernels\/gpus\/.*/d" filelist.txt

# Remove all CMakeFiles generated file
sed -i '/CMakeFiles/d' filelist.txt

# Remove all .cmake files
sed -i '/.cmake/d' filelist.txt

# Remove all .in files
sed -i '/.in$/d' filelist.txt

# Remove all clang files
sed -i '/^\.clang/d' filelist.txt

# Remove files compiled from jdf files
for jdf in `find -name "*\.jdf"`
do
    name=`basename $jdf .jdf`
    echo Remove $name
    sed -i "/$name\.[ch]/d" filelist.txt
done

# Remove installed files
sed -i '/^install.*/d' filelist.txt

# Remove original files used for precision generation
for file in `git grep "@precisions" | awk -F ":" '{ print $1 }'`
do
    sed -i "\:^$file.*:d" filelist.txt
done

# Remove external header files
for file in include/cblas.h include/lapacke.h
do
    sed -i "\:^$file.*:d" filelist.txt
done

# Remove submodules
for file in spm cmake_modules/morse_cmake
do
    sed -i "\:^$file:d" filelist.txt
done

grep '/\.c$/d' filelist.txt > filelist-c.txt

