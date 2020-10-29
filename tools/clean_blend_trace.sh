#!/usr/bin/env sh
###
#
#  @file clean_blend_trace.sh
#  @copyright 2013-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.3
#  @author Mathieu Faverge
#  @date 2019-11-12
#
###

DIRNAME=$1

FILENAME=$DIRNAME/blend.trace

grep '^3[1-2] .*$' $FILENAME | sed 's/"//g' | sed 's/ /;/g' > /tmp/all_counters.csv

THREADS=$( cut -d';' -f 4 /tmp/all_counters.csv | sort -u | xargs )
for i in $THREADS
do
    grep ";${i};"  /tmp/all_counters.csv > /tmp/counters.csv

    python3 ./clean_counters.py > /tmp/counters_${i}.csv
done

sed '/^3[0-2] .*$/d' $FILENAME > /tmp/header
cat /tmp/header /tmp/counters_T[0-9]*.csv /tmp/counters_Appli.csv > $DIRNAME/blend2.trace

rm -f /tmp/all_counters.csv /tmp/counters.csv
rm -f /tmp/header /tmp/counters_T[0-9]*.csv /tmp/counters_Appli.csv

