#!/bin/bash -f
# usage createGoodList.sh castordir outputdir

crabdir=$1
castordir=$2

rm -f list.txt
grep -R default_data.root $crabdir/res/ | awk -F"-> " '{print $2}' | awk -F "/" '{print "'"$castordir"'" "/" $NF}' | uniq | sort > list.txt
