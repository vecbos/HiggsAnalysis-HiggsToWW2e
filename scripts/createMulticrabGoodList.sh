#!/bin/bash -f
# usage createGoodList.sh multicraboutputdir castordir

crabdir=$1
castordir=$2

# do the list of the tasks
grep "\[" $crabdir/multicrab.cfg | grep -v MULTICRAB | awk -F "[" '{print $2}' | awk -F "]" '{print $1}' > tasks.txt

N=0
while read LINE ; do
    N=$((N+1))
    echo "Processing $LINE"
    tasks[${N}]=$LINE
    grep -R default_data.root $LINE/res/ | awk -F"-> " '{print $2}' | awk -F "/" '{print "'"$castordir"'" "/" $NF}' | uniq | sort > $LINE.list
done < tasks.txt

rm -f tasks.txt

