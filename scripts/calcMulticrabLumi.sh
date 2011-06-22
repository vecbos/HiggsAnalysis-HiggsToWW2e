#!/bin/bash -f
# usage createGoodList.sh multicraboutputdir

crabdir=$1

# do the list of the tasks
grep "\[" $crabdir/multicrab.cfg | grep -v MULTICRAB | awk -F "[" '{print $2}' | awk -F "]" '{print $1}' > tasks.txt

N=0
while read LINE ; do
    N=$((N+1))
    echo "Processing $LINE"
    tasks[${N}]=$LINE
    lumiCalc.py -c frontier://LumiCalc/CMS_LUMI_PROD -i $LINE/res/lumiSummary.json overview >& $LINE.lumi & 
done < tasks.txt

rm -f tasks.txt

