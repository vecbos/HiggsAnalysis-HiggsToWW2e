#!/bin/sh

source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh
eval `scramv1 runtime -sh`
source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh

BLACKLIST=`cat ~/physEGAMMA/dashboard/blacklist.txt` 
for run in `find . -type d -name "$1*"| sed -e "s%$1_%%g" | sed -e "s%./%%g"` 
do
    task=${1}_${run}
    echo "checking task $task"

    jobsSubmitting=`crab -status -c $task | grep "Submit" `
    if [ "${jobsSubmitting}AAA" != "AAA" ]; then
	echo "problems in submission for ${task}. Need to be relaunched"
	rm -rf $task; crab -create -submit --GRID.ce_black_list='ce02.lcg.cscs.ch,ce11.lcg.cscs.ch' --GRID.se_black_list=${BLACKLIST}  -cfg crab_${task}.cfg > crab_${task}_resubmit.log 2>&1 &
    fi

   jobsDone=`crab -status -c $task | grep -A 2 "Wrapper Exit Code" | grep -A 2 "Wrapper Exit Code : 0" | grep -A2 ">>>>>>>>>" | grep "List of jobs:" | sed -e "s%.*List of jobs: %%g" | awk '{printf "%s,",$1}' | sed -e 's%,$%%g'` 
    if [ "${jobsDone}AAA" != "AAA" ]; then
	echo "Jobs are done. Need to check output of jobs ${jobsDone} for task $task"
    fi

    jobsDone=`crab -status -c $task | grep -A 2 "Wrapper Exit Code" | grep -A 2 "Wrapper Exit Code : 60303" | grep -A2 ">>>>>>>>>" | grep "List of jobs:" | sed -e "s%.*List of jobs: %%g" | awk '{printf "%s,",$1}' | sed -e 's%,$%%g'` 
    if [ "${jobsDone}AAA" != "AAA" ]; then
	echo "Jobs are in status 60303 (FileAlready present in the output). Need to check output of jobs ${jobsDone} for task $task"
    fi

    jobsCrashed=`crab -status -c $task | grep -A 2 "Wrapper Exit Code" | grep -v "Wrapper Exit Code : 0" | grep -v "Wrapper Exit Code : 60303" | grep -A2 ">>>>>>>>>" | grep "List of jobs:" | sed -e "s%.*List of jobs: %%g" | awk '{printf "%s,",$1}' | sed -e 's%,$%%g'` 
    if [ "${jobsCrashed}AAA" != "AAA" ]; then
	echo "relaunching crashed jobs for task $task: ${jobsCrashed}"
	crab -getoutput ${jobsCrashed} -c $task > crab_${task}_getoutput.log ; sleep 10; crab -forceResubmit ${jobsCrashed} --GRID.ce_black_list='ce02.lcg.cscs.ch,ce11.lcg.cscs.ch' --GRID.se_black_list=${BLACKLIST}  -c $task > crab_${task}_resubmit_crashed.log 2>&1 &
    fi

    jobsAborted=`crab -status -c $task | grep -A 2 "Aborted" |  grep -A2 ">>>>>>>>>" | grep "List of jobs:" | sed -e "s%.*List of jobs: %%g" | awk '{printf "%s,",$1}' | sed -e 's%,$%%g'` 
    if [ "${jobsAborted}AAA" != "AAA" ]; then
	echo "relaunching aborted jobs for task $task: jobs ${jobsAborted}"
	crab -forceResubmit ${jobsAborted} --GRID.ce_black_list='ce02.lcg.cscs.ch,ce11.lcg.cscs.ch' --GRID.se_black_list=${BLACKLIST} -c $task > crab_${task}_resubmit_aborted.log 2>&1 &
    fi

done
