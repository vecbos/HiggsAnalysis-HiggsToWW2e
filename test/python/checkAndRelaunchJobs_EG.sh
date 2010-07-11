#!/bin/sh

###FIXME CHECK THAT 2 STRINGS ARE GIVEN IN INPUT

source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh
eval `scramv1 runtime -sh`
source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh

CASTORDIR=$2
BLACKLIST=`cat ~/physEGAMMA/dashboard/blacklist.txt` 

for run in `find . -type d -name "$1*"| grep -v ".old" | sed -e "s%$1_%%g" | sed -e "s%./%%g"` 
do
    task=${1}_${run}

    echo "++++++++++ Checking task $task ++++++ "

#    rm -rf ${task}/res/lumiSummary.json
    if [ -f ${task}/res/lumiSummary.json ]; then
	echo ">>>> Task ${task} already cleared. Moving on"
	continue
    fi
    
    #CHECKING STATUS. HAVE TO DO IT TWICE SOMETIMES CRAB SERVER IS CRAZY!!
    crab -status -c $task > crab_${task}_status.log 2>&1
    sleep 5
    crab -status -c $task > crab_${task}_status.log 2>&1

    nTotalJobs=`cat crab_${task}_status.log | grep "Total Jobs" | awk '{print $2}'`
    jobsSubmittingCreated=`cat crab_${task}_status.log | egrep "(Submit|Created|Running|Ready)" `
    if [ "${jobsSubmittingCreated}AAA" != "AAA" ]; then
	echo ">>>> ${task} still in submission/created/running/ready state. Better to wait or relaunch it.\n>>>> Do you want to relaunch it now (y/n)? 10 seconds to decide"
	read -t 10 KEYINPUT
	if [ "$KEYINPUT" == "y" ]; then
	    rfrm -r ${CASTORDIR}/${run}; mv -f $task $task.old; crab -create -submit --GRID.ce_black_list='ce02.lcg.cscs.ch,ce11.lcg.cscs.ch' --GRID.se_black_list=${BLACKLIST}  -cfg crab_${task}.cfg > crab_${task}_resubmit.log 2>&1 &
	fi
	continue
    fi

   jobsDoneOK=`cat crab_${task}_status.log | grep -A 2 "Wrapper Exit Code" | grep -A 2 "Wrapper Exit Code : 0" | grep -A2 ">>>>>>>>>" | grep "List of jobs:" | sed -e "s%.*List of jobs: %%g" | awk '{printf "%s,",$1}' | sed -e 's%,$%%g' | awk -f expandJobsList.awk | sed -e 's%,$%%g'` 
    if [ "${jobsDoneOK}AAA" != "AAA" ]; then
	echo ">>>> Jobs ${jobsDoneOK} are OK for task $task"
    fi
   jobsDoneDup=`cat crab_${task}_status.log | grep -A 2 "Wrapper Exit Code" | grep -A 2 "Wrapper Exit Code : 60303" | grep -A2 ">>>>>>>>>" | grep "List of jobs:" | sed -e "s%.*List of jobs: %%g" | awk '{printf "%s,",$1}' | sed -e 's%,$%%g' | awk -f expandJobsList.awk | sed -e 's%,$%%g' ` 
   if [ "${jobsDoneDup}AAA" != "AAA" ]; then
       echo ">>>> Jobs ${jobsDoneDup} are in status 60303 (FileAlready present in the output). In a future version can will be cleaned and relaunched"
   fi

    jobsCrashed=`cat crab_${task}_status.log | grep -A 2 "Wrapper Exit Code" | grep -v "Wrapper Exit Code : 0" | grep -v "Wrapper Exit Code : 60303" | grep -A2 ">>>>>>>>>" | grep "List of jobs:" | sed -e "s%.*List of jobs: %%g" | awk '{printf "%s,",$1}' | sed -e 's%,$%%g'` 
    if [ "${jobsCrashed}AAA" != "AAA" ]; then
	echo ">>>> Relaunching crashed jobs for task $task: ${jobsCrashed}"
	crab -getoutput ${jobsCrashed} -c $task > crab_${task}_getoutput.log ; sleep 10; crab -forceResubmit ${jobsCrashed} --GRID.ce_black_list='ce02.lcg.cscs.ch,ce11.lcg.cscs.ch' --GRID.se_black_list=${BLACKLIST}  -c $task > crab_${task}_resubmit_crashed.log 2>&1 &
    fi

    jobsAborted=`cat crab_${task}_status.log | grep -A 2 "Aborted" |  grep -A2 ">>>>>>>>>" | grep "List of jobs:" | sed -e "s%.*List of jobs: %%g" | awk '{printf "%s,",$1}' | sed -e 's%,$%%g'` 
    if [ "${jobsAborted}AAA" != "AAA" ]; then
	echo ">>>> Relaunching aborted jobs for task $task: jobs ${jobsAborted}"
	crab -forceResubmit ${jobsAborted} --GRID.ce_black_list='ce02.lcg.cscs.ch,ce11.lcg.cscs.ch' --GRID.se_black_list=${BLACKLIST} -c $task > crab_${task}_resubmit_aborted.log 2>&1 &
    fi


    ##FINAL CHECKS. IF EVERYTHING OK THEN CLEARING OUTPUT. CONSIDERING GOOD FOR THE MOMENT ALSO 60303
    jobsDone=${jobsDoneOK},${jobsDoneDup},
    if [ "${jobsDone}AAA" != ",AAA," ]; then
	echo ">>>> Checking output of jobs for task $task"
	FILENOTOK=""
	for job in `echo ${jobsDone} | sed -e 's%,% %g'`
	do
	  FILEOK=`rfdir ${CASTORDIR}/${run} | grep vecbos_${run}_${job} | awk '{print $9}'`
	  if [ "${FILEOK}AAA" == "AAA" ]; then
	      FILENOTOK=`echo "${FILENOTOK}${job},"`
	  fi
	done
	if [ "${FILENOTOK}AAA" != "AAA" ]; then
	    echo ">>>> Jobs ${FILENOTOK} do not have their output. Relaunching them"
	    jobsCrashed=`echo ${FILENOTOK} | sed -e 's%,$%%g'` 
	    crab -getoutput ${jobsCrashed} -c $task > crab_${task}_getoutput.log ; sleep 10; crab -forceResubmit ${jobsCrashed} --GRID.ce_black_list='ce02.lcg.cscs.ch,ce11.lcg.cscs.ch' --GRID.se_black_list=${BLACKLIST}  -c $task > crab_${task}_resubmit_crashed.log 2>&1 &
	else
	    nJobsDoneOK=`echo $jobsDoneOK | awk -F "," '{print NF}'`
	    nJobsDoneDup=`echo $jobsDoneDup | awk -F "," '{print NF}'`
	    (( nJobsDone = nJobsDoneOK + nJobsDoneDup ))
	    if [ "${jobsCrashed}AAA" == "AAA"  -a "${jobsAborted}AAA" == "AAA" -a "${nTotalJobs}AAA" == "${nJobsDone}AAA" ]; then
		echo ">>>> Everything OK. Clearing task ${task}"
		crab -getoutput -c $task > crab_${task}_getoutput.log ; sleep 2; crab -report -c $task > crab_${task}_report.log 2>&1 &
	    else
		echo ">>>> Still not all jobs done/retrieved/ok for ${task}. Try later"
	    fi
	    
	fi
    fi
done
