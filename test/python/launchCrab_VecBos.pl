#!/usr/bin/perl

use Getopt::Long;

$TASKNAME="MinBias-v9_V2";
$DATASET="/MinimumBias/Commissioning10-GOODCOLL-v9/RAW-RECO";
$DIRNAME="MinBias-v9_V2";
$RUNFILE="runs_vecbos.txt";
$TEMPLATECFG="";
GetOptions(
	   'help|h|?'=>\$help,
	   'dataset=s'=>\$DATASET,
	   'dirname=s'=>\$DIRNAME,
	   'taskname=s'=>\$TASKNAME,
	   'runfile=s'=>\$RUNFILE,
	   'templateCfg=s'=>\$TEMPLATECFG,
	   );

#Get blacklist for Rizzi's script
$BLACKLIST=`cat /afs/cern.ch/user/m/meridian/physEGAMMA/dashboard/blacklist.txt`;
chomp($BLACKLIST);

sub usage() {
    print <<EOF;

  usage: $0 [options] 

    options:
      --dataset DATASET: 
      Dataset to monitor; defaults to '$DATASET'

      --dirname DIRNAME: 
      Output dir in castor; will be created under /castor/cern.ch/user/m/meridian/VecBos/

      --taskname TASKNAME: 
      Reference name for the task

      --runfile RUNFILE: 
      File containing the list of runs to be processed in the format RUN NEVTSXJOB TOTEVTSXJOB

      --templateCfg TEMPLATECFG: 
      template CMSSW cfg to launch. Every line with CHANGERUNMBER will be converted in the actual processed RUN
      
EOF
exit(2);
}

usage() if ($help);
# or scalar(@ARGV) == 0);
open FILE, "${RUNFILE}" or die $!;
system "rfmkdir /castor/cern.ch/user/m/meridian/VecBos/${DIRNAME} && rfchmod 775 /castor/cern.ch/user/m/meridian/VecBos/${DIRNAME}";

while (<FILE>) { 
    next if /^\#/;
    $line=$_;
    chomp($line);
    @fields =  split(/\s+/,$line);
    $RUNNUMBER=$fields[0];
    $DATASETSHORT=${TASKNAME}."_".${RUNNUMBER};
    $EVENTSPERJOB=$fields[1];
    $NEVENTS=$fields[2];
    print "$fields[0], $DATASETSHORT, $fields[1], $fields[2] \n"; 
    system "cat ${TEMPLATECFG} | sed -e \"s%CHANGERUNNUMBER%${RUNNUMBER}%g\" > vecbos_data_${RUNNUMBER}.py";
    system "rfmkdir /castor/cern.ch/user/m/meridian/VecBos/${DIRNAME}/${RUNNUMBER} && rfchmod 775 /castor/cern.ch/user/m/meridian/VecBos/${DIRNAME}/${RUNNUMBER}";
    open CRABCFG, ">crab_$DATASETSHORT.cfg" or die $!;
print CRABCFG 
    <<EOF;
[CMSSW]
lumis_per_job=$EVENTSPERJOB
total_number_of_lumis=$NEVENTS
# number_of_jobs=1
#dbs_url=http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet 
datasetpath=${DATASET}
pset=vecbos_data_${RUNNUMBER}.py
#increment_seeds=generator,g4SimHits
runselection=${RUNNUMBER}
#/MinimumBias/BeamCommissioning09-PromptReco-v4/RECO
#use_parents=1
output_file=vecbos_${RUNNUMBER}.root
#split_by_run=1

[USER]
return_data=0
email=Paolo.Meridiani@cern.ch
ui_working_dir=${DATASETSHORT}
copy_data=1
space_token = CMS_T3
return_data = 0
copy_data = 1

storage_element = srm-cms.cern.ch
storage_path=/srm/managerv2?SFN=/castor/cern.ch
user_remote_dir=/user/m/meridian/VecBos/${DIRNAME}/${RUNNUMBER}
additional_input_files=/afs/cern.ch/user/e/emanuele/public/4Likelihood/PDFsSQLite/CMSSW_3_2_X/electronIdLikelihoodTkIsolated.db
check_user_remote_dir=0

[CRAB]
scheduler=glite
jobtype=cmssw
use_server=1
#server_name=ucsd

[GRID]
ce_black_list=ce02.lcg.cscs.ch,ce11.lcg.cscs.ch
se_black_list=${BLACKLIST}

EOF
system(". /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh; eval `scramv1 runtime -sh`; . /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh; crab -create -submit -cfg crab_${DATASETSHORT}.cfg  > crab_${DATASETSHORT}.log 2>&1 &");
sleep 5
}
