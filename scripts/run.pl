#!/usr/bin/perl

#--------------------------------------------------------------
# History     30/Apr/2007  First shot
# -------
#             
#             01/May/2007  Added wait option and randomization
#             05/June/2007 Added max Datasets option
#--------------------------------------------------------------
# Emanuele Di Marco  <emanuele@cern.ch>
#--------------------------------------------------------------

print "Starting...\n";

use Time::Local;
use Getopt::Std;
getopts('d:c:b:q:s:C:g:w:l:z:nrh');

if(!$opt_d) {help();}
if(!$opt_c) {help();}
$datasetsFile = $opt_d;
$templateCfg = $opt_c;
if($opt_b){$baseName = $opt_b;}
else{$baseName = "analysis";}
if($opt_q){$queue = $opt_q;}
else{$queue="1nh";}
if($opt_s){$scratch=$opt_s;}
else{$scratch="./"}
if($opt_g){$group=$opt_g;}
else{$group=1;}
if($opt_w){$workdir=$opt_w;}
else{$workdir="/afs/cern.ch/user/e/emanuele/work/HtoWWAnalysis/src/";}
if($opt_z){$wait=$opt_z;}
else{$wait=9999999;}
if($opt_C){$castorDir=$opt_C;}
if($opt_l){$limitToDataset=$opt_l;}
else{$limitToDataset=9999999;}

my $logdir = "$scratch/log";
my $rootdir = "$scratch/output";
my $confdir = "$scratch/conf";
my $scriptdir = "$scratch/script";

# -- Create directories if not yet existent
if (-d "$logdir") {
    # do nothing
} else {
    system("/bin/mkdir $logdir"); 
    system("chmod 755 $logdir");
    if (-d "$logdir") {
	print " -> created $logdir\n";
    } else {
	die "run: cannot create $logdir\n";
    }
}
if (-d "$rootdir") {
    # do nothing
} else {
    system("/bin/mkdir $rootdir"); 
    system("chmod 755 $rootdir");
    if (-d "$rootdir") {
	print " -> created $rootdir\n";
    } else {
	die "run: cannot create $rootdir\n";
    }
}
if (-d "$confdir") {
    # do nothing
} else {
    system("/bin/mkdir $confdir"); 
    system("chmod 755 $confdir");
    if (-d "$confdir") {
	print " -> created $confdir\n";
    } else {
	die "run: cannot create $confdir\n";
    }
}
if (-d "$scriptdir") {
    # do nothing
} else {
    system("/bin/mkdir $scriptdir"); 
    system("chmod 755 $scriptdir");
    if (-d "$scriptdir") {
	print " -> created $scriptdir\n";
    } else {
	die "run: cannot create $scriptdir\n";
    }
}

$logdir = "$logdir";
$confdir = "$confdir";
$scriptdir = "$scriptdir";
if(!$opt_C){$rootdir = "$rootdir";}
else {$rootdir = "rfio:".$castorDir;}

# get some useful environment variables
my $user = $ENV{USER};
my $hostname = $ENV{HOST};
my $pwd = $ENV{PWD};


# open the datasets file 
# the plain rootfiles list from DBS page:
# http://cmsdbs.cern.ch/discovery/
open(DATASETS,"$datasetsFile");
@datasets=<DATASETS>;
$i+=0;
$nfiles=$#datasets;
print "your database contains $nfiles collections\n";


$cfgIndex+=0;
# loop on datasets
for($i=0; $i<($#datasets+1) && $i<$limitToDataset; $i+=$group) {
    $cfgIndex++;
    # open the template cfg file
    open (TPLCFGFILE,"$templateCfg");
    @tplcfgfile=<TPLCFGFILE>;

    # write the cfg file
    $cfgfile = $confdir."/".$baseName."\-$cfgIndex".".cfg";
    open (CFGFILE,">$cfgfile");
    # some standard stuff
    my $now = localtime time;
    print CFGFILE "##\n## This file was generated automatically on $now\n";
    print CFGFILE "## by user $user on host $hostname from $pwd\n";
    print CFGFILE "##\n\n";
    $j+=0;
    for($j=0; $j<($#tplcfgfile+1); $j++) {
	if($tplcfgfile[$j] =~ /\s+\"\/store\/\S+\"/) {
	    $g+=0;
	    for($g=0;$g<$group-1;$g++){
		$dataset=$datasets[$i+$g];
		chop $dataset;
		print CFGFILE "\"$dataset\",\n";
	    }
	    $dataset=$datasets[$i+$group-1];
	    chop $dataset;
	    print CFGFILE "\"$dataset\"\n";
	    print "Writing cfg file $cfgfile... \n";
	}
	elsif($tplcfgfile[$j] =~ /\s+replace\streeDumper.nameFile\s\=\s\"default\.root\"/) {
	    $nameRoot=$rootdir."/".$baseName."\-$cfgIndex".".root";
	    print CFGFILE "replace treeDumper.nameFile = \"$nameRoot\"\n"
	} 
	else {
	    $line=$tplcfgfile[$j];
	    print CFGFILE "$line";
	}
    }
}

# now loop on created cfgfiles, and submit the jobs
if(!$opt_n) {
    @LIST=qx(/bin/ls -1 "$confdir" | /bin/grep "$baseName");
    
    my $subJobs = 0;
    my $totalJobs = $#LIST;
    my $totSubJobs = 0;
    foreach $cfgfile(@LIST) {
	# -- Cut off trailing extension and possible directories in front
	($barefile = $cfgfile) =~ s/\.cfg//g;
	$rest = substr($barefile, 0, rindex($barefile, '/')+1); 
	$barefile =~ s/$rest//;
	chop $barefile;
	$scriptfile = $scriptdir."/".$barefile.".csh";
	open(SCRIPTFILE,">$scriptfile");
	print SCRIPTFILE "#!/bin/csh\n\n";
	print SCRIPTFILE "setenv workdir \"$workdir\"\n";
	print SCRIPTFILE "cd \$workdir\n\n";
	print SCRIPTFILE "eval \`scramv1 ru \-csh\`\n";
	print SCRIPTFILE "project CMSSW\n";
	print SCRIPTFILE "cd \-\n";
	print SCRIPTFILE "cmsRun $confdir\/$cfgfile\n";
	system("/bin/chmod 777 $scriptfile");
	$logfile = $logdir."/".$barefile.".log";

	# at most $wait jobs should be pending, so before submitting more jobs, check that enough have been processed
      retry:
	if($subJobs>$wait) {
	    my @jobs = `/usr/bin/bjobs -u $user -q $queue | /bin/grep PEND`;
	    my $pendJobs = $#jobs;
	    if($pendJobs>$wait){
		my $dateAct = `/bin/date`; chop $dateAct;
		print "total: $totalJobs, submitted: $totSubJobs, limit: $wait, pending: $pendJobs ($dateAct)\n";
		sleep 600;
		goto retry;
	    } 
	    else {
		$subJobs = 0;
	    }
	}
	$totSubJobs++;
	$subJobs++;
	if($opt_r) {
	    system("/usr/bin/bsub -q $queue -J `/afs/cern.ch/user/e/emanuele/public/rline` -C 0 -o $logfile $scriptfile");
	}
	else {
	    system("bsub -q $queue -J $barefile -C 0 -o $logfile $scriptfile");
	}
    }
}


sub help(){
    print "Usage: run.pl -d <datasetsFile> [-l <maxDatasets>] -c <templateCfg> -w <workdir> [-g <group>] [-b <basename>] [-q <queue>] [-s <scratch>] [-C <castorDir>] [-z <maxJobs>] [-n] [-r]\n";
    print "Exemplum: run.pl -d HiggsFileList.txt -l 100 -c hToWWAnalysis.cfg -w ~/work/HtoWWAnalysis/src -g 2 -b qqH160 q 1nh -s ~/scratch0 -z 20 -r \n";
    print "Options:\n";
    print "-d datsetsFile:   choose the txt file with the list of source files as it comes from cms DBS page\n";
    print "-l <maxDatasets>: choose the maximum number of datasets over which to run\n";
    print "-c cfgFile:       choose the template cfg file\n";
    print "-w /afs/.../src:  choose the workdir (the dir where you do eval...)\n";
    print "-g <group>:       choose the grouping: process <group> datasets/job (default is 1)\n";
    print "-b <basename>:    choose the basename for the cfg files and output ROOT files (default is \"analysis\")\n";
    print "-q <queue>:       choose the batch queue to submit (default is \"1nh\")\n";
    print "-s <scratch>:     choose the area where to put the output ROOT files (default is current dir)\n";
    print "-C /castor/...:   choose the CASTOR dir where put the output rootfiles. If not given use <scratch>/output\n";
    print "-z <maxJobs>:     if more than <maxJobs> are pending, wait with submission\n";
    print "-n:               if given, create the cfg files, but not the scriptfiles and do not submit the jobs\n";
    print "-r:               if given, randomize the names of the jobs\n";
    die "enjoy, wasting CPU is the spice of life...\n";
}
