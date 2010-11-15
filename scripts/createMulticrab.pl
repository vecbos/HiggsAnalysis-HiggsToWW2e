#!/usr/bin/perl

use Time::Local;
use Getopt::Std;

getopts('d:m:');
$datasetsFile = $opt_d;
if($opt_m) {$multicrabFile = $opt_m;}
else {$multicrabFile = "multicrab.cfg";}

open(DATASETS,"$datasetsFile");
@datasets=<DATASETS>;
$i+=0;
$ndatasets=$#datasets;
print "your file contains $ndatasets datasets\n";

open (MULTICRABFILE,">$multicrabFile");
print MULTICRABFILE "[MULTICRAB]\n";
print MULTICRABFILE "cfg=crab.cfg\n\n";

# loop on datasets
for($i=0; $i<($#datasets+1); $i++) {
    $dataset=$datasets[$i];
    chop $dataset;
    print "processing $dataset ...";
    if($dataset =~ /\/(\S+)\/\S+\/\S+/) {        
        $dir = $1;
        print MULTICRABFILE "[$dir]\n";
        print MULTICRABFILE "CMSSW.datasetpath=$dataset\n\n";
    }
}

print "Done. multicrab file in: $multicrabFile\n";

