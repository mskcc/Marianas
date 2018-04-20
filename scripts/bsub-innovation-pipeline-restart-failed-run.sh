#! /bin/bash

# restart a failed run.
# The title file, impact.conf and VariantCallingListOfFiles.txt, RecalibrationInputBams.list, MergedPool_title.txt must be chnaged manually before this script is run. (eg. change the steps to run)

# $1 - project/job name


# define the variables needed to run the pipeline
jobName=$1
workingDir=`pwd`
main=~/software/IMPACT-Pipeline/bin/RunIlluminaProcess.pl
configFile=$workingDir/impact.conf
SVConfigFile=~/resources/impact-pipeline/impact-sv.conf
perl=/opt/common/CentOS_6/perl/perl-5.20.2/bin/perl


# most of these things are required to be in the path
export PATH=/opt/common/CentOS_6/R/R-3.1.2/bin:/opt/common/CentOS_6/bedtools/bedtools-2.17.0:/opt/common/CentOS_6/java/jdk1.8.0_31/bin:/opt/common/CentOS_6/samtools/samtools-1.2:$PATH
export TMPDIR=/ifs/work/scratch
export _JAVA_OPTIONS=-Djava.io.tmpdir=/ifs/work/scratch

## Run IMPACT-Pipeline on LSF
commandString="bsub -q sol -cwd $workingDir -J $jobName -eo $jobName.stderr -oo $jobName.stdout -We 96:00 -R \"rusage[mem=16]\" -M 20 $perl $main -c $configFile -d $workingDir -o $workingDir"

echo $commandString
eval $commandString
