#! /bin/bash

# run the SV pipeline.

# Usually needed when an integrated run hangs and is killed.
# Run directly in the project working directory. It expects that the files are not "organized" and are at the top level.

# $1 - job name


workingDir=`pwd`

# make this more generic - right now it is dependent on what is organized by organize-impact-files.sh
mv FinalBams/* $workingDir
rmdir FinalBams

mv Results/* $workingDir
rmdir Results

mv AllSampleResults/* $workingDir
rmdir AllSampleResults

mv StdLogFiles/* $workingDir
rmdir StdLogFiles

# define the needed variables
SVConfigFile=/home/patelju1/resources/impact-pipeline/impact-sv.conf
perl=/opt/common/CentOS_6/perl/perl-5.20.2/bin/perl
SVScript=/home/patelju1/software/IMPACT-SV/bin/RunStructuralVariantPipeline_Delly.pl


# most of these things are required to be in the path
export PATH=/opt/common/CentOS_6/R/R-3.1.2/bin:/opt/common/CentOS_6/bedtools/bedtools-2.17.0:/opt/common/CentOS_6/java/jdk1.8.0_31/bin:/opt/common/CentOS_6/samtools/samtools-1.2:$PATH
export TMPDIR=/ifs/work/scratch

jobName=$1.RunSV


# run SV directly
commandString="bsub -q sol -cwd $workingDir -J $jobName -o $jobName.stdout -e $jobName.stderr -We 24:00 -R \"rusage[mem=16]\" -R \"rusage[iounits=0]\" -M 20 $perl $SVScript -c $SVConfigFile -d $workingDir -o $workingDir"

echo $commandString
eval $commandString
