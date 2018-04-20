#! /bin/bash

# run it from the directory that you would like to be your outDir


# $1 - the project diretory. This directory will be traversed to look for fastq.gz files and corresponding sample sheets
# $2 - the title file. Must not be in the current directory (ie the working directory for this script)
# $3 - project/job name
# $4 - panel type: impact or hemepact
# $5 - run SV? yes or no.



workingDir=`pwd`

# choose the right file for the given panel type and copy it to working dir
case $4 in
impact) 
	configFile=$workingDir/impact.conf
	cp ~/resources/impact-pipeline/impact.conf $configFile
	;;
hemepact)
	configFile=$workingDir/hemepact.conf
	cp ~/resources/impact-pipeline/hemepact.conf $configFile
	;;
*) 
	echo -e "panel type not recognized: $4"
	exit
	;;
esac

# modify the config file as needed


# call the low level script
~/software/bin/bsub-innovation-pipeline.sh $1 $2 $3 $5 $configFile $workingDir

