#! /bin/bash

# run it from the directory that you would like to be your outDir


# $1 - the project diretory. This directory will be traversed recursively to look for fastq.gz files and corresponding sample sheets
# $2 - the title file. Must not be in the current directory (ie the working directory for this script)
# $3 - project/job name
# $4 - panel type: impact or hemepact. default is impact



workingDir=`pwd`

# choose the right file for the given panel type and copy it to working dir
case $4 in
impact|"") 
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
sed  -i 's/1,2,3,4,5,6,7/1,2,3,4/g' $configFile

# call the low level script
~/software/bin/bsub-innovation-pipeline.sh $1 $2 $3 no $configFile $workingDir

