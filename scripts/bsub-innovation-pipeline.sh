#! /bin/bash

# the low level script that actually runs the pipeline.
# this script must not be called directly. Must be called from one of the higher level scripts.

# $1 - the project diretory. This directory will be searched for fastq.gz files and sample sheets corresponding to the samples mentioned in the title file
# $2 - the title file. Must not be in the current directory (ie the working directory for this script)
# $3 - project/job name
# $4 - run SV? yes or no
# $5 - full path to config file
# $6 - full path to working dir

# copy the title file ici
cp $2 title_file.txt

resources=~/resources/impact-pipeline
barcodeKeyFile=$resources/barcodeKey96.txt

# just get the header from the first sample's sample sheet
sampleID=`tail -n +2 title_file.txt | head -1 | cut -f 3`
head -1 $1/Sample_$sampleID/SampleSheet.csv > SampleSheet.csv


# make the combined sample sheet
for sampleID in  `tail -n +2 title_file.txt | cut -f 3` #`find $1 -name SampleSheet.csv`
do
	f=$1/Sample_$sampleID/SampleSheet.csv
	tail -n +2 $f | awk 'BEGIN{FS=OFS=","}{gsub(/_/, "-", $3); print}' >> SampleSheet.csv
done


# make symlinks for all the fastq files using IMPACT-valid names
# That includes taking care of - vs _ and replacing the S.. thing with proper barcode.
# assumption 1: the proper sample name itself will never have _.
# assumption 2: the sample id has 5 undescore-separated parts before the S.. thing starts
# (May not require these assumptions anymore).

# traverse Sample_ directories
for sampleID in `tail -n +2 title_file.txt | cut -f 3`
do
	directory=$1/Sample_$sampleID

	# basic sanity check
	if [ ! -d "$directory" ] || [ ! "`ls -A $directory`" ]
	then
		echo -e "No fastqs for $directory"
		echo -e "ABORTING"
		exit
	fi

	#check if barcodes in the title file and the sample sheet match
	# get the barcode for the index in the sample sheet
	index=`head $directory/SampleSheet.csv | tail -1 | cut -d ',' -f 5`

	# compare barcodes from title file and $barcodeKeyFile
	expectedBarcode=`grep -w $index $barcodeKeyFile | cut -f 2`
	actualBarcode=`grep -w $sampleID title_file.txt | cut -f 1`
	echo -e "$expectedBarcode and $actualBarcode"

	if [ "$expectedBarcode" == "" ] || [ "$actualBarcode" == "" ]
	then
		echo -e "empty barcode."
		echo -e "ABORTING"
		exit
	elif [ "$expectedBarcode" != "$actualBarcode" ]
	then
		echo -e "$sampleID does not have correct barcode in the title file."
		echo -e "$expectedBarcode != $actualBarcode"
		echo -e "ABORTING"
		exit
	fi

	# get the barcode index
	fastqStem=`awk 'BEGIN{FS=","}{if(NR==1) next; gsub(/_/, "-", $3); print $3"_"$5}' $directory/SampleSheet.csv`
	for fastq in $directory/*.fastq.gz
	do
		fileName=`basename $fastq`
		# this part works irrespective of whether you have the _S thingy or the actual barcode in fastq file name
		linkName=`echo $fileName | awk -v stem=$fastqStem '{split($0, b, "_"); sThingy="_"b[length(b)-2]"_"; split($0, c, sThingy); print stem"_"c[2]}'`
		resolvedPath=`readlink -f $fastq`
		ln -s $resolvedPath $linkName
	done
done


# change "NormalPool" to "PoolNormal", only accept the first pooled normal
awk 'BEGIN{FS=OFS="\t"; gotPooledNormal=0;}{if($6=="NormalPool" || $6=="PoolNormal"){if(gotPooledNormal==1){next;} else{$6 = "PoolNormal"; gotPooledNormal=1; print $0;}} else {print $0}}' title_file.txt > t
mv t title_file.txt

# make the valid title file from the given title file
# replace _ with - 2nd line onwards
sed -i -e '2,$s/_/-/g' title_file.txt

# sort by barcode
sort -k 1,1 title_file.txt > t
mv t title_file.txt

# define the needed variables
main=~/software/IMPACT-Pipeline/bin/RunIlluminaProcess.pl
SVConfigFile=~/resources/impact-pipeline/impact-sv.conf
perl=/opt/common/CentOS_6/perl/perl-5.20.2/bin/perl



# most of these things are required to be in the path
export PATH=/opt/common/CentOS_6/R/R-3.1.2/bin:/opt/common/CentOS_6/bedtools/bedtools-2.17.0:/opt/common/CentOS_6/java/jdk1.8.0_31/bin:/opt/common/CentOS_6/samtools/samtools-1.2:$PATH
export TMPDIR=/ifs/work/scratch
export _JAVA_OPTIONS=-Djava.io.tmpdir=/ifs/work/scratch

jobName=$3
configFile=$5
workingDir=$6


## Run IMPACT-Pipeline on LSF
if [ "$4" == "yes" ]
then
	commandString="bsub -q sol -cwd $workingDir -J $jobName -eo $jobName.stderr -oo $jobName.stdout -We 96:00 -R \"rusage[mem=16]\" -M 20 $perl $main -c $configFile -sc $SVConfigFile -d $workingDir -o $workingDir"
else
	commandString="bsub -q sol -cwd $workingDir -J $jobName -eo $jobName.stderr -oo $jobName.stdout -We 96:00 -R \"rusage[mem=16]\" -M 20 $perl $main -c $configFile -d $workingDir -o $workingDir"
fi

echo $commandString
eval $commandString
