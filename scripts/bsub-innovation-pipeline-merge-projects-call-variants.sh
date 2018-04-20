#! /bin/bash

# merge 2 or more runs that have completed step 3 and run steps 4-7 on them together.
# run this script from an empty folder.


# $1 - project/job name
# $2 - panel type: impact or hemepact.
# $3 - custom capture? yes or no
# $4 - run SV? yes or no.
# $5 onwards - a list of working directories of projects to be combined with the current project. The DMP working directory must be the last one in the list. It is possible to provide a non-existing DMP directory to avoid including DMP samples in the run.


# Process:
# for each listed working directory
# 1. add the contents of the title file to the new title file
# 2. assign a new barcode to the samples
# 3. make links to bams and other required files with new names
# 4. add new bam names to list of files


workingDir=`pwd`
resources=~/resources/impact-pipeline
listOfFiles=VariantCallingListOfFiles.txt
rm -f $listOfFiles


# Create title file
printf "Barcode\tPool\tSample_ID\tCollab_ID\tPatient_ID\tClass\tSample_type\tInput_ng\tLibrary_yield\tPool_input\tBait_version\tGender\tPatientName\tMAccession\tExtracted_DNA_Yield\n" > title_file.txt
# create dummy sample sheet
touch SampleSheet.csv

declare -A seenBarcodes
barcodeNumber="5000"
gotPooledNormal=""
jobName=$1
panel=$2
customCapture=$3
runSV=$4
shift 4

newPool="MergedPool"
DMPDir="${@: -1}"

# iterate over given working directories
for projectDir in "$@"
do
	if [[ "$projectDir" == "$DMPDir" ]]
	then
		echo -e "DMP directory: $DMPDir"

		# skip the DMP directory if it's not valid
		if [ ! -d "$projectDir" ]
		then
			echo -e "Invalid DMP directory: $projectDir"
			echo -e "Will continue without DMP samples."
			break
		fi

		# we have a proper DMP directory
		tail -n +2 $projectDir/title_file.txt > headless_title_file.txt

		# We are processing the DMP title file, remove all samples from those patients whose samples are only present in DMP bams
		# It is required that the DMP working directory be the last one in the list of working directories.
		# Also remove matched normals for those samples which have at least one matched nomral from the plates
		awk 'BEGIN{FS=OFS="\t"; while(getline < "title_file.txt"){ if(index($3, "DummyPool")==0){pid[$5]=$5; if($6=="Normal"){ matchedNormals[$5]=$3;}}}}{ if(pid[$5]==null || ($6=="Normal" && matchedNormals[$5] != null)){ next;} print}' headless_title_file.txt > t
		mv t headless_title_file.txt
	else

		# non-DMP, plate directory. just process
		tail -n +2 $projectDir/title_file.txt > headless_title_file.txt
	fi

  # read each line from the title file, split on tab (bash is pathetic!)
	while IFS=$'\t' read -r -a array
	do
		oldBarcode=${array[0]}
		oldSampleID=${array[2]}
		oldPool=${array[1]}
		newSampleID="$oldSampleID-$oldPool"

		# Change "NormalPool" to "PoolNormal"
		if [ "${array[5]}" == "NormalPool" ]
		then
			array[5]="PoolNormal"
		fi

		# Add only the very first pooled normal
		if [ "${array[5]}" == "PoolNormal" ]
		then
			if [ "$gotPooledNormal" == "true" ]
			then
				continue
			else
				gotPooledNormal="true"
			fi
		fi

		# change barcode
		if [ "${seenBarcodes[$oldBarcode]}" ]
		then
			# use a new barcode with 4 didgits in it
			newBarcode="bc$barcodeNumber"
			((barcodeNumber++))
		else
			# or use old barcode left-padded with 0's
			# pad with appropriate number of 0's to make sure it has 4 digits
			if [ "${#oldBarcode}" == "4" ]
			then
				newBarcode="bc00${oldBarcode/bc/}"
			else
				newBarcode="bc0${oldBarcode/bc/}"
			fi
		fi

		# add old barcode to seen barcodes
		seenBarcodes[$oldBarcode]=$oldBarcode

		echo -e "Barcode: $oldBarcode -> $newBarcode AND Sample: $oldSampleID -> $newSampleID"

		# add proper line to the title file
		array[0]=$newBarcode
		array[1]=$newPool
		array[2]=$newSampleID
		(IFS=$'\t'; echo -e "${array[*]}") >> title_file.txt

		# make proper links

		# bam
		file="$projectDir/FinalBams/${oldSampleID}_${oldBarcode}_${oldPool}_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam"
		if [ ! -e "$file" ]
		then
			echo -e "File does not exist: $file"
			echo -e "ABORTING"
			exit
		fi

		fullPath=`readlink -f $file`
		linkName="${newSampleID}_${newBarcode}_${newPool}_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam"
		ln -s $fullPath $linkName

		# also add bam file to the list of files to be processed
		echo -e $workingDir/$linkName >> $listOfFiles

		# bai
		file="$projectDir/FinalBams/${oldSampleID}_${oldBarcode}_${oldPool}_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bai"
		if [ ! -e "$file" ]
		then
			echo -e "File does not exist: $file"
			echo -e "ABORTING"
			exit
		fi

		fullPath=`readlink -f $file`
		linkName="${newSampleID}_${newBarcode}_${newPool}_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bai"
		ln -s $fullPath $linkName

		# stats
		file="$projectDir/AllSampleResults/${oldSampleID}_${oldBarcode}_${oldPool}_L000_R1_mrg_cl.stats"
		if [ ! -e "$file" ]
		then
			echo -e "File does not exist: $file"
			echo -e "ABORTING"
			exit
		fi

		fullPath=`readlink -f $file`
		linkName="${newSampleID}_${newBarcode}_${newPool}_L000_R1_mrg_cl.stats"
		ln -s $fullPath $linkName

		file="$projectDir/AllSampleResults/${oldSampleID}_${oldBarcode}_${oldPool}_L000_R2_mrg_cl.stats"
		if [ ! -e "$file" ]
		then
			echo -e "File does not exist: $file"
			echo -e "ABORTING"
			exit
		fi

		fullPath=`readlink -f $file`
		linkName="${newSampleID}_${newBarcode}_${newPool}_L000_R2_mrg_cl.stats"
		ln -s $fullPath $linkName

		# MD.metrics
		file="$projectDir/AllSampleResults/${oldSampleID}_${oldBarcode}_${oldPool}_L000_mrg_cl_aln_srt_MD.metrics"
		if [ ! -e "$file" ]
		then
			echo -e "File does not exist: $file"
			echo -e "ABORTING"
			exit
		fi

		fullPath=`readlink -f $file`
		linkName="${newSampleID}_${newBarcode}_${newPool}_L000_mrg_cl_aln_srt_MD.metrics"
		ln -s $fullPath $linkName

	done < headless_title_file.txt

	# not copying because we will get the bams for the same patient from multiple plates
	# and we want to account for all of them.
	# copy the "targetsToRealign" files from the project directory
	# these files are on per patient basis, just copy them in bulk, patient ids don't change
	# cp $projectDir/*_targetsToRealign_covered_srt.bed .
	# cp $projectDir/*_targetsToRealign_covered_srt.list .

done


rm headless_title_file.txt

# sort title file by barcode
sort -k 1,1 title_file.txt > t
mv t title_file.txt


# choose the right config file for the given panel type and copy it to working dir
case $panel in
impact)
	configFile=$workingDir/impact.conf
	cp $resources/impact.conf $configFile
	;;
hemepact)
	configFile=$workingDir/hemepact.conf
	cp $resources/hemepact.conf $configFile
	;;
*)
	echo -e "panel type not recognized: $4"
	exit
	;;
esac

# modify the config file as needed
sed  -i 's/1,2,3,4,5,6,7/4,5,6,7/g' $configFile
awk -v list=$listOfFiles -v project=$newPool '{if(index($1, "ListOfFiles")!=0){print "ListOfFiles = "list} else if(index($1, "ProjectName")!=0){print "ProjectName = "project} else {print}}' $configFile > t
mv t $configFile


# define the variables needed to run the pipeline
# main=~/software/IMPACT-Pipeline/bin/RunIlluminaProcess.pl
main=~/software/IMPACT-Pipeline/bin/RunIlluminaProcess.pl
SVConfigFile=~/resources/impact-pipeline/impact-sv.conf
perl=/opt/common/CentOS_6/perl/perl-5.20.2/bin/perl


# most of these things are required to be in the path
export PATH=/opt/common/CentOS_6/R/R-3.1.2/bin:/opt/common/CentOS_6/bedtools/bedtools-2.17.0:/opt/common/CentOS_6/java/jdk1.8.0_31/bin:/opt/common/CentOS_6/samtools/samtools-1.2:$PATH
export TMPDIR=/ifs/work/scratch
export _JAVA_OPTIONS=-Djava.io.tmpdir=/ifs/work/scratch

## Run IMPACT-Pipeline on LSF
if [ "$runSV" == "yes" ]
then
	commandString="bsub -q sol -cwd $workingDir -J $jobName -eo $jobName.stderr -oo $jobName.stdout -We 96:00 -R \"rusage[mem=16]\" -M 20 $perl $main -c $configFile -sc $SVConfigFile -d $workingDir -o $workingDir"
else
	commandString="bsub -q sol -cwd $workingDir -J $jobName -eo $jobName.stderr -oo $jobName.stdout -We 96:00 -R \"rusage[mem=16]\" -M 20 $perl $main -c $configFile -d $workingDir -o $workingDir"
fi

echo $commandString
eval $commandString

#echo bsub -q sol -cwd "." -eo %J.e -oo %J.o -We 24:00 -R "rusage[mem=2]" -M 4 $perl $main -c $configFile -d . -o .
# bsub -q sol -cwd "." -eo %J.e -oo %J.o -We 24:00 -R "rusage[mem=2]" -M 4 $perl $main -c $configFile -d . -o .
