#! /bin/bash


# Create title file and other dummy files
# Run Innovation Impact pipeline from step 3 to step 4.
# Steps 3-4 include Mark Duplicates, Indel Realignment, Base Quality Recalibration, Metrics Calculation and QC Report Genaration

# $1 - tab delimited file with full paths, patient ids and panel name for bams to be processed
# $2 - project/job name
# $3 - panel type: impact or hemepact. default is impact


# Create a dummy title file for the given list of bam files
printf "Barcode\tPool\tSample_ID\tCollab_ID\tPatient_ID\tClass\tSample_type\tInput_ng\tLibrary_yield\tPool_input\tBait_version\tGender\tPatientName\tMAccession\tExtracted_DNA_Yield\n" > title_file.txt

bamFileNameSuffix="L000_mrg_cl_aln_srt.bam"
barcodeNumber="501"

resources=~/resources/impact-pipeline
workingDir=`pwd`
listOfFiles=JuberListOfFiles.txt
rm -f $listOfFiles


# read each line that has full bam file path and patient id
while read -a array
do
	bamFile=${array[0]}
	patientID=${array[1]}
	# replace all _ with - in patientID
	patientID="${patientID//_/-}"
	baitVersion=${array[2]}

	# add proper line to the title file for each bam file
	# get the sample id from bam header
	sampleID=`samtools view -H $bamFile | grep "@RG" | cut -f 2 | sed 's/ID://'`
	pool="DummyPool"
	barcode=bc$barcodeNumber

	# get sample class
	if [[ $sampleID == *-T ]]
	then
		class="Tumor"
	else
		class="Normal"
	fi

	printf "$barcode\t$pool\t$sampleID\t-\t$patientID\t$class\t-\t-\t-\t-\t$baitVersion\t-\t-\t-\t-\n" >> title_file.txt

	# create local links for the bam file and its index
	linkName=$sampleID"_"$barcode"_"$pool"_"$bamFileNameSuffix
	fullPathBam=`readlink -f $bamFile`
	fullPathBai=`readlink -f ${bamFile/.bam/.bai}`
	if [ ! -e "$fullPathBam" ] || [ ! -e "$fullPathBai" ]
	then
		echo -e "Bam or index does not exist"
		echo -e "$fullPathBam"
		echo -e "ABORTING"
		exit
	fi
	ln -s $fullPathBam $linkName
	ln -s $fullPathBai ${linkName/.bam/.bai}

	# add the validly named bam to the list of files that will be supplied to the pipeline
	echo -e $workingDir/$linkName >> $listOfFiles

	# for each bam file, copy the dummy files needed for the successful run from step 3 to step 4
	# cp $resources/H-1-IGO-05500-AF-27_bc81_Proj-5500-AF_L000_mrg_cl_aln_srt_MD.metrics ${linkName/.bam/_MD.metrics}
	cp $resources/H-1-IGO-05500-AF-27_bc81_Proj-5500-AF_L000_R1_mrg_cl.stats ${linkName/_mrg_cl_aln_srt.bam/_R1_mrg_cl.stats}
	cp $resources/H-1-IGO-05500-AF-27_bc81_Proj-5500-AF_L000_R2_mrg_cl.stats ${linkName/_mrg_cl_aln_srt.bam/_R2_mrg_cl.stats}

	((barcodeNumber++))
done < $1

# sort title file by barcode
sort -k 1,1 title_file.txt > t
mv t title_file.txt

# choose the right config file for the given panel type and copy it to working dir
case $3 in
impact|"")
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
sed  -i 's/1,2,3,4,5,6,7/3,4/g' $configFile
awk -v list=$listOfFiles -v project=$pool '{if(index($1, "ListOfFiles")!=0){print "ListOfFiles = "list} else if(index($1, "ProjectName")!=0){print "ProjectName = "project} else {print}}' $configFile > t
mv t $configFile

# create dummy sample sheet
touch SampleSheet.csv


# define the variables needed to run the pipeline
main=~/software/IMPACT-Pipeline/bin/RunIlluminaProcess.pl
svConfigFile=$resources/template_sv.conf
perl=/opt/common/CentOS_6/perl/perl-5.20.2/bin/perl



# most of these things are required to be in the path
export PATH=/opt/common/CentOS_6/R/R-3.1.2/bin:/opt/common/CentOS_6/bedtools/bedtools-2.17.0:/opt/common/CentOS_6/java/jdk1.8.0_31/bin:/opt/common/CentOS_6/samtools/samtools-1.2:$PATH
export TMPDIR=/ifs/work/scratch
export _JAVA_OPTIONS=-Djava.io.tmpdir=/ifs/work/scratch

jobName=$2

## Run IMPACT-Pipeline on LSF
commandString="bsub -q sol -cwd $workingDir -J $jobName -eo $jobName.stderr -oo $jobName.stdout -We 48:00 -R \"rusage[mem=16]\" -M 20 $perl $main -c $configFile -d $workingDir -o $workingDir"
echo $commandString
eval $commandString

#
