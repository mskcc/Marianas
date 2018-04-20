#! /bin/bash

# make a title file from a project fastq directory

# $1 - project dir that has sample dirs that in turn contain samples fastq.gz and SampleSheet.csv files


pool=`basename $1`
titleFile=$pool-title-file.txt
class="Tumor"
baitVersion="IMPACT410"
patientNumber="1"
#BarcodeKey="/ifs/e63data/bergerm1/Resources/DMP/pubdata/illumina/VERSIONS/v2.0/barcodeKey96.txt"
BarcodeKey="$HOME/resources/impact-pipeline/barcodeKey96.txt"

# Create a title file
printf "Barcode\tPool\tSample_ID\tCollab_ID\tPatient_ID\tClass\tSample_type\tInput_ng\tLibrary_yield\tPool_input\tBait_version\tGender\tPatientName\tMAccession\tExtracted_DNA_Yield\n" > $titleFile


# for each sample, add a line to title file
for sampleDir in `ls -d $1/Sample*`
do
	sampleID=${sampleDir/Sample_/}
	sampleID=`basename $sampleID`
	barcodeSequence=`tail -1 $sampleDir/SampleSheet.csv | cut -d ',' -f 5`
	barcode=`grep -w $barcodeSequence $BarcodeKey | cut -f 2`
	patientID="Patient-${patientNumber}"

	printf "$barcode\t$pool\t$sampleID\t-\t$patientID\t$class\t-\t-\t-\t-\t$baitVersion\t-\t-\t-\t-\n" >> $titleFile

	((patientNumber++))
done
