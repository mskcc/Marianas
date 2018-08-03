#! /bin/bash

# get the read count for each sample in each project in the given run level folder

# $1 - run level folder that contains projects that contain sample fastq folders

barcodeKeyFile=~/resources/impact-pipeline/barcodeKey96.txt

runName=`basename $1`
outFile=${runName}-barcodes.txt
rm -f $outFile

for sampleDirectory in `ls -d $1/Project*/Sample_*`
do
  sample=`basename $sampleDirectory`
  sample=${sample/Sample_/}
  echo -e $sample
  barcode=`tail -1 $sampleDirectory/SampleSheet.csv | cut -d ',' -f 5`
	barcodeName=`grep -w $barcode $barcodeKeyFile | cut -f 2`

  # go through all available lanes
  for fastq in `ls $sampleDirectory/*R1*fastq.gz`
  do
    lane=`echo $fastq | awk '{split($0, a, "_"); for(i in a){if(index(a[i], "L0")==1) {print a[i]; exit}}}'`
    reads=`zcat $fastq | wc -l`
    reads=$((reads/4))
    printf "$sample\t$barcode\t$runName-$lane\t$reads\t$barcodeName\n" >> $outFile
  done
done

echo -e "Done!"


#
