#! /bin/sh

# $1 - bam folder
# $2 - tumor-normal pairings file (typically MergedPool_NormalUsedInMutationCalling.txt)


python=/opt/common/CentOS_6-dev/python/python-2.7.11/bin/python
iCallSV=/home/shahr2/git/iCallSV/iCallSV/iCallSV.py
configFile=~/resources/icallsv/icallsv.conf


mkdir iCallSVOutput
workingDir=`readlink -f iCallSVOutput`


while IFS=$'\t' read -r -a array
do
  tumorStem=${array[0]}
  normalStem=${array[1]}
  jobName="$tumorStem-iCallSV"
  tumorBam=`ls $1/$tumorStem*.bam`
  tumorBam=`readlink -f $tumorBam`
  normalBam=`ls $1/$normalStem*.bam`
  normalBam=`readlink -f $normalBam`
  tumorSample=`samtools view -H $tumorBam | grep "@RG" | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF; i++) {if(index($i, "SM:")==1){split($i, a, ":"); print a[2]; exit 0}}}'`
  normalSample=`samtools view -H $normalBam | grep "@RG" | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF; i++) {if(index($i, "SM:")==1){split($i, a, ":"); print a[2]; exit 0}}}'`


  commandString="bsub -q sol -cwd $workingDir -J $jobName -eo $jobName.stderr -oo $jobName.stdout -We 48:00 -M 16 -n 9 $python $iCallSV -v -sc $configFile -abam $tumorBam -bbam $normalBam -aId $tumorSample -bId $normalSample -o $workingDir -op $tumorStem"

  echo $commandString
  eval $commandString

done < $2
