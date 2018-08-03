#! /bin/sh

# $1 - bam folder
# $2 - tumor-normal pairings file (typically MergedPool_NormalUsedInMutationCalling.txt)


python=/opt/common/CentOS_6-dev/python/python-2.7.11/bin/python
configManta=/home/patelju1/Downloads/manta-1.0.1.centos5_x86_64/bin/configManta.py
referenceFasta=/ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta
bamFolder=`readlink -f $1`
pairingFile=`readlink -f $2`


mkdir MantaOutput
workingDir=`readlink -f MantaOutput`
cd MantaOutput

while IFS=$'\t' read -r -a array
do
  tumorStem=${array[0]}
  normalStem=${array[1]}
  jobName="$tumorStem-Manta"
  tumorBam=`ls $bamFolder/$tumorStem*.bam`
  tumorBam=`readlink -f $tumorBam`
  normalBam=`ls $bamFolder/$normalStem*.bam`
  normalBam=`readlink -f $normalBam`
  tumorSample=`samtools view -H $tumorBam | grep "@RG" | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF; i++) {if(index($i, "SM:")==1){split($i, a, ":"); print a[2]; exit 0}}}'`
  normalSample=`samtools view -H $normalBam | grep "@RG" | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF; i++) {if(index($i, "SM:")==1){split($i, a, ":"); print a[2]; exit 0}}}'`

  mkdir $tumorStem
  runDir=`readlink -f $tumorStem`
  $python $configManta --exome --normalBam $normalBam --tumorBam $tumorBam --referenceFasta $referenceFasta --runDir $runDir

  commandString="bsub -q sol -cwd $tumorStem -J $jobName -eo $jobName.stderr -oo $jobName.stdout -We 48:00 -M 16 $python $runDir/runWorkflow.py -m local -j 8"

  echo $commandString
  eval $commandString

done < $pairingFile
