#! /bin/bash

# process UMI fastq (put UMI in read name, filter etc.), and create sample folder in the current directory.

# $1 - UMI style. Currently supported: loeb and loop

# $2 - new project directory

# $3 - R1 fastq to be processed. Must have _R1_ in the name. Corresponding R2 fastq must be present in the same directory


java=/opt/common/CentOS_6/java/jdk1.8.0_31/bin/java

export TMPDIR=/ifs/work/scratch

if [ "$1" == "loeb" ]
then
  $java -server -Xms8g -Xmx8g -cp ~/software/Marianas.jar org.mskcc.marianas.umi.duplex.fastqprocessing.ProcessLoebUMIFastq $3 10 TGACT $2
elif [ "$1" == "loop" ]
then
  $java -server -Xms8g -Xmx8g -cp ~/software/Marianas.jar org.mskcc.marianas.umi.duplex.fastqprocessing.ProcessLoopUMIFastq $3 3 $2
else
  echo -e "Unrecognized UMI style: $1"
  echo -e "Aborting."
fi






#
