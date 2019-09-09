#! /bin/bash

# process UMI fastq (put UMI in read name, filter etc.), and create sample folder in the current directory.

# $1 - UMI style. Currently supported: loop

# $2 - R1 fastq to be processed.

# $3 - R2 fastq to be processed.


java=/opt/common/CentOS_6/java/jdk1.8.0_31/bin/java

export TMPDIR=/ifs/work/scratch

if [ "$1" == "loop" ]
then
  $java -server -Xms8g -Xmx8g -cp ~/software/Marianas.jar org.mskcc.marianas.umi.duplex.fastqprocessing.ProcessLoopUMIFastq $2 $3 3
else
  echo -e "Unrecognized UMI style: $1"
  echo -e "Aborting."
fi






#
