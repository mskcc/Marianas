#! /bin/bash

# $1 - ref fasta
# $2 - query fasta


bsub -eo %J.e -oo %J.o -cwd "." -R "rusage[mem=16]" -We 0:60 /opt/common/CentOS_6-dev/blast/ncbi-blast-2.3.0+/bin/blastn 












#
