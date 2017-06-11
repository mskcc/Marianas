#!/bin/bash

# submit any script and its arguments on the luna cluster with the following job parameters

# $* - path to the script followed by script arguments


bsub -q sol -cwd "." -R "rusage[mem=8]" -We 24:00 -o %J.o -e %J.e $*


#
