#! /bin/bash

# for a hung/crashed but still practically complete stage-1 run, move the files to their proper folders

# No arguments, will work on pwd

echo -e "FinalBams..."
mkdir FinalBams
mv *.bam *.bai FinalBams

echo -e "Results..."
mkdir Results
mv *.html *.rhtml figure *_ALL_* *_copynumber_segclusp.nvn.pdf *_copynumber_segclusp.pdf *_loessnorm.pdf *_normvsnorm_segclusp.pdf Results

echo -e "AllSampleResults..."
mkdir AllSampleResults
mv *_L000_R1_mrg_cl.stats *_L000_R2_mrg_cl.stats *_L000_mrg_cl_aln_srt_MD.metrics *.pdf AllSampleResults

echo -e "StdLogFiles..."
mkdir StdLogFiles
mv *.stdout *.stderr StdLogFiles

