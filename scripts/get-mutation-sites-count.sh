#!/bin/sh


awk 'BEGIN{FS=OFS="\t"}{if($4<500) next; max=1; count=0; for(i=5; i<=8; i++){if($i>$max) max=i} for(i=5; i<=8; i++){if(i==max) continue; f=$i/$4; if($i>=2 && f >= 0.001 && f < 0.15) count++} if(count>0) print}' /home/patelju1/projects/Juber/bam-metrics/custom-panel-5500-CZ-collapsed-bams/PC41-PC55-1to10-IGO-05500-CZ-2_bc210_Pool-05500-CZ-Tube1-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt | wc -l
