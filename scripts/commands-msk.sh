
# find largest files and folders on the disk
du -am ~ > t
sort -k 1,1 -nr t > large-files.txt



# divide a bed file into 2, containing and not containg given genes
awk 'BEGIN{FS=OFS="\t"; while(getline < "custom-1-gene.txt"){a[$1]=$1}}{split($5, b, "_"); gene=b[1]; if(a[gene]!=null) print $0}' ../picard_targets.interval_list_IMPACT341_sorted.bed > custom-1-gene.bed

awk 'BEGIN{FS=OFS="\t"; while(getline < "custom-1-gene.txt"){a[$1]=$1}}{split($5, b, "_"); gene=b[1]; if(a[gene]==null) print $0}' ../picard_targets.interval_list_IMPACT341_sorted.bed > impact-minus-custom-1-gene.bed


# run the collate script to get the results
~/software/bin/custom-panel-collate-stats.sh custom-1-gene impact-minus-custom-1-gene custom-11-gene impact-minus-custom-11-gene

# add MD tag to the bam file
  samtools calmd -b  MOMO_0089_AC6G07ANXX___P5500_Y___PC11_5_A___hg19___MD.sorted.bam ~/resources/hg19-ucsc/human_hg19.fa > a.bam
  samtools index a.bam

# use samtools to create a pileup, including duplicate reads
# it still requires the reference fasta irrespective of MD tag
samtools mpileup -C50 -O -A -d 1000000 --excl-flags UNMAP,SECONDARY  -f ~/resources/hg19-ucsc/human_hg19.fa  a.bam > a-pileup-samtools.txt

# get BRAF part of pileup
cat a-pileup-samtools.txt | awk '$1=="chr7" && $2>=140430000 && $2<=140650000' | less

# sort the pileup file created by LocusPocus Metrics module
sort -k 1,1 -k 2,2n a-pileup.txt > t





# full and partial runs
from fastqs
from bams
from standard bams
from-fastqs-to-standard-bams
from-bams-to-standard-bams




# get a selected region from a bam file as another bam file
samtools view -b a.bam chr1 > b.bam


# get a nice summary of IMPACT_local.log
grep -i -e starting -e started -e finished -e running -e complete -e completed  -e Thanks IMPACT_local.log


# get the most relevant stderr files first
for f in `ls -t *.stderr`; do cat $f >> error.txt; done
for f in `ls -t StdLogFiles/*.stderr`; do cat $f >> error.txt; done


# kill all jobs
for jobID in `bjobs | tail -n +2 | cut -d ' ' -f 1`; do bkill $jobID; done
bkill -b -u 0
bkill -b -u patelju1

# modify a pending notify job so that it is finished and the pipeline moves on
bmod -wn 3466588

# use scp to copy needed files from luna
scp -r patelju1@luna.mskcc.org:/home/patelju1/projects/impact-5500-af-hodgkins-runs/Hodgkins9/\{title_file.txt,*.conf,SampleSheet.csv,IMPACT_local.log,Results,error.txt\} .


# compare BRAF coverage between igvtools and LocuPocus
# run IGVTools count function
java -Xmx8G -Xms8G -Djava.awt.headless=true -jar ~/IGVTools/igvtools.jar count -w 1 --bases --includeDuplicates --query chr7:1-150000000 ../a.bam a.wig hg19
awk '$1>=140419127 && $1<=140624564' a.wig > sample.txt
awk '$1=="chr7"' a-pileup.txt > double-capture/a-pileup-chr7.txt
awk 'BEGIN{FS=OFS="\t"; while(getline < "sample.txt"){a[$1]=$2; c[$1]=$3; g[$1]=$4; t[$1]=$5} }{if(a[$2]==null){next;} print $2, $5-a[$2], $6-c[$2], $7-g[$2], $8-t[$2], "|", $5, $6, $7, $8  }' a-pileup-chr7.txt > t


# wdiff word by word diff
/home/patelju1/software/wdiff-1.2.2/src/wdiff -3  t1 t2 | grep -e "\[-" -e "\-]" -e "\{+" -e "\+}"



# Join the PatientID-DMPID and DMPID-BamStem files to get all the necessary bam file stems and corresponding pateint ids that is the input for bam -> standard bam script
awk 'BEGIN{FS=OFS="\t"; while(getline < "PatientID-DMPID.txt"){if(NR==1) continue; split($2, a, ", "); for(i=1; i<=length(a); i++){split(a[i], b, "-"); DMPNumber=b[2]; pid[DMPNumber]=$1;}} FS=","; OFS="\t"}{split($1, b, "-"); if(b[4]=="IM5"){panel="IMPACT410"} else {panel="IMPACT341"} if(pid[b[2]]!=null){print $2, pid[b[2]], panel, $1}}' /ifs/dmprequest/12-245/key.txt | sort -k 2,2 > t

# for DMP normals, prefer IM5 bams
awk 'BEGIN{FS=OFS="\t"; while(getline < "t"){split($4, a, "-"); if(a[4]=="IM5"){IM5[$2]=$2}}}{ split($4, a, "-"); if(index($1, "-N") != 0 && a[4]=="IM3" && IM5[$2] != null ) next; print }' t > aggregated-info.txt
# create a tab delimited file of full paths to both tumor and normal bams and pateint ids from an aggregated info file as generated above
awk 'BEGIN{FS=OFS="\t"; path="/ifs/dmpshare/share/irb12_245"}{ if($1=="#N/A" || seen[$1] != null) next; print path"/"$1".bam", $2, $3; seen[$1] = $1; }' aggregated-info.txt > file-list.txt


# install/clone IMPACT Pipeline
git clone --recursive https://github.com/rhshah/IMPACT-Pipeline.git


# Changes made to IMPACT Pipeline
# change timings in the IMPACT Pipeline code.
# Only need to change files in bin and support-scripts directories
sed -i 's/0:59/24:00/g' bin/*
sed -i 's/0:59/24:00/g' support-scripts/*
# also check if there are other timings and change them as well
grep -R -- "-We " * | grep -v "24:00" | grep -- "-We"
# change $pcount comparison value from 4 to 8/10 in mutation calling??
$pcount
# set iounits to 0??
-R "rusage[iounits=0]" for some of the jobs



# in a title file, remove those samples that are only in DMP
awk 'BEGIN{FS=OFS="\t"; while(getline < "title_file.txt"){ if(index($3, "DummyPool")==0){pid[$5]=$5} }}{ if(index($3, "DummyPool")!=0 && pid[$5]==null){next;} print}' title_file.txt > t


# in a title file, print all the patient ids that do not have a matched normal
awk 'BEGIN{FS=OFS="\t";}{if(p[$5]==null || p[$5]=="Tumor"){ p[$5]=$6;}}END{ for(i in p){ if(p[i]!="Normal") print i  }}' run-5500-AB/title_file.txt

# print inter-exonic spacing in a bed file, taking into consideration + and - strand
awk '{if(chr==$1){spacing=x1-$3; if(spacing<=0){spacing=$2-y1;} print spacing"\t"interval1" AND "$5} x1=$2; y1=$3; chr=$1; interval1=$5 }' ../bedFiles/STAG2-CDKN2A-TP53-EWSR1.bed

# print the intervals in a bed file that are < 100bp apart
awk 'BEGIN{FS=OFS="\t"}{if(chr!=$1) {chr=$1; next;} spacing=x1-$3; if(spacing<0 || spacing>100) spacing=$2-y1; if(spacing >= 0 && spacing <= 100) {print spacing"\t"interval1" AND "$5} x1=$2; y1=$3; chr=$1; interval1=$5 }' ../bedFiles/AL-Switch.bed


# print a subsequence from reference fasta
samtools faidx /home/patelju1/resources/impact-GRCh37/Homo_sapiens_assembly19.fasta 1:166244805-166246460


# git commands
# view git file diffs in gui
git difftool -g
# get submodule update from upstream
git submodule update --remote Innovation-IMPACT-Pipeline
# commit from inside submodule
git commit; git push origin master
# then commit the main projects
git commit; git push origin master



# Processing IMPACT output

# make a list of mutations from the annotated_exonic_variants.txt file where each mutation is repsresented exactly once. This mutations list is used by Waltz Genotyping module
tail -n +2 /home/patelju1/projects/Cercek-Colon/run-X-DMP-1/annotated_exonic_variants.txt | awk 'BEGIN{FS=OFS="\t"}{ key=$3"\t"$4"\t"$5"\t"$6; if(m[key]!=null) next; m[key]=key; name=$8""substr($14, 2); print key"\t"name  }' > Cercek-Colon.txt

# just add patient id to annotated_exonic_variants.txt
awk 'BEGIN{FS=OFS="\t"; while(getline < "title_file.txt"){pid[$3]=$5}}{if(NR==1){$1="Sample\tPatient"; print; next;} $1=$1"\t"pid[$1]; for(i=37; i<=NF; i++){split($i, a, ";"); split(a[3], c, "="); if(c[2]==0) {$i="-"}} print}' annotated_exonic_variants.txt > annotated_exonic_variants-with-patient-id.txt


# pull genotypes for various samples for given patient from annotated_exonic_variants-with-patient-id.txt file
awk 'BEGIN{ patientName="AC-pt1"; FS=OFS="\t"; while(getline < "title_file.txt"){patient[$3]=$5}}{if(NR==1){line="Mutation"; for(i=37; i<=NF; i++){sampleName[i]=$i; if(patient[$i]==patientName){line=line"\t"$i}} print line} else {if($2!=patientName) next; line=$9""substr($15, 2); for(i=37; i<=NF; i++){sample=sampleName[i]; if(index(sample, "05500-AJ-P-1")!=0 || patient[sample]!=patientName) continue; split($i, a, ";"); split(a[1], b, "="); dp=b[2]; split(a[3], b, "="); ad=b[2]; split(a[4], b, "="); vf=b[2];  col=vf" ("ad"/"dp")"; line=line"\t"col;} print line }}' annotated_exonic_variants-with-patient-id.txt


# find called mutations in given sample (these would get highlighted in the presentation)
awk 'BEGIN{sampleName="AC-pt1-3-5500-X"; FS=OFS="\t"}{if($1==sampleName) print $9""substr($15, 2)}' annotated_exonic_variants-with-patient-id.txt




# run bcl2fastq2
~/software/bin/bsub-script.sh ~/software/bcl2fastq-ryan --sample-sheet ../../SampleSheet-hiseq.csv --runfolder-dir ../../../bcls/170214_PITT_0102_AHHN7NBBXX/ --output-dir . -d 4 -p 4 --barcode-mismatches 0

# option to get the UMI in read2
~/software/bin/bsub-script.sh ~/software/bcl2fastq-juber --sample-sheet SampleSheet-hiseq.csv --runfolder-dir ../../../bcls/170214_PITT_0102_AHHN7NBBXX/ --output-dir . -d 4 -p 4 --barcode-mismatches 0 --use-bases-mask y125,i4,y129

# get cluster count from umi-summary
cat BC11_bc11_Project-5500-AV_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam.umi-summary | cut -f 4 | sort | uniq -c > BC11-distinct-umis.txt

# sort UMIs by number of clusters
cat BC11-distinct-umis.txt | sort -k 1,1nr > BC11-distinct-umis-count-sorted.txt



# inspect genotypes in multiple bam files
~/software/bin/inspect-genotypes.sh ~/projects/PUMA/run-5500-AR/FinalBams ~/temp/contamination.bed ~/temp/mutations.txt

# copy contents from terminal without errors
cat locuspocus-genotypes-samples.txt | pbcopy
pbpaste > t (or command v)


# run pairwise, 2-sequence blast
/opt/common/CentOS_6-dev/blast/ncbi-blast-2.3.0+/bin/blastn -query t.fasta -subject adapters.fasta -outfmt 7 -num_alignments 1 -dust no -word_size 8 -evalue 1000 > t.txt








# go over assigned and unassigned barcode read count files and tabulate per barcode read numbers for a lane
awk 'BEGIN{FS=OFS="\t"; lane="PITT_0087_AHFMV3BBXX-L002"; while(getline < "CG-L002-unassigned-standard-barcodes.txt"){counts[$3]=$2;} while(getline < "PITT_0087_AHFMV3BBXX-barcodes.txt"){if($3!=lane) continue; counts[$5]=$4;}}{if(counts[$1]==null){count=0}else{count=counts[$1]} print count}' barcodes.txt


# print the read sequences that violate IDT Loop UMI-Constant rules
gunzip -c Sample_SK-PB-191-G-30-Loop_IGO_05500_CR_9/SK-PB-191-G-30-Loop_IGO_05500_CR_9_S75_L004_R1_001.fastq.gz | awk '{getline seq; getline; getline; count++; last=substr(seq, 3,1); c1=(last=="A" || last=="T") && substr(seq, 5,1)!="T"; c2=(last=="G" || last=="C") && substr(seq, 4,1)!="T"; if(c1 || c2){print count, seq}}' | head


# calculate substitution rates from a pileup file
awk 'BEGIN{FS=OFS="\t"; base[5]="A"; base[6]="C"; base[7]="G"; base[8]="T"}{max=-1; garbage=0; genotype=-1; for(i=5; i<=8; i++){if($i>max){max=$i; genotype=i}} for(i=5; i<=8; i++){if(i==genotype) continue; if($i>max*.1) next}   for(i=5; i<=8; i++){if(i==genotype) continue; substitution=base[genotype]">"base[i]; genotypeCount[substitution]+=max; altCount[substitution]+=$i;}}END{for(s in genotypeCount){print s, genotypeCount[s], altCount[s], 100*altCount[s]/(genotypeCount[s]+altCount[s]); g+=genotypeCount[s]; a+=altCount[s]} print "XAggregate", g, a, 100*a/(g+a) }' SK-PB-191-G-30-Loop-IGO-05500-CR-9_bc212_05500-CR_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt | sort -k 1,1 > G-30-collapsed-substitution-errors.txt









# ccc
