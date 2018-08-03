
#run custom walker/annotator
java -Xms1g -Xmx1g  -cp ../../gatk/GenomeAnalysisTK-3.0-0/GenomeAnalysisTK.jar:bin/  org.broadinstitute.sting.gatk.CommandLineGATK -T PedigreeAnnotator -R ~/resources/hg/human_g1k_v37.fasta -sn F_C13_32F -sn M_C13_32M -sn SA_C13_32a1 -env --variant trio/C13_32F_C13_32M_C13_32a1.variantCalls.vcf -o test-out.vcf

or

java -Xms1g -Xmx1g  -cp ../../gatk/GenomeAnalysisTK-3.0-0/GenomeAnalysisTK.jar:pedigree-annotator.jar  org.broadinstitute.sting.gatk.CommandLineGATK -T PedigreeAnnotator -R ~/resources/hg/human_g1k_v37.fasta -sn F_C13_32F -sn M_C13_32M -sn SA_C13_32a1 -env --variant trio/C13_32F_C13_32M_C13_32a1.variantCalls.vcf -o test-out.vcf


#print lines where PedCat field is missing
awk 'BEGIN{FS=OFS="\t"}{if(index($0, "#")==1) next; split($9, c, ":"); split($12, b, ":"); if(length(b)!=length(c)) print $0}' test-out.vcf | head


#SNAP indexing


#SNAP alignment
snap paired hg-index/ R1.fastq R2.fastq -o out.sam -t 40

bwa mem -t 40 -M human_g1k_v37.fasta R1.fastq R2.fastq




#HaplotypeCaller qsub command
qsub -pe smp 8 -q prod.q /nethome/avabhyankar/resources/scripts/jointcalls/HaplotypeCaller_SingleSample_GVCF.v2.sh /data/analysis/NYGC/Project_NYGC_01354_WGS/Sample_NA12878_Nano_350_H0111ALXX_L004/analysis/NA12878_Nano_350_H0111ALXX_L004.final.bam NA12878_Nano_350_H0111ALXX_L004

#Create list of raw gvcfs to be genotyped
ls /data/analysis/RubinMA/Project_RUB_01238_WGS/Sample_*/analysis/*.haplotypeCalls.raw.gvcf.vcf.gz > gvcf.list.2014-08-04

#Joint Genotyping command. List must have extension .list 
qsub -pe smp 10 -q prod.q /nethome/avabhyankar/resources/scripts/jointcalls/GenotypeGVCFs.withList.v3.sh /data/analysis/RubinMA/Project_RUB_01238_WGS/gvcf.list.2014-08-04.list 2014-08-04_RubinMA

#melt the recalibrated vcfs back into respective subfolders
/nethome/jpatel/workspace/util/concordance/bulk-melt-vcfs.sh /data/analysis/avabhyankar/HaplotypeCaller/GenotypedGVCFs/genotyped/2014-08-04_RubinMA.recalibrated_variants.vcf.gz .




#run a java class from util.jar using proper classpath for mutliple jar files
java -cp util.jar:/Users/jpatel/git/htsjdk/htsjdk.jar:/Users/jpatel/git/htsjdk/lib/apache-ant-1.8.2-bzip2.jar org.nygenome.variantquality.AddStrandRDT

java -Xms2g -Xmx2g -cp util.jar:/Users/jpatel/git/htsjdk/htsjdk.jar:/Users/jpatel/git/htsjdk/lib/apache-ant-1.8.2-bzip2.jar org.nygenome.variantquality.AddStrandRDToVCF "/Users/jpatel/git/util/util/COL0028.variantCalls.vcf"  "/Users/jpatel/git/util/util/COL0028.variantCalls.StrandDP.vcf"  "/Users/jpatel/resources/hg/human_g1k_v37.fasta"   "/Users/jpatel/git/util/util/COL0028.final.bam"


#run single jar application util.jar
java -Xms2g -Xmx2g -cp util.jar org.nygenome.variantquality.AddStrandDPToVCF "/Users/jpatel/git/util/util/COL0028.variantCalls.vcf"  "/Users/jpatel/git/util/util/COL0028.variantCalls.StrandDP.vcf"  "/Users/jpatel/resources/hg/human_g1k_v37.fasta"   "/Users/jpatel/git/util/util/COL0028.final.bam"


#run julia script to extract fastq from fast5
julia ../fast5tofastq.jl *.fast5 > a.fastq


#ONT read alignment with usearch
mkdir splits
cd splits/
split ../lambda-lowtemp.fastq -l 4 -a 5 sample-reads.part
for f in `ls *sample-reads*`; do awk '{n++; if(n%4==1) print ">"$0; else if(n%4==2) print}' $f > $f.fasta; usearch --usearch_global $f.fasta --maxaccepts 1 --id 0.20 --minseqlength 64 --wordlength 4 --uc $f.uc --maxrejects 32 --strand both --db ../../../lambda-reference/lambda-reference.udb --log $f.log --alnout $f.aln; done
cd ..
cat splits/*.uc > lambda-lowtemp.uc
cat splits/*.aln > lambda-lowtemp.aln



 
#mount isilon locally using sshfs
sshfs jpatel@pennstation.nygenome.org:/ isilon/


#bash read a file
while read sampleName
do
	echo $sampleName
done < "sample_list.txt"


#see all the samples in a chip file
awk 'BEGIN{FS=OFS="\t"}{if(a[$2]==null) {print $2; a[$2]=$2 }}' /data/sequence/ILLUMINA_HISCAN/Concordance/CLIA_34samples_Exomechip_FinalReport.txt

#single-sample self-concordance
java -server -Xms3g -Xmx3g -cp ~/software/util.jar org.nygenome.concordance.SeqChipConcordance ../split-genotypes/09C79932.HumanExome-12v1.txt ../vcfs/09C79932-concordance.vcf.gz ~/resources/hg/human_g1k_v37.fasta /nethome/jpatel/resources/markers/exome-chip-markers.interval_list

#multi-sample self-concordance


#all-to-all concordance



#run verifyBamID to estmate contamination
qsub -cwd -p -90 -b y -q prod.q  ~/software/verifyBamID_1.1.0/verifyBamID/bin/verifyBamID --vcf ~/resources/verifybamid/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf.gz --bam /data/analysis/ReichardtL/Project_REI_01344_WGS/Sample_SSC00090/analysis/SSC00090.old.final.bam --out SSC00090.old.final.bam --verbose --ignoreRG


#run Vlada's script to split GenomeStudio chip file
~/software/bin/split-genotypes-by-sample.py /data/sequence/ILLUMINA_HISCAN/Concordance/chip-file.txt


#get contamination values from VerifyBamID output for all samples
cat *.selfSM | awk '{if(index($0, "#")==1) next; print $1"\t"$7}' | sort -k 2,2 -nr > contamination.txt

#print contamination values in sample list order
awk 'BEGIN{FS=OFS="\t"; while(getline < "contamination.txt") conc[$1]=$2; }{ print conc[$1] }' sample-lis

#collect self concordance in a tabular format from self concordance
cat *.self-concordance.txt | awk 'BEGIN{FS=OFS="\t"; print "Comparison\tTotalConcordance\tVariantsConcordance\tRelationship"}{ if($1=="Total(including non-variants)") {comparison=$2; totalConcordance=$6; relationship="Self"} else if($1=="PASS(Variants Only)-cumulative"  || $1==".(Variants Only)-cumulative"){print comparison"\t"totalConcordance"\t"$6"\t"relationship } }' | sort -k 3,3 -nr > self-concordance.dat

#print concordance values in sample list order
awk 'BEGIN{FS=OFS="\t"; while(getline < "self-concordance.dat") conc[$1]=$3; }{ gsub(/^[ \t]+|[ \t]+$/, "", $1); print conc[$1] }' sample-list


#collect concordance in a tabular format from pairwise concordance
cat pair*-concordance.txt | awk 'BEGIN{FS=OFS="\t"; print "Chip\tSeq\tConcordance"}{ if($1=="Total(including non-variants)") {chip=$2; seq=$3} else if($1=="PASS(Variants Only)-cumulative"  || $1==".(Variants Only)-cumulative"){print chip"\t"seq"\t"$6} }' | sort -k 3,3 -nr > pairwise-concordance.dat


#combine concordance and contamination results and print them in sample list order
awk 'BEGIN{FS=OFS="\t"; while(getline < "self-concordance.dat") conc[$1]=$3; while(getline < "contamination.txt") cont[$1]=$2; }{ print $1, conc[$1], cont[$1] }' sample-list-1.txt


#collect concordance in a tabular format from all-to-all concordance
cat *-to-all-concordance.txt > all-to-all-concordance-file.txt
awk 'BEGIN{FS=OFS="\t"; print "Comparison\tTotalConcordance\tVariantsConcordance\tRelationship"}{ if($1=="Total(including non-variants)") {comparison=$2"/"$3; totalConcordance=$6; if($2==$3) relationship="Self"; else relationship="Unrelated"} else if($1=="PASS(Variants Only)-cumulative"  || $1==".(Variants Only)-cumulative"){print comparison"\t"totalConcordance"\t"$6"\t"relationship } }' all-to-all-concordance-file.txt > all-to-all-concordance.dat




#check if final.bam exists for all samples
for d in `ls -d Sample_SSC*`; do if ! [ -e $d/analysis/*.final.bam  ]; then echo $d; fi; done


#rename multiple files using a pattern
for f in `ls Sample_PI-5340*/analysis/*recalibrated_variants.vcf.gz`; do mv $f ${f/.recalibrated_variants.vcf.gz}.recalibrated.haplotypeCalls.vcf.gz; done
for f in `ls Sample_PI-5340*/analysis/*recalibrated_variants.vcf.gz.tbi`; do mv $f ${f/.recalibrated_variants.vcf.gz.tbi}.recalibrated.haplotypeCalls.vcf.gz.tbi; done


#get numbers for expected and reported by verifyBamID contamination on the ladder of contamination
echo -e "file\texpected\tverifyBamID" > verifyBamID-performance.dat
for fileName in `ls *.selfSM`; do expected=`echo ${fileName/.final.bam.selfSM} | cut -d'_' -f3`; verifyBamID=`awk '{if(index($0, "#")==1) next; print $7}' $fileName`;  echo -e "$fileName\t$expected\t$verifyBamID";  done > verifyBamID-performance.dat

#SGE qsub alter queue details
qalter -u jpatel -q reanal.q


#replace text in a file
sed -i 's/old-word/new-word/g' a.txt
#on mac -i needs argument
sed -i '' 's/old-word/new-word/g' a.txt


#Descendence run

#create compact genotypes collection
time java -server -Xms20g -Xmx20g -cp ~/software/descendence.jar org.nygenome.descendence.genotypes.GenotypesProcessor "pedigrees/REI-40-JG.fam" "/nethome/jpatel/resources/1000-genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz" 0.3 0.000000001

#run Descendence on the compact genotypes collection
time java -server -Xms20g -Xmx20g -cp ~/software/descendence.jar org.nygenome.descendence.DescendenceStandalone "pedigrees/REI-40-JG.fam" "CompactGenotypesCollection.ser" 8 6.0 6.0 6.0 300 0.000000001


#seems the best so far, but wait for latest pedigree information
time java -server -Xms120g -Xmx120g -cp ~/software/descendence.jar org.nygenome.descendence.DescendenceStandalone "pedigrees/REI-160-JG.fam" "/nethome/jpatel/resources/1000-genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz" 10 0.3 40 0 0 80 0.000000001 0.00000001

#Descendence remote debugging
time java -Xdebug -agentlib:jdwp=transport=dt_socket,address=8000,server=y,suspend=y -server -Xms20g -Xmx20g -cp ~/software/descendence.jar org.nygenome.descendence.DescendenceStandalone "pedigrees/REI-40-JG.fam" "/nethome/jpatel/resources/1000-genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz" 32 0.3 50 0 0 30 0.000000001 0.00000001


#remove bad markers in chip seq concordance
cat mismatches-* > all-mismatches.txt

awk 'BEGIN{}{ if($1!="M") next; id=$2":"$3; if(count[id]==null) count[id]=1; else count[id]=count[id]+1; }END{for (i in count) print i"\t"count[i]}' all-mismatches.txt | sort -k 2,2 -nr  > M-frequency.txt

awk 'BEGIN{ while(getline < "M-frequency.txt") { split($1, b, ":"); bad[$1"-"b[2]] = $2;}}{ if(bad[$1]<130) print}' ~/resources/markers/exome-chip-markers.interval_list > selected-marker-2.interval_list




# git overwrite local repo with the latest copy of the remote repo
git fetch --all
git reset --hard origin/master


# awk trim white spaces
awk '{ gsub(/^[ \t]+|[ \t]+$/, ""); print }' sample-list.txt


# add pedigree information to all-to-all-concordance.dat 
# use Jave class AddPedigreeInfoToConcordance in Descendence


# run LocusPocus via qsub
qsub -pe smp 8 -cwd -p -90 -b y -q reanal.q java -server -Xms16g -Xmx16g -cp ~/software/LocusPocus.jar org.nygenome.locuspocus.LocusPocus /data/analysis/ReichardtL/Project_REI_01344_WGS/Sample_SSC12236/analysis/SSC12236.final.bam 8


# search with GenomeWhiff

# create signature from VCF
java -server -cp ~/software/GenomeWhiffServer.jar org.nygenome.genomewhiff.search.CreateSignatureFromVCF UnifiedSignatureType.ser /data/benchmark/Contamination/Project_Contamination_Genomes/Sample_NA12891_NA12892_0.04/analysis/NA12891_NA12892_0.04.haplotypeCalls.raw.gvcf.aa.vcf.gz

# create signature from GenomeStudio file
java -server -cp ~/software/GenomeWhiffServer.jar org.nygenome.genomewhiff.search.CreateSignatureFromGSFile UnifiedSignatureType.ser /data/sequence/ILLUMINA_HISCAN/Concordance/GLE_10473_Exome_4samples_12June2015.txt

# simple search
java -server -Xms8g -Xmx8g -cp ~/software/GenomeWhiffServer.jar org.nygenome.genomewhiff.search.Search DigitalKindredIndex-Unified.ser UnifiedSignatureType.ser partial/SSC05726-0.06.gws simple true 30 3000

# composite search
java -server -Xms8g -Xmx8g -cp ~/software/GenomeWhiffServer.jar org.nygenome.genomewhiff.search.Search DigitalKindredIndex-Unified.ser UnifiedSignatureType.ser /data/benchmark/Contamination/Project_Contamination_Genomes/Sample_NA12891_NA12892_0.04/analysis/NA12891_NA12892_0.04.haplotypeCalls.raw.gvcf.aa.vcf.gz.gws composite 15892 100



# add the last dummy column to a pedigree ped/fam file
cd /nethome/jpatel/workspace/Malefactor
awk 'BEGIN{FS=OFS="\t"}{print $0"\txxx"}' /data/analysis/GleesonJ/Project_GLE_10546_Exome/GLE_10546.ped > clia-trios.ped

# run Malefactor annotation/prioritization
java -server -Xms4g -Xmx4g -cp ~/software/Malefactor.jar org.nygenome.malefactor.annotator.PedigreeAnnotator /data/analysis/avabhyankar/HaplotypeCaller/GenotypedGVCFs/genotyped/2015-06-03_OST-CLIA-GLE_Exome.recalibrated_variants.vcf.gz GLE_10546.ped intervals.interval_list trio




# get valid vcfs from total-indexed-samples.txt
awk 'BEGIN{FS="|"} {if(index($1, "/data/analysis")==1 && index($1, "WGS")!=0){ if($1 ~ /vcf$/) {print $1} else if($1 ~ /\.gz$/){gsub(/\.gz$/, "", $1); print $1}}}' total-indexed-samples.txt > t


# .bashrc useful alias
alias grep='grep --color'
alias ll='ls -lh -G'

# kill a unix process
kill -9 pid

# replace or delete tokens in files/strings (tr can be very useful for deleting, replacing, squeezing etc.)
tr '[:digit:]' ' ' a.txt

# from a list of files, check if files exist
while IFS=$'\t' read -r -a ar; do if [ ! -e ${ar[0]} ]; then echo -e "${ar[0]} does not exist"; fi; done < t1




