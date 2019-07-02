# Marianas: Deep Next Generation Sequencing, Molecular Barcoding


This software was developed at the Innovation Lab, Center for Molecular Oncology, Memorial Sloan Kettering Cancer Center (MSKCC).

Marianas is a Java software that processes NGS sequencing data where reads contain molecular barcodes (UMIs). This processing involves 3 steps:
1. UMI clipping: Identify UMIs and bridge bases. ADD UMI pair to fastq read name, discard bridge bases
2. Collapsing: Process standard bam file containing aligned, UMI-clipped reads. Identify PCR duplicates by using mapping position and UMIs. Correct PCR and sequencing errors by taking consensus among PCR duplicates. Create collapsed fastq files
3. Bam separation: From the collapsed bam, create separate bam files contianing duplex-only and simplex-only reads



## Java

Java 1.8 or above is required.

## Dependencies (bundled with the release jar)

1. BioinfoUtils
2. HTSJDK
3. Google Guava
4. Apache Commons IO


## Using Marianas

**1. UMI clipping**

java -server -Xms8g -Xmx8g -cp Marianas.jar org.mskcc.marianas.umi.duplex.fastqprocessing.ProcessLoopUMIFastq read1.fastq.gz read2.fastq.gz umi-length

This will create a pair of _umi-clipped.fastq.gz files in the current working directory


**2. Collapsing and making consensus reads**


Collapsing involves 3 sub-steps: First pass that collapses the "left" read cluster of the read pairs, followed by a sort step on the first pass output, followed by second pass that collapses the "right" read cluster of the read pairs.

  
java -server -Xms8g -Xmx8g -cp Marianas.jar org.mskcc.marianas.umi.duplex.DuplexUMIBamToCollapsedFastqFirstPass bam-file pileup-file minMappingQuality minBaseQuality mismatches wobble minConsensusPercent reference-fasta-file

sort -S 8G -k 6,6n -k 8,8n first-pass.txt > first-pass.mate-position-sorted.txt

java -server -Xms8g -Xmx8g -cp Marianas.jar org.mskcc.marianas.umi.duplex.DuplexUMIBamToCollapsedFastqSecondPass bam-file pileup-file minMappingQuality minBaseQuality mismatches wobble minConsensusPercent reference-fasta-file


**3. Separating simplex and duplex bams from original collapsed bam**

java -server -Xms8g -Xmx8g -cp Marianas.jar org.mskcc.marianas.umi.duplex.postprocessing.SeparateBams collapsed.bam





