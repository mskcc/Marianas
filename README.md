# Marianas: Deep Next Generation Sequencing, Molecular Barcoding


This software was developed at the Innovation Lab, Center for Molecular Oncology, Memorial Sloan Kettering Cancer Center.

Marianas is a Java software that processes NGS sequencing data where reads contain molecular barcodes (UMIs). This processing involves 3 steps:
1. UMI clipping: Identify UMIs and bridge bases. ADD UMI pair to fastq read name, discard bridge bases
2. Collapsing: Process standard bam file containing aligned, UMI-clipped reads. Identify PCR duplicates by using mapping position and UMIs. Correct PCR and sequencing errors by taking consensus among PCR duplicates. Create collapsed fastq files
3. Bam separation: From the collapsed bam, create separate bam files contianing duplex-only and simplex-only reads


## Dependencies (bundled with the release)

1. BioinfoUtils
2. HTSJDK
3. Google Guava
4. Apache Commons IO
