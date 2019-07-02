/**
 * 
 */
package org.mskcc.marianas.umi.duplex;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.mskcc.marianas.util.ClusterCollectionBuilder;
import org.mskcc.marianas.util.StaticResources;

import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

/**
 * @author Juber Patel
 * 
 *         Generate UMI metrics from a bam file
 *
 */
public class DuplexUMIBamToCollapsedFastqSecondPass
{
	/**
	 * @param args
	 *            args[0] - bam file
	 *            args[1] - pileup file
	 *            args[2] - min mapping quality
	 *            args[3] - min base quality
	 *            args[4] - UMI allowed mismatches
	 *            args[5] - UMI allowed wobble
	 *            args[6] - min consensus percent
	 *            args[7] - reference fasta
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception
	{
		File bamFile = new File(args[0]);
		File pileupFile = new File(args[1]);
		int minMappingQuality = Integer.parseInt(args[2]);
		int minBaseQuality = Integer.parseInt(args[3]);
		int UMIMismatches = Integer.parseInt(args[4]);
		int wobble = Integer.parseInt(args[5]);
		int minConsensusPercent = Integer.parseInt(args[6]);
		// set the reference fasta
		File refFastaFile = new File(args[7]);
		File refFastaIndexFile = new File(args[7] + ".fai");
		FastaSequenceIndex refFastaIndex = new FastaSequenceIndex(
				refFastaIndexFile);
		IndexedFastaSequenceFile refFasta = new IndexedFastaSequenceFile(
				refFastaFile, refFastaIndex);
		new StaticResources(refFasta, refFastaIndex, pileupFile, null);

		// no args after this point

		long start = System.currentTimeMillis();

		File firstPassFile = new File("first-pass.mate-position-sorted.txt");
		File altAlleleFile = new File("second-pass-alt-alleles.txt");
		File insertionsFile = new File("second-pass-insertions.txt");

		System.out.println("Marianas " + StaticResources.version);
		System.out.println("Second Pass");
		System.out.println("Processing " + firstPassFile.getAbsolutePath()
				+ " to produce fastqs");

		secondPass(bamFile, minMappingQuality, minBaseQuality, UMIMismatches,
				wobble, minConsensusPercent, firstPassFile, altAlleleFile,
				insertionsFile);

		long end = System.currentTimeMillis();
		System.out.println("Finished processing in " + ((end - start) / 1000)
				+ " seconds.");
	}

	private static void secondPass(File bamFile, int minMappingQuality,
			int minBaseQuality, int mismatches, int wobble,
			int minConsensusPercent, File firstPassFile, File altAlleleFile,
			File insertionsFile) throws Exception
	{
		BufferedReader firstPassReader = new BufferedReader(
				new FileReader(firstPassFile));

		BufferedWriter fastq1 = new BufferedWriter(
				new FileWriter(new File("collapsed_R1_.fastq")));
		BufferedWriter fastq2 = new BufferedWriter(
				new FileWriter(new File("collapsed_R2_.fastq")));
		BufferedWriter altAlleleWriter = new BufferedWriter(
				new FileWriter(altAlleleFile));
		BufferedWriter insertionsWriter = new BufferedWriter(
				new FileWriter(insertionsFile));

		ClusterCollectionBuilder clusterBuilder = new ClusterCollectionBuilder(
				bamFile, minMappingQuality, minBaseQuality, mismatches, wobble,
				minConsensusPercent, false);

		DuplicateReadClusterCollection clusterCollection = clusterBuilder
				.next();
		String line = null;
		String[] words = null;
		int wantedContigIndex = -1;
		int wantedStartPosition = -1;
		int totalClusters = 0;
		int badClusters = 0;
		int goodClusters = 0;
		String consensusSequenceInfo = null;

		// go over first pass clusters one by one
		// find corresponding second pass clusters
		while ((line = firstPassReader.readLine()) != null)
		{
			// words[2] = UMI, words[5] = mateContigIndex, words[6] =
			// mateContig, words[7] =
			// mateStartPosition
			words = line.split("\t");
			wantedContigIndex = Integer.parseInt(words[5]);
			wantedStartPosition = Integer.parseInt(words[7]);
			totalClusters++;

			// not applying filters at this stage anymore
			// if (!goodToWriteDuplex(words))
			// {
			// badClusters++;
			// continue;
			// }

			// this is specifically being done for the the subset bams used
			// in testing that may have only one read from pairs near the
			// end. Not needed for full bams
			// end of bam file
			if (clusterCollection == null)
			{
				break;
			}

			int comp = compareGenomicPositions(wantedContigIndex,
					wantedStartPosition, clusterCollection.getContigIndex(),
					clusterCollection.getStartPosition());

			// pass1 cluster mate position is greater than current pass2 cluster
			// position
			// skip forward
			while (comp > 0)
			{
				clusterCollection = clusterBuilder.next();

				// this is specifically being done for the the subset bams used
				// in testing that may have only one read from pairs near the
				// end. Not needed for full bams
				// end of bam file
				if (clusterCollection == null)
				{
					break;
				}

				comp = compareGenomicPositions(wantedContigIndex,
						wantedStartPosition, clusterCollection.getContigIndex(),
						clusterCollection.getStartPosition());
			}

			// positions match!
			if (comp == 0)
			{
				DuplicateReadCluster[] clusters = clusterCollection
						.getProcessedClusters();

				// iterate and try to match UMI
				for (int i = 0; i < clusters.length; i++)
				{
					// UMIs match!
					if (clusters[i] != null
							&& clusters[i].getUMI().equals(words[2]))
					{
						try
						{
							consensusSequenceInfo = clusters[i].collapseMe(
									altAlleleWriter, insertionsWriter, false);
						}
						catch (Exception e)
						{
							System.err.println(
									"Problem creating consensus sequence info for cluster:");
							System.err.println(clusters[i].getContig() + ":"
									+ clusters[i].getStartPosition() + ":"
									+ clusters[i].getUMI());
							e.printStackTrace();
							continue;
							// System.exit(1);
						}

						// there is meaningful consensus sequence info
						if (consensusSequenceInfo != null)
						{
							String[] words2 = consensusSequenceInfo.split("\t");
							writeReadPair(words, words2, fastq1, fastq2);
							goodClusters++;
							break;
						}
					}
				}
			}
			else if (comp < 0)
			{
				// missed the cluster !!!

				// TODO Check how many we are missing !!!
				// TODO throw an exception

				int a = 5;
			}
		}

		clusterBuilder.close();
		firstPassReader.close();
		fastq1.close();
		fastq2.close();
		altAlleleWriter.close();
		insertionsWriter.close();

		clusterBuilder.printNumbers();

		System.out.println("Total first pass clusters: " + totalClusters);
		System.out.println("Filtered out: " + badClusters + " ("
				+ (badClusters * 1.0 / totalClusters) + ")");
		System.out.println("Good clusters: " + goodClusters + " ("
				+ (goodClusters * 1.0 / totalClusters) + ")");

	}

	/**
	 * apply some basic filters to the current cluster
	 * 
	 * 1. singletons (only 1 read support on only one strand) not allowed
	 * 2. if read support is on only 1 strand, it must be at least 3 reads
	 * 
	 * @param words
	 * @return
	 */
	private static boolean goodToWrite(String[] words)
	{
		int psSupport = Integer.parseInt(words[3]);
		int nsSupport = Integer.parseInt(words[4]);
		int other = 0;

		if (psSupport > 0 && nsSupport > 0)
		{
			return true;
		}

		if (psSupport == 0)
		{
			other = nsSupport;
		}

		if (nsSupport == 0)
		{
			other = psSupport;
		}

		return other >= 3;
	}

	private static boolean goodToWriteDuplex(String[] words)
	{
		int psSupport = Integer.parseInt(words[3]);
		int nsSupport = Integer.parseInt(words[4]);

		// true-duplex filtering

		if (psSupport >= 1 && nsSupport >= 1)
		{
			return true;
		}

		return false;
	}

	/**
	 * 
	 * @param mateContigIndex
	 * @param mateStartPosition
	 * @param contigIndex
	 * @param startPosition
	 * @return negative if the first position is smaller
	 *         positive if the first position is larger
	 *         0 if the two positions are equal
	 */
	private static int compareGenomicPositions(int index1, int position1,
			int index2, int position2)
	{
		if (index1 != index2)
		{
			return index1 - index2;
		}

		return position1 - position2;
	}

	private static void writeReadPair(String[] read1, String[] read2,
			BufferedWriter fastq1, BufferedWriter fastq2) throws IOException
	{
		String readName = makeReadName(read1, read2);
		fastq1.write(readName);
		fastq1.write("\n");
		fastq1.write(read1[8]);
		fastq1.write("\n+\n");
		fastq1.write(read1[9]);
		fastq1.write("\n");

		fastq2.write(readName);
		fastq2.write("\n");
		fastq2.write(read2[8]);
		fastq2.write("\n+\n");
		fastq2.write(read2[9]);
		fastq2.write("\n");

		fastq1.flush();
		fastq2.flush();

	}

	private static String makeReadName(String[] read1Info, String[] read2Info)
	{
		StringBuilder builder = new StringBuilder("@Marianas:");

		builder.append(read1Info[2]).append(':').append(read1Info[0])
				.append(':').append(read1Info[1]).append(':')
				.append(read1Info[3]).append(':').append(read1Info[4])
				.append(':').append(read2Info[0]).append(':')
				.append(read2Info[1]).append(':').append(read2Info[3])
				.append(':').append(read2Info[4]);

		return builder.toString();
	}
}
