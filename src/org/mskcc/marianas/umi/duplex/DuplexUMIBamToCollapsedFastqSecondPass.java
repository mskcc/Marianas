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
import java.util.List;

import org.mskcc.marianas.util.ClusterCollectionBuilder;
import org.mskcc.marianas.util.StaticResources;
import org.mskcc.marianas.util.Util;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;

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
	 *            args[0] - bam file; args[1] - bed file; args[2] - UMI allowed
	 *            mismatches; args[3] - UMI allowed wobble; args[4] - reference
	 *            fasta; args[5] - output folder
	 * 
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception
	{
		File bamFile = new File(args[0]);
		File bedFile = new File(args[1]);
		int UMIMismatches = Integer.parseInt(args[2]);
		int wobble = Integer.parseInt(args[3]);
		// set the reference fasta
		new StaticResources(new IndexedFastaSequenceFile(new File(args[4])));
		File outputFolder = new File(args[5]);

		// no args after this point

		long start = System.currentTimeMillis();

		// preliminaries
		String fileName = bedFile.getName();
		String intervalsLabel = fileName.substring(0, fileName.indexOf(".bed"));

		// go through the intervals in the given bed file and collect numbers
		List<Interval> intervals = Util.loadIntervals(bedFile);

		File firstPassFile = new File(outputFolder,
				"first-pass.mate-position-sorted.txt");

		System.out.println("Marianas Loeb UMI Second Pass");
		System.out.println("Processing " + firstPassFile.getAbsolutePath()
				+ " to produce fastqs");

		secondPass(bamFile, intervals, UMIMismatches, wobble, outputFolder,
				firstPassFile);

		long end = System.currentTimeMillis();
		System.out.println("Finished processing in " + ((end - start) / 1000)
				+ " seconds.");
	}

	private static void secondPass(File bamFile, List<Interval> intervals,
			int mismatches, int wobble, File outputFolder, File firstPassFile)
			throws Exception
	{
		BufferedReader firstPassReader = new BufferedReader(
				new FileReader(firstPassFile));

		BufferedWriter fastq1 = new BufferedWriter(
				new FileWriter(new File(outputFolder, "collapsed_R1_.fastq")));

		BufferedWriter fastq2 = new BufferedWriter(
				new FileWriter(new File(outputFolder, "collapsed_R2_.fastq")));

		ClusterCollectionBuilder clusterBuilder = new ClusterCollectionBuilder(
				bamFile, wobble, mismatches, false);

		DuplicateReadClusterCollection clusterCollection = clusterBuilder
				.next();
		String line = null;
		String[] words = null;
		int wantedContigIndex = -1;
		int wantedStartPosition = -1;
		int totalClusters = 0;
		int badClusters = 0;
		int goodClusters = 0;

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

			// apply filters
			if (!goodToWrite(words))
			{
				badClusters++;
				continue;
			}

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
						String[] words2 = clusters[i]
								.consensusSequenceInfo(false).split("\t");

						writeReadPair(words, words2, fastq1, fastq2);
						goodClusters++;
						break;
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
