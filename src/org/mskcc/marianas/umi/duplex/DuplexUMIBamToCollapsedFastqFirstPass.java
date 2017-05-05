/**
 * 
 */
package org.mskcc.marianas.umi.duplex;

import java.io.BufferedWriter;
import java.io.File;
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
public class DuplexUMIBamToCollapsedFastqFirstPass
{
	/**
	 * @param args
	 *            args[0] - bam file; args[1] - pileup file; args[2] - UMI
	 *            allowed
	 *            mismatches; args[3] - UMI allowed wobble; args[4] - reference
	 *            fasta; args[5] - R1 fastq name; args[6] - output folder
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception
	{
		File bamFile = new File(args[0]);
		File pileupFile = new File(args[1]);
		int UMIMismatches = Integer.parseInt(args[2]);
		int wobble = Integer.parseInt(args[3]);
		// set the reference fasta
		new StaticResources(new IndexedFastaSequenceFile(new File(args[4])),
				pileupFile, null);
		File outputFolder = new File(args[5]);

		// no args after this point

		long start = System.currentTimeMillis();

		File firstPassFile = new File(outputFolder, "first-pass.txt");

		System.out.println("Marianas Loeb UMI First Pass");
		System.out.println("Processing " + bamFile.getName());

		firstPass(bamFile, UMIMismatches, wobble, firstPassFile);

		long end = System.currentTimeMillis();
		System.out.println("Finished processing in " + ((end - start) / 1000)
				+ " seconds.");
	}

	private static void firstPass(File bamFile, int mismatches, int wobble,
			File firstPassFile) throws Exception
	{

		BufferedWriter firstPassWriter = new BufferedWriter(
				new FileWriter(firstPassFile));

		ClusterCollectionBuilder clusterBuilder = new ClusterCollectionBuilder(
				bamFile, wobble, mismatches, true);

		DuplicateReadClusterCollection clusterCollection = null;

		// iterate
		while (true)
		{
			clusterCollection = clusterBuilder.next();

			// end of bam file
			if (clusterCollection == null)
			{
				break;
			}

			recordPositiveStrandClusterCollection(clusterCollection,
					firstPassWriter);
		}

		clusterBuilder.close();
		firstPassWriter.close();

		clusterBuilder.printNumbers();
	}

	private static void recordPositiveStrandClusterCollection(
			DuplicateReadClusterCollection clusterCollection,
			BufferedWriter firstPassWriter) throws IOException
	{
		// TODO write all the needed information !!!

		// write clusters
		DuplicateReadCluster[] processedClusters = clusterCollection
				.getProcessedClusters();

		for (int i = 0; i < processedClusters.length; i++)
		{
			DuplicateReadCluster cluster = processedClusters[i];
			if (cluster == null)
			{
				continue;
			}

			firstPassWriter
					.write(cluster.consensusSequenceInfo(null, true) + "\n");
		}

		firstPassWriter.flush();
	}

	//////////////////////////
	////////////////////////////
	//////////////////////////////
	/////////////////////////////

	// old code
	//////////////////////////////////////////

	/*
	 * public void old()
	 * {
	 * while (reader.hasNext())
	 * {
	 * String readName = record.getReadName();
	 * String UMIString = readName
	 * .substring(readName.lastIndexOf(':') + 1);
	 * 
	 * // at least one read of the pair is unmapped
	 * if (record.getReadUnmappedFlag() || record.getMateUnmappedFlag())
	 * {
	 * Integer count = unmapped.get(UMIString);
	 * if (count == null)
	 * {
	 * count = 0;
	 * }
	 * 
	 * count++;
	 * unmapped.put(UMIString, count);
	 * continue;
	 * }
	 * 
	 * // only look at the reads mapping on the positive strand, that
	 * // gives us sufficient information
	 * if (record.getReadNegativeStrandFlag())
	 * {
	 * continue;
	 * }
	 * 
	 * totalUMIPairs++;
	 * 
	 * // only process read pairs that are concordant ie have proper
	 * // orientation
	 * int concordance = getConcordance(record);
	 * 
	 * // int maxGs = (UMIString.length() / 4) + 1;
	 * int maxNucleotides = 10;
	 * if (Util.poly(UMIString, 'A', maxNucleotides))
	 * {
	 * polyGUMIs++;
	 * continue;
	 * }
	 * 
	 * // discordant read pair
	 * if (concordance == 0)
	 * {
	 * invalidFragments++;
	 * continue;
	 * }
	 * 
	 * validFragments++;
	 * 
	 * UMIString = getStandardForm(UMIString);
	 * 
	 * boolean positiveStrand;
	 * 
	 * // positive strand ie read1 mapped on positive strand and read2
	 * // mapped on negative strand
	 * if (concordance == 1)
	 * {
	 * positiveStrand = true;
	 * }
	 * else
	 * {
	 * // negative strand ie read1 mapped on negative strand and
	 * // read2 mapped on positive strand
	 * positiveStrand = false;
	 * }
	 * 
	 * // add the UMI to the correct cluster collection
	 * clustersByPosition[alignmentStart - intervalStart].add(UMIString,
	 * positiveStrand);
	 * }
	 * 
	 * iterator.close();
	 * 
	 * processIntervalUMIs(clustersByPosition, UMIMismatches, wobble, sampleID,
	 * firstPassWriter);
	 * 
	 * firstPassWriter.close();
	 * reader.close();
	 * 
	 * System.out.println("Total UMI Pairs: " + totalUMIPairs);
	 * System.out.println("Poly-G Fragments: " + polyGUMIs + " ("
	 * + ((1.0 * polyGUMIs) / totalUMIPairs) + ")");
	 * System.out.println(
	 * "Non-Poly-G but Invalid Fragments: " + invalidFragments + " ("
	 * + ((1.0 * invalidFragments) / totalUMIPairs) + ")");
	 * System.out.println("Valid Fragments: " + validFragments);
	 * 
	 * }
	 */

	/**
	 * 
	 * @param clustersByPosition
	 * @param uMIMismatches
	 *            number of mismatches allowed in UMIs
	 * @param wobble
	 *            number of positions to consider to look for UMI matches
	 * @param sampleID
	 * @param summaryWriter
	 * @throws IOException
	 */
	/*
	 * private static void processIntervalUMIs(
	 * DuplicateReadClusterCollection[] clustersByPosition,
	 * int UMIMismatches, int wobble, String sampleID,
	 * BufferedWriter summaryWriter) throws IOException
	 * {
	 * 
	 * System.out.println("sorting by cluster count...");
	 * 
	 * // probably take another tack: sort all clusters by count and then for
	 * // each cluster find other clusters that are within wobble and mismatch
	 * // distance
	 * DuplicateReadCluster[] clusters = sortClustersByCount(
	 * clustersByPosition);
	 * 
	 * System.out.println("merging...");
	 * // process the cluster collection array
	 * for (int i = 0; i < clusters.length; i++)
	 * {
	 * if (clusters[i].startPosition == 0)
	 * {
	 * continue;
	 * }
	 * 
	 * for (int j = i + 1; j < clusters.length; j++)
	 * {
	 * if (clusters[j].startPosition == 0)
	 * {
	 * continue;
	 * }
	 * 
	 * if (sameCluster(clusters[i], clusters[j], wobble,
	 * UMIMismatches))
	 * {
	 * clusters[i].psReadCount += clusters[j].psReadCount;
	 * clusters[i].nsReadCount += clusters[j].nsReadCount;
	 * clusters[j].startPosition = 0;
	 * clusters[j].psReadCount = 0;
	 * clusters[j].nsReadCount = 0;
	 * }
	 * }
	 * }
	 * 
	 * System.out.println("sorting by position...");
	 * // sort by position and count
	 * sortClustersByPositionAndCount(clusters);
	 * 
	 * System.out.println("writing...");
	 * 
	 * // write clusters
	 * for (int i = 0; i < clusters.length; i++)
	 * {
	 * DuplicateReadCluster cluster = clusters[i];
	 * if (cluster.startPosition == 0)
	 * {
	 * continue;
	 * }
	 * 
	 * summaryWriter.write(sampleID + "\t" + cluster.contig + "\t"
	 * + cluster.startPosition + "\t" + cluster.UMI + "\t"
	 * + cluster.psReadCount + "\t" + cluster.nsReadCount + "\n");
	 * }
	 * 
	 * summaryWriter.flush();
	 * }
	 */

	/**
	 * determine whether the 2 clusters should be considered as same given the
	 * allowed position wobble and UMI base mismatches
	 * 
	 * @param cluster1
	 * @param cluster2
	 * @param wobble
	 * @param uMIMismatches
	 * @return
	 */
	/*
	 * private static boolean sameCluster(DuplicateReadCluster cluster1,
	 * DuplicateReadCluster cluster2, int wobble, int UMIMismatches)
	 * {
	 * if (Math.abs(cluster1.startPosition - cluster2.startPosition) <= wobble
	 * && Util.distance(cluster1.UMI, cluster2.UMI) <= UMIMismatches)
	 * {
	 * return true;
	 * }
	 * 
	 * return false;
	 * }
	 * 
	 * private static DuplicateReadCluster[] sortClustersByCount(
	 * DuplicateReadClusterCollection[] clustersByPosition)
	 * {
	 * List<DuplicateReadCluster> clusters = new
	 * ArrayList<DuplicateReadCluster>();
	 * 
	 * for (int i = 0; i < clustersByPosition.length; i++)
	 * {
	 * clusters.addAll(clustersByPosition[i].clusters.values());
	 * }
	 * 
	 * clusters.sort(new Comparator<DuplicateReadCluster>()
	 * {
	 * 
	 * @Override
	 * public int compare(DuplicateReadCluster c1, DuplicateReadCluster c2)
	 * {
	 * return (c2.psReadCount + c2.nsReadCount)
	 * - (c1.psReadCount + c1.nsReadCount);
	 * }
	 * });
	 * 
	 * return clusters.toArray(new DuplicateReadCluster[0]);
	 * }
	 * 
	 * private static void sortClustersByPositionAndCount(
	 * DuplicateReadCluster[] clusters)
	 * {
	 * Arrays.sort(clusters, new Comparator<DuplicateReadCluster>()
	 * {
	 * 
	 * @Override
	 * public int compare(DuplicateReadCluster c1, DuplicateReadCluster c2)
	 * {
	 * int diff = c1.startPosition - c2.startPosition;
	 * if (diff != 0)
	 * {
	 * return diff;
	 * }
	 * else
	 * {
	 * return (c2.psReadCount + c2.nsReadCount)
	 * - (c1.psReadCount + c1.nsReadCount);
	 * }
	 * }
	 * });
	 * }
	 */

}
