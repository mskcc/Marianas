/**
 * 
 */
package org.mskcc.marianas.umi.duplex;

import java.io.BufferedWriter;
import java.io.File;
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
public class DuplexUMIBamToCollapsedFastqFirstPass
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

		// no args[] after this point

		long start = System.currentTimeMillis();

		File firstPassFile = new File("first-pass.txt");
		File altAlleleFile = new File("first-pass-alt-alleles.txt");
		File insertionsFile = new File("first-pass-insertions.txt");

		System.out.println("Marianas " + StaticResources.version);
		System.out.println("First Pass");
		System.out.println("Processing " + bamFile.getName());

		firstPass(bamFile, minMappingQuality, minBaseQuality, UMIMismatches,
				wobble, minConsensusPercent, firstPassFile, altAlleleFile,
				insertionsFile);

		long end = System.currentTimeMillis();
		System.out.println("Finished processing in " + ((end - start) / 1000)
				+ " seconds.");
	}

	private static void firstPass(File bamFile, int minMappingQuality,
			int minBaseQuality, int mismatches, int wobble,
			int minConsensusPercent, File firstPassFile, File altAlleleFile,
			File insertionsFile) throws Exception
	{

		BufferedWriter firstPassWriter = new BufferedWriter(
				new FileWriter(firstPassFile));
		BufferedWriter altAlleleWriter = new BufferedWriter(
				new FileWriter(altAlleleFile));
		BufferedWriter insertionsWriter = new BufferedWriter(
				new FileWriter(insertionsFile));

		ClusterCollectionBuilder clusterBuilder = new ClusterCollectionBuilder(
				bamFile, minMappingQuality, minBaseQuality, mismatches, wobble,
				minConsensusPercent, true);

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
					firstPassWriter, altAlleleWriter, insertionsWriter);
		}

		clusterBuilder.close();
		firstPassWriter.close();
		altAlleleWriter.close();
		insertionsWriter.close();

		clusterBuilder.printNumbers();
	}

	private static void recordPositiveStrandClusterCollection(
			DuplicateReadClusterCollection clusterCollection,
			BufferedWriter firstPassWriter, BufferedWriter altAlleleWriter,
			BufferedWriter insertionsWriter) throws IOException
	{
		// TODO write all the needed information !!!

		String consensusSequenceInfo = null;
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

			try
			{
				consensusSequenceInfo = cluster.collapseMe(altAlleleWriter,
						insertionsWriter, true);
			}
			catch (Exception e)
			{
				System.err.println(
						"Problem creating consensus sequence info for cluster:");
				System.err.println(cluster.getContig() + ":"
						+ cluster.getStartPosition() + ":" + cluster.getUMI());
				e.printStackTrace();
				continue;
				// System.exit(1);
			}

			// there is meaningful consensus sequence info
			if (consensusSequenceInfo != null)
			{
				firstPassWriter.write(consensusSequenceInfo + "\n");
			}
		}

		firstPassWriter.flush();
		altAlleleWriter.flush();
		insertionsWriter.flush();
	}
}
