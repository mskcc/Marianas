/**
 * 
 */
package org.mskcc.marianas.umi.duplex.debug;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.mskcc.marianas.umi.duplex.DuplicateReadCluster;
import org.mskcc.marianas.umi.duplex.DuplicateReadClusterCollection;
import org.mskcc.marianas.util.ClusterCollectionBuilder;
import org.mskcc.marianas.util.StaticResources;

import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

/**
 * @author Juber Patel
 *
 */
public class ClusterDetailsAtPosition
{

	/**
	 * @param args
	 *            args[0] - bam file; args[1] - position string; args[2] - UMI
	 *            allowed
	 *            mismatches; args[3] - UMI allowed wobble; args[4] - reference
	 *            fasta;
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception
	{
		File bamFile = new File(args[0]);
		File pileupFile = new File(args[1]);
		String position = args[2];
		int UMIMismatches = Integer.parseInt(args[3]);
		int wobble = Integer.parseInt(args[4]);
		// set the reference fasta
		File refFastaFile = new File(args[5]);
		File refFastaIndexFile = new File(args[5] + ".fai");
		FastaSequenceIndex refFastaIndex = new FastaSequenceIndex(
				refFastaIndexFile);
		IndexedFastaSequenceFile refFasta = new IndexedFastaSequenceFile(
				refFastaFile, refFastaIndex);
		new StaticResources(refFasta, refFastaIndex, pileupFile, position);
		File outputFile = new File(args[6]);

		// no args after this point

		long start = System.currentTimeMillis();

		System.out.println("Marianas Debug Module");
		System.out.println("Position: " + position);
		System.out.println("Processing " + bamFile.getName());

		firstPass(bamFile, UMIMismatches, wobble, outputFile);

		long end = System.currentTimeMillis();
		System.out.println("Finished processing in " + ((end - start) / 1000)
				+ " seconds.");
	}

	private static void firstPass(File bamFile, int mismatches, int wobble,
			File outputFile) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));

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

			recordPositiveStrandClusterCollection(clusterCollection, writer);
		}

		writer.close();
		clusterBuilder.close();
		clusterBuilder.printNumbers();
	}

	private static void recordPositiveStrandClusterCollection(
			DuplicateReadClusterCollection clusterCollection,
			BufferedWriter writer) throws IOException
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

			// TODO decide what else you want to write for first pass
			cluster.collapseMe(writer, true);
		}
	}
}
