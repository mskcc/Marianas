/**
 * 
 */
package org.mskcc.marianas.umi.simplex;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.mskcc.juber.util.CustomCaptureException;
import org.mskcc.juber.util.Util;
import org.mskcc.marianas.umi.duplex.DuplicateReadCluster;
import org.mskcc.marianas.umi.duplex.DuplicateReadClusterCollection;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;

/**
 * @author Juber Patel
 * 
 *         Generate UMI metrics from a bam file
 *
 */
public class ProcessSimplexUMIBam
{

	/**
	 * 
	 * @param
	 */
	/**
	 * @param args
	 *            args[0] - bam file; args[1] - bed file; args[2] - UMI length
	 * @throws IOException
	 * @throws CustomCaptureException
	 */
	public static void main(String[] args)
			throws IOException, CustomCaptureException
	{
		File bamFile = new File(args[0]);
		File bedFile = new File(args[1]);
		int UMILength = Integer.parseInt(args[2]);
		int UMIMismatches = Integer.parseInt(args[3]);
		int wobble = Integer.parseInt(args[4]);

		long start = System.currentTimeMillis();

		// preliminaries
		String fileName = bedFile.getName();
		String intervalsLabel = fileName.substring(0, fileName.indexOf(".bed"));
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(bamFile);
		SAMFileHeader header = reader.getFileHeader();
		reader.close();
		// It is assumed that there is only one read group in the bam file!!!
		String sampleID = header.getReadGroups().get(0).getSample();

		// go through the intervals in the given bed file and collect numbers
		List<Interval> intervals = Util.loadIntervals(bedFile);

		System.out.println("Processing " + bamFile.getName() + " at "
				+ intervalsLabel + " intervals...");

		processBamAtIntervals(bamFile, intervals, UMILength, UMIMismatches,
				wobble, sampleID);

		long end = System.currentTimeMillis();
		System.out.println("Finished processing in " + ((end - start) / 1000)
				+ " seconds.");
	}

	/**
	 * Strategy for generating consensus sequences:
	 * 
	 * A molecule is identified by position AND UMI. read1's are clustered by
	 * these 2 things. read2's are clustered by these same 2 things from read1.
	 * No attributes from read2's are used to cluster read2's. That means read1
	 * and read2 clusters must be built at the same time, fetching read2 for
	 * each read1, using a separate reader. This will be quite slow but seems to
	 * be the best way.
	 * 
	 * 
	 * @param bamFile
	 * @param intervals
	 * @param UMILength
	 * @param wobble
	 * @param uMIMismatches
	 * @param sampleID
	 * @throws CustomCaptureException
	 * @throws IOException
	 */
	private static void processBamAtIntervals(File bamFile,
			List<Interval> intervals, int UMILength, int UMIMismatches,
			int wobble, String sampleID)
			throws CustomCaptureException, IOException
	{
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(bamFile);
		SAMFileHeader header = reader.getFileHeader();

		File summaryFile = new File(bamFile.getName() + ".umi-summary");
		BufferedWriter summaryWriter = new BufferedWriter(
				new FileWriter(summaryFile));

		// UMI -> unmapped count map
		Map<String, Integer> unmapped = new HashMap<String, Integer>();
		long totalUMIs = 0;
		long polyGUMIs = 0;
		int maxGs = (UMILength / 4) + 1;

		for (Interval interval : intervals)
		{
			String contigName = Util.guessCorrectContigName(header,
					interval.getContig());
			int intervalStart = interval.getStart();
			int intervalEnd = interval.getEnd();

			System.out.println("reading records for " + interval);
			// create and initialize cluster collection for each position in the
			// current interval
			DuplicateReadClusterCollection[] clustersByPosition = new DuplicateReadClusterCollection[intervalEnd
					- intervalStart + 1];
			for (int i = 0; i < clustersByPosition.length; i++)
			{
				clustersByPosition[i] = new DuplicateReadClusterCollection(
						contigName, intervalStart + i);
			}

			// iterate
			SAMRecordIterator iterator = reader.query(contigName, intervalStart,
					intervalEnd, false);
			while (iterator.hasNext())
			{
				SAMRecord record = iterator.next();

				int alignmentStart = record.getAlignmentStart();

				// only process those alignments whose start is in the current
				// interval
				if (alignmentStart < intervalStart
						|| alignmentStart > intervalEnd)
				{
					continue;
				}

				String readName = record.getReadName();
				String UMI = readName.substring(readName.length() - UMILength);

				// unmapped read
				if (record.getReadUnmappedFlag())
				{
					Integer count = unmapped.get(UMI);
					if (count == null)
					{
						count = 0;
					}

					count++;
					unmapped.put(UMI, count);
					continue;
				}

				// ignore read2 alignments, since read2 start positions are not
				// informative and read2s should be clustered by corresponding
				// read1 start position and UMI
				if (!record.getFirstOfPairFlag())
				{
					continue;
				}

				totalUMIs++;
				if (Util.polyG(UMI, maxGs))
				{
					polyGUMIs++;
					continue;
				}

				// add the UMI to the correct cluster collection
				clustersByPosition[alignmentStart - intervalStart].add(UMI);
			}

			iterator.close();

			processIntervalUMIs(clustersByPosition, UMIMismatches, wobble,
					sampleID, summaryWriter);
		}

		summaryWriter.close();
		reader.close();

		System.out.println("Poly-G: " + (1.0 * polyGUMIs) / totalUMIs);
	}

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
	private static void processIntervalUMIs(
			DuplicateReadClusterCollection[] clustersByPosition,
			int UMIMismatches, int wobble, String sampleID,
			BufferedWriter summaryWriter) throws IOException
	{

		System.out.println("sorting by cluster count...");
		// probably take another tack: sort all clusters by count and then for
		// each cluster find other clusters that are within wobble and mismatch
		// distance
		DuplicateReadCluster[] clusters = sortClustersByCount(
				clustersByPosition);

		System.out.println("merging...");
		// process the cluster collection array
		for (int i = 0; i < clusters.length; i++)
		{
			if (clusters[i].startPosition == 0)
			{
				continue;
			}

			for (int j = i + 1; j < clusters.length; j++)
			{
				if (clusters[j].startPosition == 0)
				{
					continue;
				}

				if (sameCluster(clusters[i], clusters[j], wobble,
						UMIMismatches))
				{
					clusters[i].count += clusters[j].count;
					clusters[j].startPosition = 0;
					clusters[j].count = 0;
				}
			}
		}

		System.out.println("sorting by position...");
		// sort by position and count
		sortClustersByPositionAndCount(clusters);

		System.out.println("writing...");
		// write clusters
		for (int i = 0; i < clusters.length; i++)
		{
			DuplicateReadCluster cluster = clusters[i];
			if (cluster.startPosition == 0)
			{
				continue;
			}

			summaryWriter.write(sampleID + "\t" + cluster.contig + "\t"
					+ cluster.startPosition + "\t" + cluster.UMI + "\t"
					+ cluster.count + "\n");
		}

		summaryWriter.flush();
	}

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
	private static boolean sameCluster(DuplicateReadCluster cluster1,
			DuplicateReadCluster cluster2, int wobble, int UMIMismatches)
	{
		if (Math.abs(cluster1.startPosition - cluster2.startPosition) <= wobble
				&& Util.distance(cluster1.UMI, cluster2.UMI) <= UMIMismatches)
		{
			return true;
		}

		return false;
	}

	private static DuplicateReadCluster[] sortClustersByCount(
			DuplicateReadClusterCollection[] clustersByPosition)
	{
		List<DuplicateReadCluster> clusters = new ArrayList<DuplicateReadCluster>();

		for (int i = 0; i < clustersByPosition.length; i++)
		{
			clusters.addAll(clustersByPosition[i].clusters.values());
		}

		clusters.sort(new Comparator<DuplicateReadCluster>()
		{
			@Override
			public int compare(DuplicateReadCluster c1, DuplicateReadCluster c2)
			{
				return c2.count - c1.count;
			}
		});

		return clusters.toArray(new DuplicateReadCluster[0]);
	}

	private static void sortClustersByPositionAndCount(
			DuplicateReadCluster[] clusters)
	{
		Arrays.sort(clusters, new Comparator<DuplicateReadCluster>()
		{
			@Override
			public int compare(DuplicateReadCluster c1, DuplicateReadCluster c2)
			{
				int diff = c1.startPosition - c2.startPosition;
				if (diff != 0)
				{
					return diff;
				}
				else
				{
					return c2.count - c1.count;
				}
			}
		});
	}

}
