/**
 * 
 */
package org.mskcc.marianas.util;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.pool2.ObjectPool;
import org.apache.commons.pool2.PooledObjectFactory;
import org.apache.commons.pool2.impl.GenericObjectPool;
import org.apache.commons.pool2.impl.GenericObjectPoolConfig;
import org.mskcc.marianas.umi.loeb.withcorrection.DuplicateReadCluster;
import org.mskcc.marianas.umi.loeb.withcorrection.DuplicateReadClusterCollection;

import htsjdk.samtools.SAMRecord;

/**
 * @author Juber Patel
 * 
 *         A class that builds and presents UMI left clusters at each position.
 *         Left
 *         clusters are the ones that come first while going 5' -> 3' on the
 *         positive strand (and while going 3' -> 5' on the negative strand).
 *         The corresponding right clusters can be pulled using
 *         the read names and UMI's from the left clusters.
 * 
 * 
 *         This class is built on top of SlidingWindowBamReader and gives UMI
 *         clusters at one position at a time. It will take into account wobble
 *         and allowed mismatches while constructing clusters.
 *
 */
public class FivePrimeClusterCollectionBuilder
{
	private int wobble;
	private int mismatches;
	private SlidingWindowBamReader reader;
	private DuplicateReadClusterCollection[] currentClusterWindow;
	private List<SAMRecord>[] currentRecordWindow;

	/**
	 * window index of the position that has been processed most recently. It
	 * indirectly points to the genomic position.
	 */
	private int currentPositionIndex;

	// UMI -> unmapped count map
	private Map<String, Integer> unmapped;
	private long totalUMIPairs;
	private long polyGUMIs;
	private long validFragments;
	private long invalidFragments;
	private ObjectPool<DuplicateReadCluster> clusterPool;
	private Comparator<DuplicateReadCluster> clusterCountComparator;

	public FivePrimeClusterCollectionBuilder(File bamFile, int wobble,
			int mismatches)
	{
		this.wobble = wobble;
		this.mismatches = mismatches;
		this.reader = new SlidingWindowBamReader(bamFile, wobble * 2 + 1);
		this.currentClusterWindow = new DuplicateReadClusterCollection[wobble
				* 2 + 1];
		makeClusterPool();

		for (int i = 0; i < currentClusterWindow.length; i++)
		{
			currentClusterWindow[i] = new DuplicateReadClusterCollection(
					clusterPool);
		}

		// initial value to ensure that sliding starts promptly
		this.currentPositionIndex = reader.windowSize - 1;
		// this.currentPositionIndex = -1;

		this.unmapped = new HashMap<String, Integer>();
		this.clusterCountComparator = new Comparator<DuplicateReadCluster>()
		{
			@Override
			public int compare(DuplicateReadCluster c1, DuplicateReadCluster c2)
			{
				return c2.getTotalCount() - c1.getTotalCount();
			}
		};
	}

	private void makeClusterPool()
	{
		// TODO make proper config !!!

		GenericObjectPoolConfig config = new GenericObjectPoolConfig();

		PooledObjectFactory<DuplicateReadCluster> factory = new DuplicateReadClusterFacory();
		this.clusterPool = new GenericObjectPool<DuplicateReadCluster>(factory,
				config);
	}

	/**
	 * 
	 * @return get the next cluster collection
	 * @throws Exception
	 */
	public DuplicateReadClusterCollection next() throws Exception
	{
		if (needToSlide())
		{
			// the most common case
			if (canSlideOnSameContig())
			{
				// no change in currentPositionIndex, it remains equal to wobble
				currentRecordWindow = reader.slide();

				// basic logging
				if (reader.getCurrentWindowStartPosition() % 10000000 == 0)
				{
					System.out.println("SlidingWindowBamReader: "
							+ reader.getCurrentWindowContig() + ":"
							+ reader.getCurrentWindowStartPosition());
				}

				slideClusterWindow(false);
				finalizeCurrentClusterCollection();
				return currentClusterWindow[currentPositionIndex];
			}
			else
			{
				// can't slide on the same contig. End of contig and possibly
				// end of bam file

				// consume the remaining positions in the current window before
				// going to the next window
				if (currentRecordWindow != null
						&& currentPositionIndex < currentRecordWindow.length
								- 1)
				{
					currentPositionIndex++;
					finalizeCurrentClusterCollection();
					return currentClusterWindow[currentPositionIndex];
				}
				else
				{
					// we have consumed all the positions from the current
					// window
					// so get the first window on the next contig

					// end of bam file
					if (!reader.hasNext())
					{
						return null;
					}

					currentRecordWindow = reader.slide();
					System.out.println("SlidingWindowBamReader: "
							+ reader.getCurrentWindowContig() + ":"
							+ reader.getCurrentWindowStartPosition());
					currentPositionIndex = 0;
					slideClusterWindow(true);
					finalizeCurrentClusterCollection();
					return currentClusterWindow[currentPositionIndex];
				}
			}
		}
		else
		{
			// no need to slide, just process
			// still in the beginning section of the first window on a contig
			currentPositionIndex++;
			finalizeCurrentClusterCollection();
			return currentClusterWindow[currentPositionIndex];
		}

	}

	private boolean needToSlide()
	{
		return (currentPositionIndex >= wobble);
	}

	private boolean canSlideOnSameContig()
	{
		return reader.hasNextOnSameContig();
	}

	public void close() throws IOException
	{
		currentClusterWindow = null;
		currentRecordWindow = null;
		reader.close();
	}

	/**
	 * 
	 * @param rebuildAllPositions
	 *            true: rebuild cluster collections for all the positions in the
	 *            window
	 *            false: shift clusters to the the left by one and build the new
	 *            cluster collection for the last position
	 * @throws Exception
	 */
	private void slideClusterWindow(boolean rebuildAllPositions)
			throws Exception
	{
		if (rebuildAllPositions)
		{
			for (int i = 0; i < currentClusterWindow.length; i++)
			{
				buildClusterCollection(i);
			}
		}
		else
		{
			DuplicateReadClusterCollection t = currentClusterWindow[0];

			// shift left by 1
			for (int i = 0; i < currentClusterWindow.length - 1; i++)
			{
				currentClusterWindow[i] = currentClusterWindow[i + 1];
			}

			// reuse the cluster collection that was at index 0
			currentClusterWindow[currentClusterWindow.length - 1] = t;
			buildClusterCollection(currentClusterWindow.length - 1);
		}

	}

	/**
	 * build the cluster collection for the UMI clusters at the genomic position
	 * indicated by windowIndex in currentRecordWindow and currentClusterWindow.
	 * 
	 * @throws Exception
	 **/
	private void buildClusterCollection(int windowIndex) throws Exception
	{
		// TODO Review this method !!!

		// find the UMI clusters at windowIndex

		List<SAMRecord> records = currentRecordWindow[windowIndex];

		DuplicateReadClusterCollection clusterCollection = currentClusterWindow[windowIndex];
		clusterCollection.prepareFor(reader.getCurrentWindowContig(),
				reader.getCurrentWindowStartPosition() + windowIndex);

		for (SAMRecord record : records)
		{
			String readName = record.getReadName();
			String UMIString = readName
					.substring(readName.lastIndexOf(':') + 1);

			// at least one read of the pair is unmapped
			if (record.getReadUnmappedFlag() || record.getMateUnmappedFlag())
			{
				Integer count = unmapped.get(UMIString);
				if (count == null)
				{
					count = 0;
				}

				count++;
				unmapped.put(UMIString, count);
				continue;
			}

			// this is pass 1
			// only look at the reads mapping on the positive strand
			// that is sufficient to build the left cluster
			//
			if (record.getReadNegativeStrandFlag())
			{
				continue;
			}

			totalUMIPairs++;

			// int maxGs = (UMIString.length() / 4) + 1;
			int maxNucleotides = 10;
			if (Util.poly(UMIString, 'A', maxNucleotides))
			{
				polyGUMIs++;
				continue;
			}

			// only process read pairs that are concordant ie have proper
			// orientation
			int concordance = isGoodAlignment(record);

			// discordant read pair
			if (concordance == 0)
			{
				invalidFragments++;
				continue;
			}

			validFragments++;

			UMIString = getStandardForm(UMIString);

			// positive strand: read1 mapped on positive strand and read2
			// mapped on negative strand
			// negative strand: read1 mapped on negative strand and
			// read2 mapped on positive strand
			boolean positiveStrand = concordance == 1 ? true : false;

			// add the UMI to the cluster collection
			clusterCollection.add(UMIString, record, positiveStrand);
		}
	}

	/**
	 * return the lexicographically smallest form of duplex UMI ie put
	 * lexicographically smaller UMI first
	 * 
	 * @param uMIString
	 * @return
	 */
	private String getStandardForm(String UMIString)
	{
		char UMISeparator = '+';

		// We will rearrange UMIString so that the lexicographically
		// smaller UMI comes first
		String UMI1 = UMIString.substring(0, UMIString.indexOf(UMISeparator));
		String UMI2 = UMIString.substring(UMIString.indexOf(UMISeparator) + 1);

		String UMIString1 = null;
		if (UMI1.compareTo(UMI2) < 0)
		{
			UMIString1 = UMIString;
		}
		else
		{
			UMIString1 = UMI2 + UMISeparator + UMI1;
		}

		return UMIString1;
	}

	/**
	 * 
	 * @param record
	 * @return 0 if the read pair does not have proper relative orientation
	 *         1 if the read is first of the pair ie read1, meaning the fragment
	 *         maps on the positive strand
	 *         2 if the read is not first of the pair ie read2, meaning the
	 *         fragment maps on the negative strand
	 */
	private int isGoodAlignment(SAMRecord record)
	{
		// TODO decide what quality value filter you want to apply
		if (record.getMappingQuality() < 1)
		{
			return 0;
		}

		boolean readStrand = record.getReadNegativeStrandFlag();
		boolean mateStrand = record.getMateNegativeStrandFlag();

		// read and mate mapped on the same strands
		if (!(readStrand ^ mateStrand))
		{
			return 0;
		}

		if (!readStrand && record.getFirstOfPairFlag())
		{
			return 1;
		}
		else
		{
			return 2;
		}
	}

	/**
	 * Finalize the cluster collection at currentPositionIndex by merging
	 * clusters that are within wobble and mismatch criteria
	 * 
	 * @throws Exception
	 */
	private void finalizeCurrentClusterCollection() throws Exception
	{
		// iterate over the wobble span and merge records whose UMIs are within
		// wobble and mismatch distances
		/*
		 * for (int i = currentPositionIndex - wobble; i <= currentPositionIndex
		 * + wobble; i++)
		 * {
		 * // implement ???????
		 * 
		 * // ignore invalid indexes
		 * if (i < 0 || i >= currentRecordWindow.length)
		 * {
		 * continue;
		 * }
		 * 
		 * }
		 */

		// TODO enable wobble processing !!!

		// ignore wobble right now, merge the clusters at the current position
		// that are within mismatches distance
		DuplicateReadClusterCollection clustercollection = currentClusterWindow[currentPositionIndex];

		// remember that this array is a member of the
		// DuplicateReadClusterCollection object, changes will be reflected
		// directly
		DuplicateReadCluster[] sortedClusters = asCountSortedArray(
				clustercollection);

		// process the cluster collection array
		for (int i = 0; i < sortedClusters.length; i++)
		{
			if (sortedClusters[i] == null)
			{
				continue;
			}

			for (int j = i + 1; j < sortedClusters.length; j++)
			{
				if (sortedClusters[j] == null)
				{
					continue;
				}

				if (sameCluster(sortedClusters[i], sortedClusters[j]))
				{
					sortedClusters[i].merge(sortedClusters[j]);
					sortedClusters[j] = null;
				}
			}
		}
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
	private boolean sameCluster(DuplicateReadCluster cluster1,
			DuplicateReadCluster cluster2)
	{
		if (Math.abs(cluster1.getStartPosition()
				- cluster2.getStartPosition()) <= wobble
				&& Util.distance(cluster1.getUMI(),
						cluster2.getUMI()) <= mismatches)
		{
			return true;
		}

		return false;
	}

	private DuplicateReadCluster[] asCountSortedArray(
			DuplicateReadClusterCollection clusters)
	{
		// these are not processed yet!
		DuplicateReadCluster[] clusterArray = clusters.getProcessedClusters();
		Arrays.sort(clusterArray, clusterCountComparator);
		return clusterArray;
	}

	public long getTotalUMIPairs()
	{
		return totalUMIPairs;
	}

	public long getPolyGUMIs()
	{
		return polyGUMIs;
	}

	public long getValidFragments()
	{
		return validFragments;
	}

	public long getInvalidFragments()
	{
		return invalidFragments;
	}
}
