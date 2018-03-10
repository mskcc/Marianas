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
import org.mskcc.juber.util.Util;
import org.mskcc.marianas.umi.duplex.DuplicateReadCluster;
import org.mskcc.marianas.umi.duplex.DuplicateReadClusterCollection;

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
public class ClusterCollectionBuilder
{
	private int wobble;
	private int mismatches;
	private int minConsensusPercent;
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
	private long validFragments;
	private long invalidFragments;
	private long readsOnPositiveStrand;
	private long readsOnNegativeStrand;
	private long fragmentsFromPositiveStrand;
	private long fragmentsFromNegativeStrand;
	private ObjectPool<DuplicateReadCluster> clusterPool;
	private Comparator<DuplicateReadCluster> clusterCountComparator;
	private boolean positiveStrand;
	private boolean debugMode;

	/**
	 * 
	 * @param bamFile
	 * @param wobble
	 * @param mismatches
	 * @param minConsensusPercent
	 * @param positiveStrand
	 *            the read strand (as opposed to fragment strand). We process
	 *            the reads mapping mapping 5' -> 3' on one strand at a time
	 */
	public ClusterCollectionBuilder(File bamFile, int wobble, int mismatches,
			int minConsensusPercent, boolean positiveStrand)
	{
		this.wobble = wobble;
		this.mismatches = mismatches;
		this.minConsensusPercent = minConsensusPercent;
		this.positiveStrand = positiveStrand;
		this.debugMode = StaticResources.getPositionOfInterest() != null;
		this.reader = new SlidingWindowBamReader(bamFile, wobble * 2 + 1);
		this.currentClusterWindow = new DuplicateReadClusterCollection[wobble
				* 2 + 1];
		// TODO uncomment this if you want to use object pooling. Not using
		// right now.
		// makeClusterPool();

		for (int i = 0; i < currentClusterWindow.length; i++)
		{
			currentClusterWindow[i] = new DuplicateReadClusterCollection(
					clusterPool, minConsensusPercent);
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

		PooledObjectFactory<DuplicateReadCluster> factory = new DuplicateReadClusterFacory(
				minConsensusPercent);
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
				reader.getCurrentWindowContigIndex(),
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

			// only look at the reads mapping on the specified strand
			// unless you are running in debug mode. Then look at both strands
			if (!debugMode
					&& !(positiveStrand ^ record.getReadNegativeStrandFlag()))
			{
				continue;
			}

			totalUMIPairs++;

			// TODO think what mapping quality conditions to impose on read and
			// mate
			// TODO think what other conditions could invalidate a fragment

			// silencing this part to allow detection of SVs
			// discordant, both reads map on the same strand OR
			// read maps at multiple places
			if (record.getReadNegativeStrandFlag() == record
					.getMateNegativeStrandFlag()
					|| record.getMappingQuality() < 1)
			{
				invalidFragments++;
				continue;
			}

			validFragments++;

			// record read strand
			if (positiveStrand)
			{
				readsOnPositiveStrand++;
			}
			else
			{
				readsOnNegativeStrand++;
			}

			// the strand of the fragment, not the read!
			// positive strand: read1 mapped on positive strand and read2
			// mapped on negative strand
			// negative strand: read1 mapped on negative strand and
			// read2 mapped on positive strand
			boolean fragmentPositiveStrand = !(positiveStrand
					^ record.getFirstOfPairFlag());

			// record fragment strand, only do this for read1 to avoid double
			// counting
			if (fragmentPositiveStrand)
			{
				fragmentsFromPositiveStrand++;
			}
			else
			{
				fragmentsFromNegativeStrand++;
			}

			UMIString = getStandardForm(UMIString);

			// add the UMI to the cluster collection
			clusterCollection.add(UMIString, record, fragmentPositiveStrand);
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

	public void printNumbers()
	{
		System.out.println("Total UMI Pairs: " + totalUMIPairs);
		System.out.println("Invalid Fragments: " + invalidFragments + " ("
				+ ((1.0 * invalidFragments) / totalUMIPairs) + ")");
		System.out.println("Valid Fragments: " + validFragments);

		System.out
				.println("Reads on positive strand: " + readsOnPositiveStrand);
		System.out
				.println("Reads on negative strand: " + readsOnNegativeStrand);
		System.out.println("Fragments from positive strand: "
				+ fragmentsFromPositiveStrand);
		System.out.println("Fragments from negative strand: "
				+ fragmentsFromNegativeStrand);
	}

	public long getTotalUMIPairs()
	{
		return totalUMIPairs;
	}

	public long getValidFragments()
	{
		return validFragments;
	}

	public long getInvalidFragments()
	{
		return invalidFragments;
	}

	public long getReadsOnPositiveStrand()
	{
		return readsOnPositiveStrand;
	}

	public long getReadsOnNegativeStrand()
	{
		return readsOnNegativeStrand;
	}

	public long getFragmentsFromPositiveStrand()
	{
		return fragmentsFromPositiveStrand;
	}

	public long getFragmentsFromNegativeStrand()
	{
		return fragmentsFromNegativeStrand;
	}
}
