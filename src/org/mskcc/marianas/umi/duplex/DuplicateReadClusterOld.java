/**
 * 
 */
package org.mskcc.marianas.umi.duplex;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

/**
 * @author Juber Patel
 *
 */
public class DuplicateReadClusterOld
{
	private String contig;
	private int startPosition;
	private String UMI;
	private int psReadCount;
	private int nsReadCount;
	private short[][] psBaseCounts;
	private short[][] nsBaseCounts;
	private int psLastValidIndex;
	private int nsLastValidIndex;
	private byte[] psConsensus;
	private byte[] nsConsensus;
	private byte[] consensus;

	public DuplicateReadClusterOld()
	{
		// TODO remove the second arbitrary constant
		// TODO use object pooling to avoid repeated creation of these arrays??
		// (Apache Commons Pool)
		this.psBaseCounts = new short[4][200];
		this.nsBaseCounts = new short[4][200];
		this.psConsensus = new byte[200];
		this.nsConsensus = new byte[200];
		this.consensus = new byte[200];
		this.psLastValidIndex = -1;
		this.nsLastValidIndex = -1;
	}

	public void prepareFor(String contig, int startPosition, String UMI)
	{
		this.contig = contig;
		this.startPosition = startPosition;
		this.UMI = UMI;
		this.psReadCount = 0;
		this.nsReadCount = 0;
		this.psLastValidIndex = -1;
		this.nsLastValidIndex = -1;

		// clear the 2d arrays
		for (int j = 0; j < psBaseCounts[0].length; j++)
		{
			psBaseCounts[0][j] = 0;
		}

		System.arraycopy(psBaseCounts[0], 0, nsBaseCounts[0], 0,
				psBaseCounts[0].length);

		for (int i = 1; i < psBaseCounts.length; i++)
		{
			System.arraycopy(psBaseCounts[0], 0, psBaseCounts[i], 0,
					psBaseCounts[0].length);
			System.arraycopy(psBaseCounts[0], 0, nsBaseCounts[i], 0,
					psBaseCounts[0].length);
		}

		// clear consensus arrays
		System.arraycopy(psBaseCounts[0], 0, psConsensus, 0,
				psBaseCounts[0].length);
		System.arraycopy(psBaseCounts[0], 0, nsConsensus, 0,
				psBaseCounts[0].length);
		System.arraycopy(psBaseCounts[0], 0, consensus, 0,
				psBaseCounts[0].length);

	}

	public void add(byte[] read, boolean positiveStrand)
	{
		short[][] counts = null;
		int readLastIndex = read.length - 1;

		if (positiveStrand)
		{
			psReadCount++;
			psLastValidIndex = readLastIndex > psLastValidIndex ? readLastIndex
					: psLastValidIndex;
			counts = psBaseCounts;

		}
		else
		{
			nsReadCount++;
			nsLastValidIndex = readLastIndex > nsLastValidIndex ? readLastIndex
					: nsLastValidIndex;
			counts = nsBaseCounts;
		}

		// increment the counts in the correct count array
		for (int i = 0; i <= readLastIndex; i++)
		{
			byte base = read[i];

			if (base == 'A')
			{
				counts[0][i]++;
			}
			else if (base == 'C')
			{
				counts[1][i]++;
			}
			else if (base == 'G')
			{
				counts[2][i]++;
			}
			else if (base == 'T')
			{
				counts[3][i]++;
			}
		}
	}

	public char[] consensusRead()
	{
		return null;
	}

	public void merge(DuplicateReadClusterOld other)
	{
		psReadCount += other.psReadCount;
		nsReadCount += other.nsReadCount;

		// add positive strand values
		for (int i = 0; i < psBaseCounts.length; i++)
		{
			for (int j = 0; j <= other.psLastValidIndex; j++)
			{
				psBaseCounts[i][j] += other.psBaseCounts[i][j];
			}
		}

		// add negative strand values
		for (int i = 0; i < nsBaseCounts.length; i++)
		{
			for (int j = 0; j <= other.nsLastValidIndex; j++)
			{
				nsBaseCounts[i][j] += other.nsBaseCounts[i][j];
			}
		}

		// adjust last valid indexes
		if (other.psLastValidIndex > psLastValidIndex)
		{
			psLastValidIndex = other.psLastValidIndex;
		}

		if (other.nsLastValidIndex > nsLastValidIndex)
		{
			nsLastValidIndex = other.nsLastValidIndex;
		}
	}

	public String consensus()
	{
		// calculate per strand consensus
		strandConsensus(psBaseCounts, psLastValidIndex, psConsensus);
		strandConsensus(nsBaseCounts, nsLastValidIndex, nsConsensus);

		// if a position is present in both strands

		// use RegionPileup, PositionPileup, Genotype....

		// Please use DuplicateReadCluster class
		return null;

	}

	/**
	 * compute the strand consensus sequence. Right now, just choosing the base
	 * with highest count at each position.
	 * 
	 * @param baseCounts
	 * @param lastValidIndex
	 * @param strandConsensus
	 */
	private void strandConsensus(short[][] baseCounts, int lastValidIndex,
			byte[] strandConsensus)
	{
		int index = 0;
		for (int i = 0; i <= lastValidIndex; i++)
		{
			index = baseCounts[1][i] > baseCounts[index][i] ? 1 : index;
			index = baseCounts[2][i] > baseCounts[index][i] ? 2 : index;
			index = baseCounts[3][i] > baseCounts[index][i] ? 3 : index;

			if (index == 0)
			{
				strandConsensus[i] = 'A';
			}
			else if (index == 1)
			{
				strandConsensus[i] = 'C';
			}
			else if (index == 2)
			{
				strandConsensus[i] = 'G';
			}
			else if (index == 3)
			{
				strandConsensus[i] = 'T';
			}
		}
	}

	public int getStartPosition()
	{
		return startPosition;
	}

	public String getUMI()
	{
		return UMI;
	}

	public int getTotalCount()
	{
		return psReadCount + nsReadCount;
	}
}
