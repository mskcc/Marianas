/**
 * 
 */
package org.mskcc.marianas.umi.duplex;

import org.apache.commons.math3.util.FastMath;

/**
 * @author Juber Patel
 * 
 *         Pileup at a specific genomic position (e.g. chr 2 position 3000000)
 * 
 */
public class PositionPileup
{
	private byte refBase;
	// read counts supporting various base calls
	// counts for A, C, G, T and N for this position
	private int[] baseCounts;
	private int[] baseQualities;
	private int deletions;
	private int deletionQualities;

	public PositionPileup()
	{
		// set the values that are not set by default
		// this value should never be seen anywhere
		refBase = '?';
		baseCounts = new int[5];
		baseQualities = new int[5];
	}

	public void reset(byte refBase)
	{
		this.refBase = refBase;
		for (int i = 0; i < baseCounts.length; i++)
		{
			baseCounts[i] = 0;
			baseQualities[i] = 0;

		}

		deletions = 0;
		deletionQualities = 0;
	}

	public void addBase(byte base, byte quality)
	{
		if (base == 'A' || base == 'a')
		{
			baseCounts[0]++;
			baseQualities[0] += quality;
		}
		else if (base == 'C' || base == 'c')
		{
			baseCounts[1]++;
			baseQualities[1] += quality;
		}
		else if (base == 'G' || base == 'g')
		{
			baseCounts[2]++;
			baseQualities[2] += quality;
		}
		else if (base == 'T' || base == 't')
		{
			baseCounts[3]++;
			baseQualities[3] += quality;
		}
		else
		{
			baseCounts[4]++;
			baseQualities[4] += quality;
		}
	}

	public int getCount(byte base)
	{
		if (base == 'A' || base == 'a')
		{
			return baseCounts[0];
		}
		else if (base == 'C' || base == 'c')
		{
			return baseCounts[1];
		}
		else if (base == 'G' || base == 'g')
		{
			return baseCounts[2];
		}
		else if (base == 'T' || base == 't')
		{
			return baseCounts[3];
		}
		else if (base == 'D')
		{
			return deletions;
		}
		else
		{
			return -1;
		}
	}

	public int getQualitySum(byte base)
	{
		if (base == 'A' || base == 'a')
		{
			return baseQualities[0];
		}
		else if (base == 'C' || base == 'c')
		{
			return baseQualities[1];
		}
		else if (base == 'G' || base == 'g')
		{
			return baseQualities[2];
		}
		else if (base == 'T' || base == 't')
		{
			return baseQualities[3];
		}
		else if (base == 'D')
		{
			return deletionQualities;
		}
		else
		{
			return -1;
		}
	}

	/**
	 * 
	 * @return the base that has maximum count at this position. If it's a tie
	 *         or baseCounts[4] has the max count, return -1
	 */
	public byte getMaxCountBase()
	{
		int maxCount = -1;
		for (int i = 0; i < baseCounts.length; i++)
		{
			if (baseCounts[i] > maxCount)
			{
				maxCount = baseCounts[i];
			}
		}

		if (deletions > maxCount)
		{
			return 'D';
		}

		if (baseCounts[4] == maxCount || deletions == maxCount)
		{
			return -1;
		}

		int matches = 0;
		int maxCountIndex = -1;
		byte maxCountBase = -1;
		for (int i = 0; i < baseCounts.length; i++)
		{
			if (baseCounts[i] == maxCount)
			{
				matches++;
				maxCountIndex = i;
			}
		}

		// there has to be a clear winner
		if (matches > 1)
		{
			return -1;
		}

		if (maxCountIndex == 0)
		{
			maxCountBase = 'A';
		}
		else if (maxCountIndex == 1)
		{
			maxCountBase = 'C';
		}
		else if (maxCountIndex == 2)
		{
			maxCountBase = 'G';
		}
		else if (maxCountIndex == 3)
		{
			maxCountBase = 'T';
		}

		return maxCountBase;
	}

	public int getCoverage()
	{
		// TODO double check if deletions should be part of this value
		return baseCounts[0] + baseCounts[1] + baseCounts[2] + baseCounts[3]
				+ deletions;
	}

	public void addDeletion(byte quality)
	{
		deletions++;
		deletionQualities += quality;
	}

	public void merge(PositionPileup other)
	{
		for (int i = 0; i < baseCounts.length; i++)
		{
			baseCounts[i] += other.baseCounts[i];
			baseQualities[i] += other.baseQualities[i];
		}

		deletions += other.deletions;
	}

	public byte getRefBase()
	{
		return refBase;
	}

	public int[] getBaseCounts()
	{
		return baseCounts;
	}

	/**
	 * compute consensus base quality for this strand from base quality scores.
	 * 
	 * @param consensusBase
	 * 
	 * @return
	 */
	public int getConsensusQuality(byte consensusBase)
	{
		int coverage = getCoverage();
		if (coverage == 0)
		{
			return 0;
		}

		// take average, multiply by replication factor
		int qualitySum = getQualitySum(consensusBase);
		int nonBaseQualitySum = baseQualities[0] + baseQualities[1]
				+ baseQualities[2] + baseQualities[3] + deletionQualities
				- qualitySum;
		qualitySum -= nonBaseQualitySum;

		double qualityAverage = (qualitySum * 1.0) / coverage;

		int counts = getCount(consensusBase);
		int nonBaseCounts = coverage - counts;
		counts -= nonBaseCounts;

		double replicationFactor = FastMath.log(2, counts) + 1;
		int quality = (int) (qualityAverage * replicationFactor);
		return quality;

	}
}
