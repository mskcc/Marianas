/**
 * 
 */
package org.mskcc.marianas.umi.duplex;

import org.apache.commons.math3.util.FastMath;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.array.TIntArrayList;

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
	private TIntArrayList[] baseQualities;
	private double[] LR;
	private int deletions;
	private int deletionQualities;
	private final int deletionMinConsensusPercent;

	public PositionPileup(int deletionMinConsensusPercent)
	{
		// set the values that are not set by default
		// this value should never be seen anywhere
		refBase = '?';
		baseCounts = new int[5];
		baseQualities = new TIntArrayList[5];
		for (int i = 0; i < baseQualities.length; i++)
		{
			baseQualities[i] = new TIntArrayList();
		}

		LR = new double[4];

		this.deletionMinConsensusPercent = deletionMinConsensusPercent;
	}

	public void reset(byte refBase)
	{
		this.refBase = refBase;
		for (int i = 0; i < baseCounts.length; i++)
		{
			baseCounts[i] = 0;
			baseQualities[i].clear();
		}

		for (int i = 0; i < LR.length; i++)
		{
			LR[i] = 1;
		}

		deletions = 0;
		deletionQualities = 0;
	}

	public void addBase(byte base, byte quality)
	{
		int index = toIndex(base);
		if (index == -1)
		{
			baseCounts[4]++;
			baseQualities[4].add(quality);
		}
		else
		{
			baseCounts[index]++;
			baseQualities[index].add(quality);
		}
	}

	public void addDeletion(byte quality)
	{
		deletions++;
		deletionQualities += quality;
	}

	public int getCount(byte base)
	{
		if (base == 'D')
		{
			return deletions;
		}

		int index = toIndex(base);
		if (index == -1)
		{
			return -1;
		}

		return baseCounts[index];
	}

	public double getLR(byte base)
	{
		int index = toIndex(base);
		if (index == -1)
		{
			return -1;
		}

		return LR[index];
	}
	
	public double[] getLR()
	{
		return LR;
	}

	public TIntArrayList getBaseQualities(byte base)
	{
		int index = toIndex(base);
		if (index == -1)
		{
			return null;
		}

		return baseQualities[index];
	}

	public TIntArrayList[] getBaseQualities()
	{
		return baseQualities;
	}

	public int getQualitySum(byte base)
	{
		if (base == 'D')
		{
			return deletionQualities;
		}

		int index = toIndex(base);
		if (index == -1)
		{
			return -1;
		}

		return baseQualities[index].sum();
	}

	private int toIndex(byte base)
	{
		if (base == 'A' || base == 'a')
		{
			return 0;
		}
		else if (base == 'C' || base == 'c')
		{
			return 1;
		}
		else if (base == 'G' || base == 'g')
		{
			return 2;
		}
		else if (base == 'T' || base == 't')
		{
			return 3;
		}
		else
		{
			return -1;
		}

	}

	public int getQualitySum()
	{
		return baseQualities[0].sum() + baseQualities[1].sum()
				+ baseQualities[2].sum() + baseQualities[3].sum()
				+ deletionQualities;
	}

	public byte getConsensus()
	{
		// count based approach for deletion (maybe at least 60 or 70%)
		// quality based approach for A,C,G,T

		// some threshold for too many non-canonical bases
		if (baseCounts[4] > getCoverage() / 2 || getCoverage() == 0)
		{
			return -1;
		}

		// deletion consensus is count based, it has to be at least
		// deletionMinConsensusPercent
		double deletionConsensusPercent = (deletions * 100.0d) / getCoverage();
		if (deletionConsensusPercent + 0.01 >= deletionMinConsensusPercent)
		{
			return 'D';
		}

		// A, C, G, T consensus is based on likelihoods calculated from error
		// probabilities
		// initialize
		for (int i = 0; i < LR.length; i++)
		{
			LR[i] = 1;
		}

		// accumulate likelihoods
		for (int i = 0; i < baseQualities.length - 1; i++)
		{
			TIntIterator it = baseQualities[i].iterator();
			while (it.hasNext())
			{
				double q = it.next();
				double p = FastMath.pow(10, -q / 10);
				for (int k = 0; k < LR.length; k++)
				{
					if (k == i)
					{
						LR[k] *= (1 - p);
					}
					else
					{
						LR[k] *= (p / 3);
					}
				}
			}
		}

		// compute likelhood ratios
		double t = 0;
		for (int i = 0; i < LR.length; i++)
		{
			t += LR[i];
		}

		for (int i = 0; i < LR.length; i++)
		{
			LR[i] /= t;
		}

		double maxLR = 0;
		int consensusIndex = -1;
		byte consensusBase = -1;
		for (int i = 0; i < LR.length; i++)
		{
			if (LR[i] > maxLR)
			{
				maxLR = LR[i];
				consensusIndex = i;
			}
		}

		if (consensusIndex == 0)
		{
			consensusBase = 'A';
		}
		else if (consensusIndex == 1)
		{
			consensusBase = 'C';
		}
		else if (consensusIndex == 2)
		{
			consensusBase = 'G';
		}
		else if (consensusIndex == 3)
		{
			consensusBase = 'T';
		}

		return consensusBase;

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

	public void merge(PositionPileup other)
	{
		for (int i = 0; i < baseCounts.length; i++)
		{
			baseCounts[i] += other.baseCounts[i];
			baseQualities[i].addAll(other.baseQualities[i]);
		}

		deletions += other.deletions;
		deletionQualities += other.deletionQualities;
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
		int nonBaseQualitySum = getQualitySum() - qualitySum;
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
