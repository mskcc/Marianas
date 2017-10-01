/**
 * 
 */
package org.mskcc.marianas.umi.duplex;

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
	private int insertions;
	private int deletions;

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

		insertions = 0;
		deletions = 0;
	}

	public void addBase(byte base)
	{
		if (base == 'A' || base == 'a')
		{
			baseCounts[0]++;
		}
		else if (base == 'C' || base == 'c')
		{
			baseCounts[1]++;
		}
		else if (base == 'G' || base == 'g')
		{
			baseCounts[2]++;
		}
		else if (base == 'T' || base == 't')
		{
			baseCounts[3]++;
		}
		else
		{
			baseCounts[4]++;
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
		else if (base == 'I')
		{
			return insertions;
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

	/**
	 * 
	 * @return the base that has maximum count at this position. If it's a tie
	 *         or baseCounts[4] has the max count, return -1
	 */
	public byte getMaxCountBase()
	{
		int maxCount = -1;
		for(int i=0; i<baseCounts.length; i++)
		{
			if(baseCounts[i]>maxCount)
			{
				maxCount = baseCounts[i];
			}
		}
		
		if(deletions > maxCount)
		{
			return 'D';
		}
		
		if(baseCounts[4]==maxCount || deletions == maxCount)
		{
			return -1;
		}
		
		int matches = 0;
		int maxCountIndex = -1;
		byte maxCountBase = -1;
		for(int i=0; i<baseCounts.length; i++)
		{
			if(baseCounts[i]==maxCount)
			{
				matches++;
				maxCountIndex = i;
			}
		}
		
		// there has to be a clear winner
		if(matches > 1)
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

	public void addDeletion()
	{
		deletions++;
	}

	public void addInsertion()
	{
		insertions++;
	}

	public String toString()
	{
		StringBuilder builder = new StringBuilder();

		builder.append((char) refBase).append('\t');
		builder.append(getCoverage() + baseCounts[4]).append('\t');
		builder.append(baseCounts[0]).append('\t');
		builder.append(baseCounts[1]).append('\t');
		builder.append(baseCounts[2]).append('\t');
		builder.append(baseCounts[3]).append('\t');
		builder.append(insertions).append('\t');
		builder.append(deletions).append('\t');

		return builder.toString();
	}

	public void merge(PositionPileup other)
	{
		for (int i = 0; i < baseCounts.length; i++)
		{
			baseCounts[i] += other.baseCounts[i];
			baseQualities[i] += other.baseQualities[i];
		}

		insertions += other.insertions;
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
}
