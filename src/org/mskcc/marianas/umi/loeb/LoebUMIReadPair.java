/**
 * 
 */
package org.mskcc.marianas.umi.loeb;

import org.mskcc.marianas.util.Util;

/**
 * @author Juber Patel
 *
 */
public class LoebUMIReadPair
{

	private UMIAndRead read1;
	private UMIAndRead read2;

	public LoebUMIReadPair(String seq1, String qual1, String seq2,
			String qual2, int UMILength, String constantRegion)
	{
		read1 = process(seq1, qual1, UMILength, constantRegion);
		read2 = process(seq2, qual2, UMILength, constantRegion);

	}

	/**
	 * determine if the read has valid UMI and separate the UMI and constant
	 * region
	 * from the read.
	 * 
	 * Right now, applying a strict definition of UMI presence: constant region
	 * must match identically, starting from position UMILength.
	 * 
	 * @param seq
	 * @param qual
	 * @param UMILength
	 * @param constantRegion
	 * @return
	 */
	private LoebUMIReadPair.UMIAndRead process(String seq, String qual,
			int UMILength, String constantRegion)
	{
		// check if there is the constant region at the correct position
		String region = seq.substring(UMILength,
				UMILength + constantRegion.length());
		if (Util.distance(region, constantRegion) <= 1)
		{
			UMIAndRead read = new UMIAndRead(
					seq.substring(UMILength + constantRegion.length()),
					qual.substring(UMILength + constantRegion.length()),
					seq.substring(0, UMILength), qual.substring(0, UMILength));
			return read;
		}

		return null;

	}

	public String compositeUMI()
	{
		return read1.UMI + "+" + read2.UMI;
	}

	public String compositeUMIQuals()
	{
		return read1.UMIQual + "+" + read2.UMIQual;
	}

	public boolean read1HasUMI()
	{
		return read1 != null;
	}

	public boolean read2HasUMI()
	{
		return read2 != null;
	}

	private class UMIAndRead
	{
		private String seq;
		private String qual;
		private String UMI;
		private String UMIQual;

		public UMIAndRead(String seq, String qual, String UMI, String UMIQual)
		{
			this.seq = seq;
			this.qual = qual;
			this.UMI = UMI;
			this.UMIQual = UMIQual;

		}

	}

	public String seq1()
	{
		return read1.seq;
	}

	public String qual1()
	{
		return read1.qual;
	}

	public String seq2()
	{
		return read2.seq;
	}

	public String qual2()
	{
		return read2.qual;
	}

}
