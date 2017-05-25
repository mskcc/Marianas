/**
 * 
 */
package org.mskcc.marianas.commands;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;

/**
 * @author Juber Patel
 *
 */
public class GetAverageInsertSize
{
	public static final int maxFragmentSize = 600;

	public static void getAverageInsertSize(File bamFile, Interval interval)
			throws IOException
	{
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(bamFile);

		int start = interval.getStart();
		int end = interval.getEnd();

		// adding some padding to catch the reads which are on target or
		// near target
		SAMRecordIterator iterator = reader.query(interval.getContig(),
				start - 50, end + 50, false);

		// get the read count for this interval
		long readCount = 0;
		long readCountWithoutDuplicates = 0;
		long TLEN = 0;
		long TLENWithoutDuplicates = 0;
		long ignored = 0;

		while (iterator.hasNext())
		{
			SAMRecord record = iterator.next();

			// read and its mate must form a proper pair
			if (record.getReadUnmappedFlag() || record.getMateUnmappedFlag())
			{
				continue;
			}

			int i = record.getInferredInsertSize();
			// the second read of the pair, don't count
			if (i <= 0)
			{
				continue;
			}

			if (i > maxFragmentSize)
			{
				ignored++;
				continue;
			}

			readCount++;
			TLEN += i;

			if (!record.getDuplicateReadFlag())
			{
				readCountWithoutDuplicates++;
				TLENWithoutDuplicates += i;
			}

		}

		iterator.close();
		reader.close();

		System.out.println(ignored + " pairs ignored");
		System.out.println(TLEN + "\t" + readCount + "\t" + (TLEN / readCount));
		System.out.println(
				TLENWithoutDuplicates + "\t" + readCountWithoutDuplicates + "\t"
						+ (TLENWithoutDuplicates / readCountWithoutDuplicates));

	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		File bamFile = new File(
				"/Users/patelj1/workspace/Moynahan/FinalBams/1196-2-IGO-05500-AL-21_bc37_5500-AL_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam");
		// ERBB2
		Interval interval = new Interval("17", 37855812, 37884297);

		getAverageInsertSize(bamFile, interval);

	}

}
