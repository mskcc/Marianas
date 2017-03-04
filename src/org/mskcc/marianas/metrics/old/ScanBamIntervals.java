/**
 * 
 */
package org.mskcc.marianas.metrics.old;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import org.mskcc.marianas.util.Constants;
import org.mskcc.marianas.util.CustomCaptureException;
import org.mskcc.marianas.util.Util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;

/**
 * @author Juber Patel
 *
 */
public class ScanBamIntervals
{
	/**
	 * @param args
	 * @throws IOException
	 * @throws CustomCaptureException
	 */
	public static void main(String[] args)
			throws IOException, CustomCaptureException
	{
		File bamFile = new File(args[0]);
		System.out.println("Processing " + bamFile.getName());

		long start = System.currentTimeMillis();

		for (int i = 1; i < args.length; i++)
		{
			File intervalsFile = new File(args[i]);
			List<Interval> intervals = Util.loadIntervals(intervalsFile);
			String fileName = intervalsFile.getName();
			String type = fileName.substring(0, fileName.indexOf(".bed"));
			processBamAtIntervals(bamFile, intervals, type);
		}

		// SAMRecordIterator iterator = reader.query("11", 60000, 76000, false);

		long end = System.currentTimeMillis();
		System.out.println("Finished processing in " + ((end - start) / 1000)
				+ " seconds.");
	}

	private static void processBamAtIntervals(File bamFile,
			List<Interval> intervals, String type)
					throws CustomCaptureException, IOException
	{
		System.out.println(type);

		DecimalFormat decimalFormat = new DecimalFormat("#.##");

		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(bamFile);
		SAMFileHeader header = reader.getFileHeader();

		String outFile = bamFile.getName() + "." + type;
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

		for (Interval interval : intervals)
		{
			String contigName = Util.guessCorrectContigName(header,
					interval.getContig());

			int start = interval.getStart();
			int end = interval.getEnd();
			int intervalLength = end - start + 1;

			// adding some padding to catch the reads which are on target or
			// near target
			SAMRecordIterator iterator = reader.query(contigName, start, end,
					false);

			int index = interval.getName().indexOf('_');
			String gene = interval.getName().substring(0, index);
			String probe = interval.getName().substring(index + 1);

			// get the read count for this interval
			int readCount = 0;
			int readCountWithoutDuplicates = 0;
			while (iterator.hasNext())
			{
				SAMRecord record = iterator.next();
				if (!record.getReadUnmappedFlag())
				{
					readCount++;

					if (!record.getDuplicateReadFlag())
					{
						readCountWithoutDuplicates++;
					}
				}
			}

			iterator.close();

			// normalize count by interval length
			int coverage = (readCount * Constants.readLength) / intervalLength;
			int coverageWithoutDuplicates = (readCountWithoutDuplicates
					* Constants.readLength) / intervalLength;
			double amplification = 0;

			if (coverageWithoutDuplicates != 0)
			{
				amplification = (coverage * 1.0) / (coverageWithoutDuplicates);
			}
			// write the line for this interval
			writer.write(gene + "\t" + probe + "\t" + readCountWithoutDuplicates
					+ "\t" + readCount + "\t" + intervalLength + "\t"
					+ coverageWithoutDuplicates + "\t" + coverage + "\t"
					+ decimalFormat.format(amplification) + "\n");
		}

		reader.close();
		writer.close();

	}
}
