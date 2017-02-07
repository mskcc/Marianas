/**
 * 
 */
package org.mskcc.marianas.umi.simplex;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.mskcc.marianas.umi.loeb.DuplicateReadClusterCollection;
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
 *         Generate UMI metrics from a bam file
 *
 */
public class ProcessSimplexUMIBamOld
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
		File intervalsFile = new File(args[0]);
		File bamFile = new File(args[1]);
		int UMILength = Integer.parseInt(args[2]);

		long start = System.currentTimeMillis();

		// preliminaries
		String fileName = intervalsFile.getName();
		String intervalsLabel = fileName.substring(0, fileName.indexOf(".bed"));
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(bamFile);
		SAMFileHeader header = reader.getFileHeader();
		reader.close();
		// It is assumed that there is only one read group in the bam file!!!
		String sampleID = header.getReadGroups().get(0).getSample();

		// go through the intervals in the given bed file and collect numbers
		List<Interval> intervals = Util.loadIntervals(intervalsFile);

		System.out.println("Processing " + bamFile.getName() + " at "
				+ intervalsLabel + " intervals...");

		processBamAtIntervals(bamFile, intervals, UMILength, sampleID);

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
	 * @param sampleID
	 * @throws CustomCaptureException
	 * @throws IOException
	 */
	private static void processBamAtIntervals(File bamFile,
			List<Interval> intervals, int UMILength, String sampleID)
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
		DuplicateReadClusterCollection clusters = null;

		for (Interval interval : intervals)
		{
			String contigName = Util.guessCorrectContigName(header,
					interval.getContig());
			int start = interval.getStart();
			int end = interval.getEnd();

			SAMRecordIterator iterator = reader.query(contigName, start, end,
					false);

			// int previousAlignmentStart = 0;
			while (iterator.hasNext())
			{
				SAMRecord record = iterator.next();

				int alignmentStart = record.getAlignmentStart();
				String alignmentContig = record.getContig();

				// only process those alignments whose start is in the current
				// interval
				if (alignmentStart < start || alignmentStart > end)
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

				// if(record.isSecondaryOrSupplementary() ||
				// record.getNotPrimaryAlignmentFlag())
				// {
				// continue;
				// }

				// int alignmentStart = record.getAlignmentStart();
				// if (alignmentStart < previousAlignmentStart)
				// {
				// System.out.println("PROBLEM!! " + alignmentStart + " after "
				// + previousAlignmentStart);
				// }

				// previousAlignmentStart = alignmentStart;

				// ignore read2 alignments, since read2 start positions are not
				// informative and read2s should be clustered by corresponding
				// read1 start position and UMI
				if (!record.getFirstOfPairFlag())
				{
					continue;
				}

				// current read belongs to the current cluster collection
				if (clusters != null && clusters.contig.equals(alignmentContig)
						&& clusters.startPosition == alignmentStart)
				{
					clusters.add(UMI);
				}
				else
				{
					// write current cluster collection and create a new one
					// corresponding to current record
					if (clusters != null)
					{
						clusters.write(sampleID, summaryWriter);
						// prep for reuse
						clusters.prepFor(alignmentContig, alignmentStart);

					}
					else
					{
						// first time, initialize clusters
						clusters = new DuplicateReadClusterCollection(
								alignmentContig, alignmentStart);
					}
				}
			}

			iterator.close();
		}

		summaryWriter.close();
		reader.close();
	}

}
