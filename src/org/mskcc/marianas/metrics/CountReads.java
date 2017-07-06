/**
 * 
 */
package org.mskcc.marianas.metrics;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.mskcc.juber.intervals.IntervalNameMap;
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
public class CountReads
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
		int coverageThreshold = Integer.parseInt(args[1]);
		File geneListFile = new File(args[2]);
		File intervalsFile = new File(args[3]);

		System.out.println("Scanning " + bamFile.getName());

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
		// create a ReadCounts object
		ReadCounts readCounts = new ReadCounts(bamFile.getName(), sampleID,
				intervalsLabel);
		// create a coveredRegions object
		CoveredRegions coveredRegions = new CoveredRegions(bamFile.getName(),
				coverageThreshold, geneListFile);

		// load intervals
		List<Interval> intervals = Util.loadIntervals(intervalsFile);
		IntervalNameMap intervalNameMap = toIntervalNameMap(intervals);

		// go through the bam file serially, collect general stats that do not
		// depend on bed files
		// also find regions with average coverage >= coverageThreshould
		scanBam(readCounts, coveredRegions, bamFile, intervalNameMap);

		// go through the intervals in the given bed file and collect numbers
		// processBamAtIntervals(readCounts, bamFile, intervals);

		coveredRegions.write();
		readCounts.write();

		long end = System.currentTimeMillis();
		System.out.println("Finished processing in " + ((end - start) / 1000)
				+ " seconds.");
	}

	private static void scanBam(ReadCounts readCounts,
			CoveredRegions coveredRegions, File bamFile,
			IntervalNameMap intervalNameMap) throws IOException
	{
		System.out.println("Scanning entire " + bamFile.getName());

		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(bamFile);
		SAMRecordIterator iterator = reader.iterator();
		// SAMRecordIterator iterator = reader.query("11", 60000, 76000, false);

		while (iterator.hasNext())
		{
			SAMRecord record = iterator.next();
			readCounts.totalReads++;

			if (record.getReadUnmappedFlag())
			// not applying quality filter for the time being
			// || record.getMappingQuality() < Constants.minMappingQuality)
			{
				readCounts.unmappedReads++;
				continue;
			}

			readCounts.totalMappedReads++;

			// check if on target
			List<String> intersecting = intervalNameMap.getIntersecting(
					record.getContig(), record.getAlignmentStart(),
					record.getAlignmentEnd());

			if (!intersecting.isEmpty())
			{
				readCounts.totalTargetReads++;
			}

			if (record.getDuplicateReadFlag())
			{
				readCounts.duplicateMappedReads++;

			}
			else
			{
				readCounts.uniqueMappedReads++;

				if (!intersecting.isEmpty())
				{
					readCounts.uniqueTargetReads++;
				}

				// not a duplicate read, count towards covered regions
				coveredRegions.recordAlignment(record);
			}
		}

		iterator.close();
		reader.close();
	}

	private static IntervalNameMap toIntervalNameMap(List<Interval> intervals)
	{
		// build the interval name map
		IntervalNameMap intervalNameMap = new IntervalNameMap();
		for (Interval interval : intervals)
		{
			// shorten the interval so that we only count on target read when
			// there is minReadOverlap
			int minReadOverlap = 10;
			int start = interval.getStart() + minReadOverlap;
			int end = interval.getEnd() - minReadOverlap;
			// interval is too short for further shortening
			if (start > end)
			{
				start = end = (interval.getStart() + interval.getEnd()) / 2;
			}

			intervalNameMap.add(interval.getContig(), start, end,
					interval.getName());

		}

		return intervalNameMap;
	}
}
