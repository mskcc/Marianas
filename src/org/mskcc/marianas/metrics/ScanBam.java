/**
 * 
 */
package org.mskcc.marianas.metrics;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.mskcc.marianas.util.CustomCaptureException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * @author Juber Patel
 *
 */
public class ScanBam
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

		System.out.println("Scanning " + bamFile.getName());

		long start = System.currentTimeMillis();

		// go through the bam file serially, collect general stats that do not
		// depend on bed files
		// also find regions with average coverage >= coverageThreshould
		scanBam(bamFile, coverageThreshold, geneListFile);

		long end = System.currentTimeMillis();
		System.out.println("Finished processing in " + ((end - start) / 1000)
				+ " seconds.");
	}

	private static void scanBam(File bamFile, int coverageThreshold,
			File geneListFile) throws IOException
	{
		System.out.println("Scanning bam file");

		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(bamFile);
		SAMRecordIterator iterator = reader.iterator();
		// SAMRecordIterator iterator = reader.query("11", 60000, 76000, false);

		long mapped = 0;
		long unmapped = 0;
		long duplicates = 0;
		long totalReads = 0;
		// TODO should we use basic filter instead of these piecemeal tests?
		CoveredRegions coveredRegions = new CoveredRegions(bamFile.getName(),
				coverageThreshold, geneListFile);

		while (iterator.hasNext())
		{
			SAMRecord record = iterator.next();
			totalReads++;

			if (record.getReadUnmappedFlag())
			// not applying quality filter for the time being
			// || record.getMappingQuality() < Constants.minMappingQuality)
			{
				unmapped++;
				continue;
			}

			mapped++;

			if (record.getDuplicateReadFlag())
			{
				duplicates++;
			}
			else
			{
				// not a duplicate read, count towards covered regions
				coveredRegions.recordAlignment(record);
			}
		}

		iterator.close();
		reader.close();

		String bamFileName = bamFile.getName();
		coveredRegions.write();

		BufferedWriter writer = new BufferedWriter(
				new FileWriter(bamFileName + ".stats"));
		writer.write(bamFileName + "\t");
		writer.write(totalReads + "\t");
		writer.write(unmapped + "\t");
		writer.write(mapped + "\t");
		writer.write(duplicates + "\t");
		writer.write(((duplicates * 1.0) / mapped) + "\n");

		writer.close();

	}
}
