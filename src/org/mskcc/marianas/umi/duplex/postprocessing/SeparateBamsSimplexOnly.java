/**
 * 
 */
package org.mskcc.marianas.umi.duplex.postprocessing;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * @author Juber Patel
 * 
 *         Subset a collapsed bam into multiple bam files: simplex+duplex bam
 *         and pure duplex bam
 *
 */
public class SeparateBamsSimplexOnly
{

	/**
	 * @param args
	 *            path to the collapsed bam
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		File collapsedBam = new File(args[0]);
		File simplexOnlyBam = new File(
				collapsedBam.getName().replace(".bam", "-simplex-only.bam"));
		File unfilterdOnlyBam = new File(
				collapsedBam.getName().replace(".bam", "-unfiltered-only.bam"));

		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(collapsedBam);
		SAMFileHeader header = reader.getFileHeader();
		SAMRecordIterator iterator = reader.iterator();

		SAMFileWriter simplexOnlyWriter = new SAMFileWriterFactory()
				.setCreateIndex(true)
				.makeBAMWriter(header, true, simplexOnlyBam);
		SAMFileWriter unfilteredOnlyWriter = new SAMFileWriterFactory()
				.setCreateIndex(true)
				.makeBAMWriter(header, true, unfilterdOnlyBam);

		while (iterator.hasNext())
		{
			SAMRecord record = iterator.next();

			String[] words = record.getReadName().split(":");
			int ps1 = Integer.parseInt(words[4]);
			int ns1 = Integer.parseInt(words[5]);
			int ps2 = Integer.parseInt(words[8]);
			int ns2 = Integer.parseInt(words[9]);

			// duplex
			if (ps1 >= 1 && ns1 >= 1 && ps2 >= 1 && ns2 >= 1)
			{

			}
			else if (ps1 + ns1 >= 3 && ps2 + ns2 >= 3)
			{
				// simplex
				simplexOnlyWriter.addAlignment(record);
			}
			else
			{
				// unfiltered
				unfilteredOnlyWriter.addAlignment(record);

			}
		}

		simplexOnlyWriter.close();
		unfilteredOnlyWriter.close();
		iterator.close();
		reader.close();
	}

}
