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
public class SeparateBams
{

	/**
	 * @param args
	 *            path to the collapsed bam
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		File collapsedBam = new File(args[0]);
		File simplexBam = new File(
				collapsedBam.getName().replace(".bam", "-simplex.bam"));
		File duplexBam = new File(
				collapsedBam.getName().replace(".bam", "-duplex.bam"));

		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(collapsedBam);
		SAMFileHeader header = reader.getFileHeader();
		SAMRecordIterator iterator = reader.iterator();

		SAMFileWriter simplexWriter = new SAMFileWriterFactory()
				.setCreateIndex(true)
				.makeBAMWriter(header, true, simplexBam);
		SAMFileWriter duplexWriter = new SAMFileWriterFactory()
				.setCreateIndex(true).makeBAMWriter(header, true, duplexBam);

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
				duplexWriter.addAlignment(record);
			}
			else if (ps1 + ns1 >= 3 && ps2 + ns2 >= 3)
			{
				// simplex
				simplexWriter.addAlignment(record);
			}
		}

		simplexWriter.close();
		duplexWriter.close();
		iterator.close();
		reader.close();
	}

}
