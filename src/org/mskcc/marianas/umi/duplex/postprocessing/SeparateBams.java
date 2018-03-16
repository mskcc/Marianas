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
		File simplexDuplexBam = new File(
				collapsedBam.getName().replace(".bam", "-simplex-duplex.bam"));
		File duplexBam = new File(
				collapsedBam.getName().replace(".bam", "-duplex.bam"));

		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(collapsedBam);
		SAMFileHeader header = reader.getFileHeader();
		SAMRecordIterator iterator = reader.iterator();

		SAMFileWriter simplexDuplexWriter = new SAMFileWriterFactory()
				.setCreateIndex(true)
				.makeBAMWriter(header, true, simplexDuplexBam);
		SAMFileWriter duplexWriter = new SAMFileWriterFactory()
				.setCreateIndex(true).makeBAMWriter(header, true, duplexBam);

		while (iterator.hasNext())
		{
			SAMRecord record = iterator.next();

			String[] words = record.getReadName().split(":");
			int psSupport = Integer.parseInt(words[4]);
			int nsSupport = Integer.parseInt(words[5]);

			// duplex
			if (psSupport >= 1 && nsSupport >= 1)
			{
				simplexDuplexWriter.addAlignment(record);
				duplexWriter.addAlignment(record);
			}
			else if (psSupport + nsSupport >= 3)
			{
				// simplex
				simplexDuplexWriter.addAlignment(record);
			}
		}

		simplexDuplexWriter.close();
		duplexWriter.close();
		iterator.close();
		reader.close();
	}

}
