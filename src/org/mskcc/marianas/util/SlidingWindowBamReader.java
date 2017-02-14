/**
 * 
 */
package org.mskcc.marianas.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * @author Juber Patel
 * 
 *         An efficient bam reader that presents a sliding window view of the
 *         bam file. It traverses the bam file sequentially and does not use the
 *         random access query mechanism ie query*() methods. The bam file must
 *         be sorted. The window only moves in the forward direction.
 * 
 *         The sliding behavior only happens within a contig. If a new
 *         contig is encountered, a completely new window with the first
 *         windowSize positions' records of the new contig will be created. Use
 *         hasNextOnSameContig() and hasNext() to check end of contig and end of
 *         bam file.
 * 
 *         Also, gaps in coverage are handled in different ways at the ends of
 *         contigs as opposed to between covered regions. A new sliding window
 *         begins where the coverage on the contig begins. The last sliding
 *         window ends where the coverage on the contig ends. But in between
 *         covered regions, the window slides by one position, meaning there can
 *         be many empty windows, with 0 reads at all the positions in the
 *         window.
 *
 */
public class SlidingWindowBamReader
{
	private SamReader reader;
	private SAMRecordIterator iterator;

	/**
	 * current contig of the window
	 */
	private String contig;

	/**
	 * index of the contig in the bam header
	 */
	private int contigIndex;

	/**
	 * current genomic position of the left edge of the window
	 */
	private int startPosition;

	/**
	 * current genomic position of the right edge of the window
	 */
	private int endPosition;

	/**
	 * window records, one list per position
	 */
	private List<SAMRecord>[] recordsForWindow;

	/**
	 * the unconsumed record that will be used first when the window slides the
	 * next time.
	 */
	private SAMRecord nextRecord;

	public final int windowSize;

	public SlidingWindowBamReader(File bamFile, int windowSize)
	{
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		reader = factory.open(bamFile);
		iterator = reader.iterator();

		this.windowSize = windowSize;
		recordsForWindow = new List[windowSize];
		for (int i = 0; i < recordsForWindow.length; i++)
		{
			recordsForWindow[i] = new ArrayList<SAMRecord>();
		}

		// read first record
		nextMappedRecord();
		contig = "0";
		contigIndex = -1;
		startPosition = 0;
		endPosition = 0;

	}

	private void nextMappedRecord()
	{
		while (true)
		{
			if (!iterator.hasNext())
			{
				nextRecord = null;
				return;
			}

			nextRecord = iterator.next();

			if (!nextRecord.getReadUnmappedFlag())
			{
				return;
			}
		}
	}

	/**
	 * slide the window 1 position to the right
	 * 
	 * @return the modified records lists array. Remember it is the same array
	 *         and arrayLists, so be careful with any lingering references.
	 */
	public List<SAMRecord>[] slide()
	{
		// TODO: review this method !!!

		// end of bam file
		if (!hasNext())
		{
			return null;
		}

		// new contig
		if (!hasNextOnSameContig())
		{
			// clear the lists
			for (int i = 0; i < recordsForWindow.length; i++)
			{
				recordsForWindow[i].clear();
			}

			// get new contig and positions
			contig = nextRecord.getReferenceName();
			contigIndex = nextRecord.getReferenceIndex();
			startPosition = nextRecord.getAlignmentStart();
			endPosition = startPosition + recordsForWindow.length - 1;

			while (true)
			{
				// add record
				recordsForWindow[nextRecord.getAlignmentStart() - startPosition]
						.add(nextRecord);

				// read next record
				nextMappedRecord();

				// end of bam file
				if (nextRecord == null)
				{
					return recordsForWindow;
				}

				// position is out of current window, we are done
				if (nextRecord.getAlignmentStart() > endPosition)
				{
					return recordsForWindow;
				}

				// different contig, we are done
				if (!hasNextOnSameContig())
				{
					return recordsForWindow;
				}
			}
		}
		else
		{
			// same contig

			startPosition++;
			endPosition++;

			// shift elements left by one position, reuse the list at [0] at the
			// last position
			List<SAMRecord> t = recordsForWindow[0];
			for (int i = 0; i < recordsForWindow.length - 1; i++)
			{
				recordsForWindow[i] = recordsForWindow[i + 1];
			}

			recordsForWindow[recordsForWindow.length - 1] = t;
			recordsForWindow[recordsForWindow.length - 1].clear();

			// if there are records for the new position, add them
			while (nextRecord.getAlignmentStart() == endPosition)
			{
				recordsForWindow[recordsForWindow.length - 1].add(nextRecord);

				// read next record
				nextMappedRecord();

				// end of bam file
				if (nextRecord == null)
				{
					return recordsForWindow;
				}

				// different contig, we are done
				if (!hasNextOnSameContig())
				{
					return recordsForWindow;
				}
			}

			return recordsForWindow;
		}
	}

	public boolean hasNext()
	{
		if (nextRecord == null)
		{
			return false;
		}

		return true;
	}

	public boolean hasNextOnSameContig()
	{
		if (hasNext() && nextRecord.getReferenceIndex() == contigIndex)
		{
			return true;
		}

		return false;
	}

	public void close() throws IOException
	{
		iterator.close();
		reader.close();
	}

	public String getCurrentWindowContig()
	{
		return contig;
	}

	public int getCurrentWindowContigIndex()
	{
		return contigIndex;
	}

	public int getCurrentWindowStartPosition()
	{
		return startPosition;
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		// File bamFile = new File(
		// "/Users/patelj1/workspace/Marianas/bamFiles/a.bam");
		File bamFile = new File(
				"/Users/patelj1/workspace/Marianas/bamFiles/029-C12-IGO-05500-AW-2_bc06_5500-AW_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam");

		List<SAMRecord>[] records = new ArrayList[7];

		SlidingWindowBamReader reader = new SlidingWindowBamReader(bamFile,
				records.length);

		long start = System.currentTimeMillis();
		long slideCounter = 0;
		while (reader.hasNext())
		{
			records = reader.slide();
			slideCounter++;
			if (slideCounter % 10000000 == 0)
			{
				System.out.println(slideCounter + " slides");
			}

			int i = 0;
			for (i = 0; i < records.length; i++)
			{
				if (records[i].isEmpty())
				{
					break;
				}
			}

			if (i == records.length)
			{
				int a = 5;
			}
		}

		reader.close();

		long time = (System.currentTimeMillis() - start) / 1000;
		System.out.println("Finished in " + time + " seconds.");

	}

}
