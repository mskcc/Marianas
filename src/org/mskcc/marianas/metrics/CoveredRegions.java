/**
 * 
 */
package org.mskcc.marianas.metrics;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.mskcc.juber.intervals.IntervalNameMap;
import org.mskcc.marianas.util.Constants;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;

/**
 * @author Juber Patel
 * 
 *         List of regions that have average coverage >= coverage threshold
 *
 */
public class CoveredRegions
{
	private String bamFileName;
	private IntervalNameMap exonMap;
	private int coverageThreshold;
	private Map<Interval, Integer> intervals;
	private String contig;
	private int start;
	private int end;
	private int reads;

	public CoveredRegions(String bamFileName, int coverageThreshold,
			File geneListFile) throws NumberFormatException, IOException
	{
		this.bamFileName = bamFileName;
		this.coverageThreshold = coverageThreshold;
		this.intervals = new LinkedHashMap<Interval, Integer>();

		// build the gene map
		this.exonMap = new IntervalNameMap();
		BufferedReader reader = new BufferedReader(
				new FileReader(geneListFile));
		String line = null;
		while ((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");
			int start = Integer.parseInt(parts[1]);
			int end = Integer.parseInt(parts[2]);
			exonMap.add(parts[0], start, end, parts[4]);
		}

		reader.close();
	}

	public void recordAlignment(SAMRecord record)
	{
		// very first record
		if (contig == null)
		{
			contig = record.getContig();
			start = record.getStart();
			end = record.getEnd();
			reads = 1;
		}
		// there is a gap between current and previous record
		// so record previous region, start a new region
		else if (record == null || !record.getContig().equals(contig)
				|| end < record.getStart())
		{
			int coverage = (reads * Constants.readLength) / (end - start + 1);

			if (coverage >= coverageThreshold)
			{
				Interval interval = new Interval(contig, start, end);
				intervals.put(interval, coverage);
			}

			// last region, end signal
			if (record == null)
			{
				return;
			}

			contig = record.getContig();
			start = record.getStart();
			end = record.getEnd();
			reads = 1;
		}
		else
		{
			reads++;
			int newEnd = record.getEnd();
			if (newEnd > end)
			{
				end = newEnd;
			}
		}
	}

	public void write() throws IOException
	{
		// record the last region
		recordAlignment(null);

		BufferedWriter writer = new BufferedWriter(
				new FileWriter(bamFileName + ".covered-regions"));

		for (Interval interval : intervals.keySet())
		{
			int coverage = intervals.get(interval);

			String contig = interval.getContig();
			int start = interval.getStart();
			int end = interval.getEnd();
			List<String> exons = exonMap.getIntersecting(contig, start, end);

			writer.write(contig + "\t");
			writer.write(start + "\t");
			writer.write(end + "\t");
			writer.write((end - start + 1) + "\t");
			writer.write(coverage + "\t");

			// get the names of exons overlapping the region
			String exonNames = "";
			if (!exons.isEmpty())
			{
				exonNames = exons.get(0);
				for (int i = 1; i < exons.size(); i++)
				{
					exonNames = exonNames + "," + exons.get(i);
				}
			}

			writer.write(exonNames + "\n");
		}

		writer.close();
	}

}
