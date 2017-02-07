/**
 * 
 */
package org.mskcc.marianas.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Splitter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Interval;

/**
 * @author Juber Patel
 * 
 *         Static Utility Methods
 */
public class Util
{

	public static List<Interval> loadIntervals(File intervalsFile)
			throws NumberFormatException, IOException
	{
		Splitter tabSplitter = Splitter.on('\t');
		List<Interval> list = new ArrayList<Interval>();

		BufferedReader reader = new BufferedReader(
				new FileReader(intervalsFile));

		String line = null;
		while ((line = reader.readLine()) != null)
		{
			if (line.startsWith("#"))
				continue;

			List<String> parts = tabSplitter.splitToList(line);
			boolean negative = true;
			if (parts.get(3).trim().equals("+"))
			{
				negative = false;
			}

			int start = Integer.parseInt(parts.get(1).trim());
			int end = Integer.parseInt(parts.get(2).trim());

			list.add(new Interval(parts.get(0), start, end, negative,
					parts.get(4).trim()));
		}

		reader.close();

		return list;

	}

	/**
	 * eg. can we find either "1" or "chr1" or "Chr1"?
	 * 
	 * @param header
	 * @param contig
	 * @return
	 * @throws Exception
	 */
	public static String guessCorrectContigName(SAMFileHeader header,
			String contigName) throws CustomCaptureException
	{
		int index = header.getSequenceIndex(contigName);

		if (index != -1)
			return contigName;

		String newContigName = null;

		if (contigName.startsWith("chr") || contigName.startsWith("Chr"))
		{
			newContigName = contigName.substring(3);
			index = header.getSequenceIndex(newContigName);

			if (index != -1)
			{
				return newContigName;
			}
			else
			{
				throw new CustomCaptureException(
						"Chromosome Not Found: " + contigName);
			}
		}
		else
		{
			newContigName = "chr" + contigName;
			index = header.getSequenceIndex(newContigName);

			if (index != -1)
			{
				return newContigName;
			}
			else
			{
				newContigName = "Chr" + contigName.substring(3);
				index = header.getSequenceIndex(newContigName);

				if (index != -1)
				{
					return newContigName;
				}

				throw new CustomCaptureException(
						"Chromosome Not Found: " + contigName);
			}

		}
	}

	/**
	 * find the number of places where the two given strings differ from each
	 * other
	 * 
	 * @param seq1
	 * @param seq2
	 * @return
	 */
	public static int distance(String seq1, String seq2)
	{
		int count = 0;
		for (int i = 0; i < seq1.length(); i++)
		{
			if (seq1.charAt(i) != seq2.charAt(i))
			{
				count++;
			}
		}

		return count;
	}

	public static boolean poly(String UMI, char nucleotide, int maxNucleotides)
	{
		int count = 0;
		for (int i = 0; i < UMI.length(); i++)
		{
			if (UMI.charAt(i) == nucleotide)
			{
				count++;
			}
		}

		if (count > maxNucleotides)
		{
			return true;
		}

		return false;
	}

	public static String reverseComplement(String seq)
	{
		StringBuilder s = new StringBuilder(seq.length());

		for (int j = seq.length() - 1; j >= 0; j--)
		{
			char c = seq.charAt(j);

			if (c == 'A')
			{
				s.append('T');
			}
			else if (c == 'C')
			{
				s.append('G');
			}
			else if (c == 'G')
			{
				s.append('C');
			}
			else if (c == 'T')
			{
				s.append('A');
			}
			else
			{
				s.append(c);
			}
		}

		return s.toString();
	}
}
