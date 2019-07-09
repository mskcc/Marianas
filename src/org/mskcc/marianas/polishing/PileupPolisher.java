/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.Frequency;

/**
 * @author Juber Patel
 * 
 *         Simple somatic variant caller. Assumes tumor-normal matched pair.
 * 
 *         Currently works on Waltz pileup files. Might switch to bams for more
 *         involved analysis in future.
 *
 */
public class PileupPolisher
{

	private static final double pValueCutoff = 0.0001;

	private static Map<FreqID, Frequency> afFrequencies = null;
	private static Map<FreqID, Frequency> countFrequencies = null;
	private static Map<FreqID, Integer> averageCoverages = null;

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		// must be same panel, i.e. must have same chrmosome and position on
		// corresponding lines
		File pileupFile = new File(args[0]);
		File afFrequenciesFile = new File(args[1]);
		File countFrequenciesFile = new File(args[2]);

		// no more args beyond this point

		System.out.println("Loading Noise Frequencies");
		loadNoiseFrequencies(afFrequenciesFile, countFrequenciesFile);

		System.out.println("Polishing Variants");
		polish(pileupFile);
	}

	private static void polish(File pileupFile) throws IOException
	{
		BufferedReader pileupReader = new BufferedReader(
				new FileReader(pileupFile));

		BufferedWriter polishingWriter = new BufferedWriter(new FileWriter(
				pileupFile.getName().replace("pileup", "polished-pileup")));

		DecimalFormat df = new DecimalFormat("#.####");

		String pileupLine = pileupReader.readLine();

		while ((pileupLine = pileupReader.readLine()) != null)
		{
			String[] tokens = pileupLine.split("\t");

			String chr = tokens[0];
			int position = Integer.parseInt(tokens[1]);
			char ref = tokens[2].charAt(0);

			int[] counts = new int[5];
			counts[0] = Integer.parseInt(tokens[4]);
			counts[1] = Integer.parseInt(tokens[5]);
			counts[2] = Integer.parseInt(tokens[6]);
			counts[3] = Integer.parseInt(tokens[7]);
			counts[4] = Integer.parseInt(tokens[9]);
			int total = 0;
			for (int i = 0; i < counts.length; i++)
			{
				total += counts[i];
			}

			// for each possible SNV at this position
			for (int i = 0; i < counts.length; i++)
			{
				char alt = indexToBase(i);
				double af = (counts[i] * 1.0) / total;

				FreqID id = new FreqID(chr, position, ref, alt);
				Frequency afFreq = afFrequencies.get(id);
				if(afFreq == null || af > 0.20)
				{
					continue;
				}
				
				double p = Tester.MannWhitneyUTest(afFreq, af);

				if (p >= pValueCutoff)
				{
					counts[i] = 0;
				}
			}

			// make new line and write it
			tokens[4] = Integer.toString(counts[0]);
			tokens[5] = Integer.toString(counts[1]);
			tokens[6] = Integer.toString(counts[2]);
			tokens[7] = Integer.toString(counts[3]);
			tokens[9] = Integer.toString(counts[4]);

			StringBuilder s = new StringBuilder(tokens[0]);
			for (int i = 1; i < tokens.length; i++)
			{
				s.append("\t").append(tokens[i]);
			}
			s.append("\n");

			polishingWriter.write(s.toString());

		}

		polishingWriter.close();
		pileupReader.close();
	}

	private static char indexToBase(int i)
	{
		char c;
		switch (i)
		{
		case 0:
			c = 'A';
			break;
		case 1:
			c = 'C';
			break;
		case 2:
			c = 'G';
			break;
		case 3:
			c = 'T';
			break;
		case 4:
			c = 'D';
			break;
		default:
			c = 'X';
		}

		return c;
	}

	/**
	 * load noise frequencies from given file
	 * 
	 * @param afFrequenciesFile
	 * @throws IOException
	 */
	private static void loadNoiseFrequencies(File afFrequenciesFile,
			File countFrequenciesFile) throws IOException
	{

		countFrequencies = new HashMap<FreqID, Frequency>();

		// in case the pileup file does not exist or is not readable
		if (!afFrequenciesFile.canRead())
		{
			System.err.println(
					afFrequenciesFile.getPath() + " is not readable. Exiting.");
			System.exit(1);
		}

		if (!countFrequenciesFile.canRead())
		{
			System.err.println(countFrequenciesFile.getPath()
					+ " is not readable. Exiting.");
			System.exit(1);
		}

		loadAF(afFrequenciesFile);
		loadCount(countFrequenciesFile);

	}

	private static void loadAF(File frequenciesFile)
			throws NumberFormatException, MathIllegalArgumentException,
			IOException
	{
		afFrequencies = new HashMap<FreqID, Frequency>();
		averageCoverages = new HashMap<FreqID, Integer>();

		BufferedReader reader = new BufferedReader(
				new FileReader(frequenciesFile));

		String line = null;
		while ((line = reader.readLine()) != null)
		{
			String[] tokens = line.split("\t");
			FreqID id = new FreqID(tokens[0], Integer.parseInt(tokens[1]),
					tokens[3].charAt(0), tokens[4].charAt(0));

			averageCoverages.put(id, new Integer(tokens[2]));

			Frequency frequencyTable = new Frequency();
			for (int i = 6; i < tokens.length; i += 2)
			{
				double value = Double.parseDouble(tokens[i - 1]);
				int freq = Integer.parseInt(tokens[i]);
				// add value to the frequency table
				frequencyTable.incrementValue(value, freq);
			}

			afFrequencies.put(id, frequencyTable);
		}

		reader.close();
	}

	private static void loadCount(File frequenciesFile)
			throws NumberFormatException, MathIllegalArgumentException,
			IOException
	{
		countFrequencies = new HashMap<FreqID, Frequency>();

		BufferedReader reader = new BufferedReader(
				new FileReader(frequenciesFile));

		String line = null;
		while ((line = reader.readLine()) != null)
		{
			String[] tokens = line.split("\t");
			FreqID id = new FreqID(tokens[0], Integer.parseInt(tokens[1]),
					tokens[3].charAt(0), tokens[4].charAt(0));

			Frequency frequencyTable = new Frequency();
			for (int i = 6; i < tokens.length; i += 2)
			{
				int value = Integer.parseInt(tokens[i - 1]);
				int freq = Integer.parseInt(tokens[i]);
				// add value to the frequency table
				frequencyTable.incrementValue(value, freq);
			}

			countFrequencies.put(id, frequencyTable);
		}

		reader.close();
	}
}
