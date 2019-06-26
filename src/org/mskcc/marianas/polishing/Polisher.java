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
import java.util.Iterator;
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
public class Polisher
{

	private static Map<String, String> hotspots = null;

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
		File mafFile = new File(args[0]);
		File afFrequenciesFile = new File(args[1]);
		File countFrequenciesFile = new File(args[2]);

		// no more args beyond this point

		System.out.println("Loading Noise Frequencies");
		loadNoiseFrequencies(afFrequenciesFile, countFrequenciesFile);

		System.out.println("Polishing Variants");
		polish(mafFile);
	}

	private static void polish(File mafFile) throws IOException
	{
		BufferedReader mafReader = new BufferedReader(new FileReader(mafFile));

		BufferedWriter polishingWriter = new BufferedWriter(new FileWriter(
				mafFile.getName().replace(".maf", "") + "-polishing.txt"));

		DecimalFormat df = new DecimalFormat("#.####");

		String mafLine = mafReader.readLine();

		while ((mafLine = mafReader.readLine()) != null)
		{
			String[] tokens = mafLine.split("\t");

			// look at only SNPs and single base deletions
			if (!tokens[9].equals("SNP") && !(tokens[9].equals("DEL")
					&& tokens[10].length() == 1 && tokens[12].equals("-")))
			{
				continue;
			}

			String chr = tokens[4];
			int position = Integer.parseInt(tokens[5]);
			char ref = tokens[10].charAt(0);
			char alt = tokens[12].charAt(0);

			if (tokens[9].equals("DEL"))
			{
				alt = 'D';
			}

			String sample = tokens[15];
			int depth = Integer.parseInt(tokens[144]);
			int altCount = Integer.parseInt(tokens[145]);
			double af = (altCount * 1.0) / depth;

			FreqID id = new FreqID(chr, position, ref, alt);
			Frequency afFreq = afFrequencies.get(id);
			Frequency countFreq = countFrequencies.get(id);
			Integer averageCoverage = averageCoverages.get(id);

			String p1 = null;
			String p2 = null;
			String p3 = null;
			String p4 = null;
			String cov = null;
			String freqs = null;

			if (afFreq == null)
			{
				p1 = p2 = "-";
				cov = "-";
			}
			else
			{
				p1 = Double.toString(Tester.tTest(afFreq, af));
				p2 = Double.toString(Tester.MannWhitneyUTest(afFreq, af));
				cov = Integer.toString(averageCoverage);
			}

			if (countFreq == null)
			{
				p3 = p4 = "-";
				freqs = "-";
			}
			else
			{
				p3 = Double.toString(Tester.tTest(countFreq, af));
				p4 = Double.toString(Tester.MannWhitneyUTest(countFreq, af));
				freqs = frquencyTableString(countFreq);
			}

			StringBuilder s = new StringBuilder();
			s.append(id.toString());
			s.append("\t").append(cov);
			s.append("\t").append(sample);
			s.append("\t").append(altCount);
			s.append("\t").append(depth);
			s.append("\t").append(af);
			s.append("\t").append(p1);
			s.append("\t").append(p2);
			s.append("\t").append(p3);
			s.append("\t").append(p4);
			s.append("\t").append(freqs);
			s.append("\n");

			polishingWriter.write(s.toString());
		}

		polishingWriter.close();
		mafReader.close();
	}

	private static String frquencyTableString(Frequency countFreq)
	{
		StringBuilder s = new StringBuilder();
		Iterator it = countFreq.valuesIterator();
		while (it.hasNext())
		{
			Long value = (Long) it.next();
			long freq = (long) countFreq.getCount(value);
			s.append(value).append("->").append(freq).append(";");
		}

		return s.toString();
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

	/**
	 * load the hotspots from the given file
	 * 
	 * @param hotspotsFile
	 * @throws IOException
	 */
	private static void loadHotspots(File hotspotsFile) throws IOException
	{
		hotspots = new HashMap<String, String>();

		// in case the pileup file does not exist or is not readable
		if (!hotspotsFile.canRead())
		{
			System.err.println(hotspotsFile.getPath()
					+ " is not readable. Will use reference as genotype for all positions.");
			return;
		}

		BufferedReader reader = new BufferedReader(
				new FileReader(hotspotsFile));
		for (String line = reader.readLine(); line != null; line = reader
				.readLine())
		{
			String[] tokens = line.split("\t");
			String value = null;
			if (tokens.length == 8)
			{
				value = tokens[0] + "\t" + tokens[6] + "\t" + tokens[7];
			}
			else
			{
				value = tokens[0] + "\t" + tokens[6] + "\t-";
			}

			hotspots.put(tokens[1] + "\t" + tokens[2] + "\t" + tokens[4] + "\t"
					+ tokens[5], value);
		}

		reader.close();
	}

}
