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
import java.util.LinkedHashMap;
import java.util.Map;

import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.Frequency;
import org.mskcc.juber.genotype.Genotype;

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

	private static Map<String, String> hotspots;

	private static Map<FreqID, Frequency> afFrequencies;
	private static Map<FreqID, Frequency> countFrequencies;
	private static Map<FreqID, Integer> averageCoverages;
	private static Map<String, Integer> mafColumns;
	private static String mafHeader;

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		// must be same panel, i.e. must have same chrmosome and position on
		// corresponding lines
		File mafFile = new File(args[0]);
		String depthColumnName = args[1];
		String altColumnName = args[2];
		File afFrequenciesFile = new File(args[3]);
		File countFrequenciesFile = new File(args[4]);

		// no more args beyond this point

		System.out.println("Loading Noise Frequencies");
		loadNoiseFrequencies(afFrequenciesFile, countFrequenciesFile);

		System.out.println("Polishing Variants");
		polish(mafFile, depthColumnName, altColumnName);
		System.out.println("Done");
	}

	private static void polish(File mafFile, String depthColumnName,
			String altColumnName) throws IOException
	{
		BufferedReader mafReader = new BufferedReader(new FileReader(mafFile));

		BufferedWriter writer = new BufferedWriter(new FileWriter(
				mafFile.getName().replace(".maf", "-polished.maf")));

		DecimalFormat df = new DecimalFormat("#.####");

		String mafLine = mafReader.readLine();
		processMafHeader(mafLine);
		writer.write(mafHeader + "\n");

		while ((mafLine = mafReader.readLine()) != null)
		{
			String[] tokens = mafLine.split("\t");

			String p = "-";
			String cov = "-";
			String freqs = "-";

			// look at only SNPs and single base deletions
			if (!tokens[mafColumns.get("Variant_Type")].equals("SNP")
					&& !(tokens[mafColumns.get("Variant_Type")].equals("DEL")
							&& tokens[mafColumns.get("Reference_Allele")]
									.length() == 1
							&& tokens[mafColumns.get("Tumor_Seq_Allele2")]
									.equals("-")))
			{
				// no need to get p, cov, freqs
			}
			else
			{
				String chr = tokens[mafColumns.get("Chromosome")];
				int position = Integer
						.parseInt(tokens[mafColumns.get("Start_Position")]);
				char ref = tokens[mafColumns.get("Reference_Allele")].charAt(0);
				char alt = tokens[mafColumns.get("Tumor_Seq_Allele2")]
						.charAt(0);

				if (tokens[mafColumns.get("Variant_Type")].equals("DEL"))
				{
					alt = 'D';
				}

				int depth = Integer
						.parseInt(tokens[mafColumns.get(depthColumnName)]);
				int altCount = Integer
						.parseInt(tokens[mafColumns.get(altColumnName)]);
				double af = (altCount * 1.0) / depth;

				FreqID id = new FreqID(chr, position, ref, alt);
				Frequency afFreq = afFrequencies.get(id);
				Frequency countFreq = countFrequencies.get(id);
				Integer averageCoverage = averageCoverages.get(id);

				if (afFreq == null)
				{
					p = "-";
					cov = "-";
				}
				else
				{
					p = Double.toString(Tester.tTest(afFreq, af));
					cov = Integer.toString(averageCoverage);
				}

				if (countFreq == null)
				{
					freqs = "-";
				}
				else
				{
					freqs = frquencyTableString(countFreq);
				}
			}

			StringBuilder s = new StringBuilder(mafLine);
			s.append("\t").append(cov);
			s.append("\t").append(p);
			s.append("\t").append(freqs);
			s.append("\n");

			writer.write(s.toString());
		}

		writer.close();
		mafReader.close();
	}

	private static void processMafHeader(String header)
	{
		mafColumns = new LinkedHashMap<String, Integer>();

		// Add polishing fields
		mafHeader = header + "\tPolishing_Position_Average_Coverage"
				+ "\tPolishing_P_Value" + "\tFragment_Count->Samples_Map";

		String[] parts = mafHeader.split("\t");

		for (int i = 0; i < parts.length; i++)
		{
			mafColumns.put(parts[i], i);
		}
	}

	public String getMafHeader()
	{
		return mafHeader;
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
