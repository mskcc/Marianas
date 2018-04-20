/**
 * 
 */
package org.mskcc.marianas.variantcalling;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;

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
public class VariantCaller
{

	private static Map<String, String> hotspots = null;

	private static Map<NoiseModelID, NoiseModel> noiseModels = null;

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		// must be same panel, i.e. must have same chrmosome and position on
		// corresponding lines
		File tumorPileupFile = new File(args[0]);
		File normalPileupFile = new File(args[1]);
		String sampleName = args[2];
		File hotspotsFile = new File(args[3]);
		File noiseFrequenciesFile = new File(args[4]);

		// no more args beyond this point

		System.out.println("Loading Mutation Hotspots");
		loadHotspots(hotspotsFile);

		System.out.println("Loading Noise Frequencies");
		loadNoiseFrequencies(noiseFrequenciesFile);

		System.out.println("Calling Variants");
		callSNVs(tumorPileupFile, normalPileupFile, sampleName);
	}

	private static void callSNVs(File tumorPileupFile, File normalPileupFile,
			String sampleName) throws IOException
	{
		BufferedReader tumorReader = new BufferedReader(
				new FileReader(tumorPileupFile));
		BufferedReader normalReader = new BufferedReader(
				new FileReader(normalPileupFile));

		BufferedWriter callsWriter = new BufferedWriter(
				new FileWriter(sampleName + "-calls.txt"));
		callsWriter.write(
				"Sample\tChr\tPosition\tRef\tAlt\tAltCount\tAltFraction\tPValue\tTumorTotal\tTumorA\tTumorC\tTumorG\tTumorT"
						+ "\tNormalTotal\tNormalA\tNormalC\tNormalG\tNormalT\tGene\tMutation\tExistingRecord\n");

		DecimalFormat df = new DecimalFormat("#.####");

		String[] tumorTokens = null;
		String[] normalTokens = null;

		for (String tumorLine = tumorReader
				.readLine(), normalLine = normalReader
						.readLine(); tumorLine != null; tumorLine = tumorReader
								.readLine(), normalLine = normalReader
										.readLine())
		{
			tumorTokens = tumorLine.split("\t");
			normalTokens = normalLine.split("\t");

			// tokens: 0-chr, 1-position, 2-ref, 3-total, 4-a, 5-c, 6-g, 7-t,
			// 8-insertions, 9-deletions

			double tumorTotal = (Integer.parseInt(tumorTokens[4])
					+ Integer.parseInt(tumorTokens[5])
					+ Integer.parseInt(tumorTokens[6])
					+ Integer.parseInt(tumorTokens[7])
					+ Integer.parseInt(tumorTokens[9])) * 1.0;
			double normalTotal = (Integer.parseInt(normalTokens[4])
					+ Integer.parseInt(normalTokens[5])
					+ Integer.parseInt(normalTokens[6])
					+ Integer.parseInt(normalTokens[7])
					+ Integer.parseInt(normalTokens[9])) * 1.0;

			if (noCoverage(tumorTotal, normalTotal))
			{
				continue;
			}

			for (int i = 4; i <= 7; i++)
			{
				int tumorCount = Integer.parseInt(tumorTokens[i]);
				int normalCount = Integer.parseInt(normalTokens[i]);

				if (putativeSomatic(tumorTotal, tumorCount, normalTotal,
						normalCount))
				{
					char altBase = toBase(i);
					StringBuilder mutationLine = new StringBuilder(
							sampleName + "\t");
					mutationLine.append(tumorTokens[0]).append("\t");
					mutationLine.append(tumorTokens[1]).append("\t");
					mutationLine.append(tumorTokens[2]).append("\t");
					mutationLine.append(altBase).append("\t");
					mutationLine.append((int) tumorCount).append("\t");

					double af = tumorCount / tumorTotal;
					NoiseModelID id = new NoiseModelID(tumorTokens[0],
							Integer.parseInt(tumorTokens[1]),
							tumorTokens[2].charAt(0), altBase);

					String p = null;
					NoiseModel noiseModel = noiseModels.get(id);
					if (noiseModel == null
							|| noiseModel.getNumObservations() < 10)
					{
						p = "-";
					}
					else
					{
						p = Double.toString(noiseModel.getPValueFor(af));
					}

					mutationLine.append(df.format(af)).append("\t");
					mutationLine.append(p).append("\t");

					mutationLine.append((int) tumorTotal).append("\t");
					mutationLine
							.append(df.format(Integer.parseInt(tumorTokens[4])))
							.append("\t");
					mutationLine
							.append(df.format(Integer.parseInt(tumorTokens[5])))
							.append("\t");
					mutationLine
							.append(df.format(Integer.parseInt(tumorTokens[6])))
							.append("\t");
					mutationLine
							.append(df.format(Integer.parseInt(tumorTokens[7])))
							.append("\t");

					mutationLine.append((int) normalTotal).append("\t");
					mutationLine
							.append(df
									.format(Integer.parseInt(normalTokens[4])))
							.append("\t");
					mutationLine
							.append(df
									.format(Integer.parseInt(normalTokens[5])))
							.append("\t");
					mutationLine
							.append(df
									.format(Integer.parseInt(normalTokens[6])))
							.append("\t");
					mutationLine
							.append(df
									.format(Integer.parseInt(normalTokens[7])))
							.append("\t");

					// add hotspots info
					String hotspotsInfo = hotspots
							.get(tumorTokens[0] + "\t" + tumorTokens[1] + "\t"
									+ tumorTokens[2] + "\t" + altBase);
					if (hotspotsInfo == null)
					{
						hotspotsInfo = "-\t-\t-";
					}

					mutationLine.append(hotspotsInfo);

					callsWriter.write(mutationLine.toString() + "\n");
				}
			}
		}

		callsWriter.close();
		tumorReader.close();
		normalReader.close();

	}

	private static char toBase(int i)
	{
		switch (i)
		{
		case 4:
			return 'A';
		case 5:
			return 'C';
		case 6:
			return 'G';
		case 7:
			return 'T';

		default:
			return 'X';
		}
	}

	/**
	 * This a key method that decides if an SNV is called at a position. All the
	 * filters and calling logic will go here.
	 * 
	 * 
	 * @param tumorTotal
	 * @param tumorCount
	 * @param normalTotal
	 * @param normalCount
	 * @return
	 */
	private static boolean putativeSomatic(double tumorTotal, int tumorCount,
			double normalTotal, int normalCount)
	{
		// TODO think about a nice way to decide putative somatic call
		if (normalCount == 0 && tumorCount > 0)
		{
			return true;
		}
		// else if (normalCount < 3 && tumorCount >= 3)
		// {
		// return true;
		// }
		else if (normalCount / normalTotal < 0.03
				&& tumorCount / tumorTotal > 0.1)
		{
			return true;
		}
		else if ((tumorCount / tumorTotal) > 3 * (normalCount / normalTotal))
		{
			return true;
		}

		return false;
	}

	private static boolean noCoverage(double tumorTotal, double normalTotal)
	{
		int coverageThreshold = 50;

		if (tumorTotal < coverageThreshold || normalTotal < coverageThreshold)
		{
			return true;
		}

		return false;
	}

	/**
	 * load noise frequencies from given file
	 * 
	 * @param noiseFrequenciesFile
	 * @throws IOException
	 */
	private static void loadNoiseFrequencies(File noiseFrequenciesFile)
			throws IOException
	{
		noiseModels = new HashMap<NoiseModelID, NoiseModel>();

		// in case the pileup file does not exist or is not readable
		if (!noiseFrequenciesFile.canRead())
		{
			System.err.println(noiseFrequenciesFile.getPath()
					+ " is not readable. Variants will be called without polishing (noise profiling) !!");
			return;
		}

		BufferedReader reader = new BufferedReader(
				new FileReader(noiseFrequenciesFile));
		for (String line = reader.readLine(); line != null; line = reader
				.readLine())
		{
			String[] tokens = line.split("\t");
			NoiseModelID id = new NoiseModelID(tokens[0],
					Integer.parseInt(tokens[1]), tokens[3].charAt(0),
					tokens[4].charAt(0));

			Frequency frequencyTable = new Frequency();
			for (int i = 6; i < tokens.length; i += 2)
			{
				double value = Double.parseDouble(tokens[i - 1]);
				int freq = Integer.parseInt(tokens[i]);
				// add value to the frequency table
				frequencyTable.incrementValue(value, freq);
			}

			NoiseModel noiseModel = new NoiseModel(id, frequencyTable);

			noiseModels.put(id, noiseModel);
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

	/*
	 * private static int[] getNormalGenotype(String[] normalTokens)
	 * {
	 * // TODO check if this is a good value
	 * double alleleLowerBound = 0.15;
	 * 
	 * int[] readCounts = new int[] { Integer.parseInt(normalTokens[4]),
	 * Integer.parseInt(normalTokens[5]),
	 * Integer.parseInt(normalTokens[6]),
	 * Integer.parseInt(normalTokens[7]) };
	 * 
	 * Arrays.sort(readCounts);
	 * double total = Double.parseDouble(normalTokens[3]);
	 * 
	 * int[] genotype = new int[] { -1, -1 };
	 * int alleleCount = 1;
	 * if(readCounts[1]/total >= alleleLowerBound)
	 * {
	 * alleleCount = 2;
	 * }
	 * 
	 * int genotypeIndex = 0;
	 * for(int i=4; i <= 7 && genotypeIndex < alleleCount ; i++)
	 * {
	 * if(readCounts[genotypeIndex] == )
	 * }
	 * 
	 * 
	 * 
	 * 
	 * 
	 * 
	 * for (int i = 4; i <= 7; i++)
	 * {
	 * double alleleFraction = Integer.parseInt(normalTokens[i]) / total;
	 * if (alleleFraction >= alleleLowerBound)
	 * {
	 * // TODO refine code to handle the unusual situation where more
	 * // than 2 alleles are above threshold
	 * if (genotypeIndex > 1)
	 * {
	 * 
	 * }
	 * else
	 * {
	 * genotype[genotypeIndex] = i;
	 * genotypeIndex++;
	 * }
	 * }
	 * }
	 * }
	 */

}
