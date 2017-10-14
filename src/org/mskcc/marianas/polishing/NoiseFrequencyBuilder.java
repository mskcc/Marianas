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
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * @author Juber Patel
 * 
 *         Build per position per substitution error models from the pileup
 *         files in the given directory
 * 
 *         Output: mean and sd of the normal distribution for each such model
 *
 */
public class NoiseFrequencyBuilder
{
	private static final int minCoverage = 50;

	private static Map<String, Map<Integer, PositionNoiseFrequencies>> noiseFrequencies = //
			new LinkedHashMap<String, Map<Integer, PositionNoiseFrequencies>>();

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		File directory = new File(args[0]);

		System.out.println(
				"Collecting position-susbstitution-specific allele fractions");

		// populate models
		for (File file : directory.listFiles())
		{
			if (!file.getName().endsWith("-pileup.txt"))
			{
				continue;
			}

			collectFrequencies(file);
		}

		System.out.println("Saving noise frequencies");
		saveFrequencies();

		System.out.println("Done");

	}

	private static void saveFrequencies() throws IOException
	{
		File modelsFile = new File("noise-models.txt");
		BufferedWriter writer = new BufferedWriter(new FileWriter(modelsFile));

		// DecimalFormat df = new DecimalFormat("#.####");

		for (String chr : noiseFrequencies.keySet())
		{
			Map<Integer, PositionNoiseFrequencies> chrModels = noiseFrequencies
					.get(chr);

			for (Integer position : chrModels.keySet())
			{
				Map<Substitution, NoiseFrequencyCollector> positionNoiseModels = chrModels
						.get(position).getModelsMap();

				for (Substitution substitution : positionNoiseModels.keySet())
				{
					NoiseFrequencyCollector noiseModel = positionNoiseModels
							.get(substitution);

					writer.write(chr + "\t");
					writer.write(position + "\t");
					writer.write(noiseModel.getAverageCoverage() + "\t");
					writer.write(substitution.getRef() + "\t");
					writer.write(substitution.getAlt() + "\t");

					writer.write(noiseModel.getFrequencyString());

					writer.write("\n");
				}
			}
		}

		writer.close();

	}

	private static void collectFrequencies(File file) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));

		String[] words = null;

		for (String line = reader.readLine(); line != null; line = reader
				.readLine())
		{
			words = line.split("\t");

			// apply coverage threshold
			int coverage = Integer.parseInt(words[3]);
			if (coverage < minCoverage)
			{
				continue;
			}

			Map<Integer, PositionNoiseFrequencies> chrFrequencies = noiseFrequencies
					.get(words[0]);

			// add chrModels if doesn't exist already
			if (chrFrequencies == null)
			{
				chrFrequencies = new LinkedHashMap<Integer, PositionNoiseFrequencies>();
				noiseFrequencies.put(words[0], chrFrequencies);
			}

			PositionNoiseFrequencies positionNoiseFrequencies = chrFrequencies
					.get(Integer.parseInt(words[1]));

			// add position noise frequency collectors if doesn't exist already
			if (positionNoiseFrequencies == null)
			{
				positionNoiseFrequencies = new PositionNoiseFrequencies();
				chrFrequencies.put(Integer.parseInt(words[1]),
						positionNoiseFrequencies);
			}

			// iterate over read counts and update models
			int indexToSkip = baseToPileupIndex(words[2].charAt(0));
			for (int i = 4; i < 8; i++)
			{
				if (i == indexToSkip)
				{
					continue;
				}

				int count = Integer.parseInt(words[i]);
				char alt = pileupIndexToBase(i);

				Substitution substitution = Substitution.get(words[2].charAt(0),
						alt);

				NoiseFrequencyCollector frequencyCollector = positionNoiseFrequencies
						.getFrequencyCollectorFor(substitution);
				frequencyCollector.add(count, coverage);
			}
		}

		reader.close();

	}

	private static char pileupIndexToBase(int i)
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

	private static int baseToPileupIndex(char base)
	{
		switch (base)
		{
		case 'A':
			return 4;
		case 'C':
			return 5;
		case 'G':
			return 6;
		case 'T':
			return 7;
		default:
			return -1;
		}
	}

}
