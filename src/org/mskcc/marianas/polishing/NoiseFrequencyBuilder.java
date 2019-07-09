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
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

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
				"Collecting position-susbstitution-specific allele fractions and allele counts");

		// collect af and fragment count frequencies
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
		File afFile = new File("af-frequencies.txt");
		File countFile = new File("count-frequencies.txt");

		BufferedWriter afWriter = new BufferedWriter(new FileWriter(afFile));
		BufferedWriter countWriter = new BufferedWriter(
				new FileWriter(countFile));

		// DecimalFormat df = new DecimalFormat("#.####");

		for (String chr : noiseFrequencies.keySet())
		{
			Map<Integer, PositionNoiseFrequencies> chrFreqs = noiseFrequencies
					.get(chr);

			for (Integer position : chrFreqs.keySet())
			{
				Map<Substitution, NoiseFrequencyCollector> positionFreqs = chrFreqs
						.get(position).getModelsMap();

				for (Substitution substitution : positionFreqs.keySet())
				{
					NoiseFrequencyCollector freqCollector = positionFreqs
							.get(substitution);

					// write AF frequencies
					afWriter.write(chr + "\t");
					afWriter.write(position + "\t");
					afWriter.write(freqCollector.getAverageCoverage() + "\t");
					afWriter.write(substitution.getRef() + "\t");
					afWriter.write(substitution.getAlt() + "\t");
					afWriter.write(freqCollector.getAFFrequencyString());
					afWriter.write("\n");

					// write count frquencies
					countWriter.write(chr + "\t");
					countWriter.write(position + "\t");
					countWriter
							.write(freqCollector.getAverageCoverage() + "\t");
					countWriter.write(substitution.getRef() + "\t");
					countWriter.write(substitution.getAlt() + "\t");
					countWriter.write(freqCollector.getCountFrequencyString());
					countWriter.write("\n");
				}
			}
		}

		afWriter.close();
		countWriter.close();

	}

	private static void collectFrequencies(File pileupFile) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(pileupFile));

		String[] words = null;
		Set<String> seenPositions = new HashSet<String>();

		for (String line = reader.readLine(); line != null; line = reader
				.readLine())
		{
			words = line.split("\t");

			String position = words[0] + "\t" + words[1];

			// duplicate position in the same pileup file
			if (seenPositions.contains(position))
			{
				continue;
			}

			seenPositions.add(position);

			// A, C, G, T, D
			int[] counts = new int[5];
			int coverage = 0;
			for (int i = 0; i < 4; i++)
			{
				counts[i] = Integer.parseInt(words[i + 4]);
				coverage += counts[i];
			}

			// add D count
			counts[4] = Integer.parseInt(words[9]);
			coverage += counts[4];

			// apply coverage threshold
			if (coverage < minCoverage)
			{
				continue;
			}

			Map<Integer, PositionNoiseFrequencies> chrFrequencies = noiseFrequencies
					.get(words[0]);

			// add chrFrequencies if doesn't exist already
			if (chrFrequencies == null)
			{
				chrFrequencies = new LinkedHashMap<Integer, PositionNoiseFrequencies>();
				noiseFrequencies.put(words[0], chrFrequencies);
			}

			PositionNoiseFrequencies positionNoiseFrequencies = chrFrequencies
					.get(Integer.parseInt(words[1]));

			// add position noise frequency collector if doesn't exist already
			if (positionNoiseFrequencies == null)
			{
				positionNoiseFrequencies = new PositionNoiseFrequencies();
				chrFrequencies.put(Integer.parseInt(words[1]),
						positionNoiseFrequencies);
			}

			// iterate over read counts and update models
			char ref = words[2].charAt(0);
			int indexToSkip = baseToIndex(ref);
			for (int i = 0; i < counts.length; i++)
			{
				if (i == indexToSkip)
				{
					continue;
				}

				int count = counts[i];
				char alt = indexToBase(i);

				Substitution substitution = Substitution.get(ref, alt);

				NoiseFrequencyCollector frequencyCollector = positionNoiseFrequencies
						.getFrequencyCollectorFor(substitution);
				frequencyCollector.add(count, coverage);
			}
		}

		reader.close();

	}

	private static char indexToBase(int i)
	{
		switch (i)
		{
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
		case 4:
			return 'D';
		default:
			return 'X';
		}

	}

	private static int baseToIndex(char base)
	{
		switch (base)
		{
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		case 'D':
			return 4;
		default:
			return -1;
		}
	}

}
