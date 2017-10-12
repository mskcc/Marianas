/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
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
public class NoiseModelBuilder
{
	private static final int minCoverage = 50;

	private static Map<String, Map<Integer, PositionNoiseModels>> noiseModels = //
			new HashMap<String, Map<Integer, PositionNoiseModels>>();

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		File directory = new File(args[0]);

		// populate models
		for (File file : directory.listFiles())
		{
			if (!file.getName().endsWith("-pileup.txt"))
			{
				continue;
			}

			updateModels(file);
		}

		// calculate parameters and write them out

		int a = 5;

	}

	private static void updateModels(File file) throws IOException
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

			Map<Integer, PositionNoiseModels> chrModels = noiseModels
					.get(words[0]);

			// add chrModels if doesn't exist already
			if (chrModels == null)
			{
				chrModels = new HashMap<Integer, PositionNoiseModels>();
				noiseModels.put(words[0], chrModels);
			}

			PositionNoiseModels positionNoiseModels = chrModels
					.get(Integer.parseInt(words[1]));

			// add position noise models if doesn't exist already
			if (positionNoiseModels == null)
			{
				positionNoiseModels = new PositionNoiseModels();
				chrModels.put(Integer.parseInt(words[1]), positionNoiseModels);
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
				double af = (count * 1.0) / coverage;
				char alt = pileupIndexToBase(i);

				Substitution substitution = Substitution.get(words[2].charAt(0),
						alt);

				NoiseModel model = positionNoiseModels.getModel(substitution);
				model.add(af);
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
