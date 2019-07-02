/**
 * 
 */
package org.mskcc.marianas.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

/**
 * @author Juber Patel
 * 
 *         A class to hold static resources that can be accessed by any class
 *         that needs them
 *
 */
public class StaticResources
{
	public static final String version = "1.8.1";
	private static IndexedFastaSequenceFile referenceFasta;
	private static FastaSequenceIndex referenceFastaIndex;
	private static Map<String, Map<Integer, Byte[]>> genotypes;
	private static String positionOfInterest;

	public static final Random rng = new Random();

	public static IndexedFastaSequenceFile getReference()
	{
		return referenceFasta;
	}

	public static FastaSequenceIndex getReferenceIndex()
	{
		return referenceFastaIndex;
	}

	public static Map<String, Map<Integer, Byte[]>> getGenotypes()
	{
		return genotypes;
	}

	public static String getPositionOfInterest()
	{
		return positionOfInterest;
	}

	public StaticResources(IndexedFastaSequenceFile referenceFasta,
			FastaSequenceIndex refFastaIndex, File pileupFile,
			String positionOfInterest) throws IOException
	{
		StaticResources.referenceFasta = referenceFasta;
		StaticResources.referenceFastaIndex = refFastaIndex;
		StaticResources.positionOfInterest = positionOfInterest;

		readGenotypes(pileupFile);

		rng.nextInt(10000);
	}

	/**
	 * read this person's genotypes at each position in the pileup file. These
	 * genotypes will be used as fall back genotypes in case no strong evidence
	 * for alternate genotypes is found. A persons genotype at each position is
	 * simply the base with max read count.
	 * 
	 * @param pileupFile
	 * @throws IOException
	 */
	private void readGenotypes(File pileupFile) throws IOException
	{
		genotypes = new HashMap<String, Map<Integer, Byte[]>>();

		// in case the pileup file does not exist or is not readable
		if (!pileupFile.canRead())
		{
			System.err.println(pileupFile.getPath()
					+ " is not readable. Will use reference as genotype for all positions.");
			return;
		}

		BufferedReader reader = new BufferedReader(new FileReader(pileupFile));

		String line = null;
		String[] words = null;
		int[] counts = new int[5];

		Map<Double, Byte> selected = new TreeMap<Double, Byte>();
		while ((line = reader.readLine()) != null)
		{
			words = line.split("\t");
			int position = Integer.parseInt(words[1]);

			counts[0] = Integer.parseInt(words[4]);
			counts[1] = Integer.parseInt(words[5]);
			counts[2] = Integer.parseInt(words[6]);
			counts[3] = Integer.parseInt(words[7]);
			counts[4] = Integer.parseInt(words[9]);

			int total = 0;
			for (int i = 0; i < counts.length; i++)
			{
				total += counts[i];
			}

			// minimum coverage required to read genotype for this position
			if (total < 50)
			{
				continue;
			}

			selected.clear();
			double proportionThreshold = 0.15;

			double proportion = (counts[0] * 1.0) / total;
			if (proportion >= proportionThreshold)
			{
				selected.put(proportion, (byte) 'A');
			}

			proportion = (counts[1] * 1.0) / total;
			if (proportion >= proportionThreshold)
			{
				selected.put(proportion, (byte) 'C');
			}

			proportion = (counts[2] * 1.0) / total;
			if (proportion >= proportionThreshold)
			{
				selected.put(proportion, (byte) 'G');
			}

			proportion = (counts[3] * 1.0) / total;
			if (proportion >= proportionThreshold)
			{
				selected.put(proportion, (byte) 'T');
			}
			
			proportion = (counts[4] * 1.0) / total;
			if (proportion >= proportionThreshold)
			{
				selected.put(proportion, (byte) 'D');
			}

			// store the genotype in the map
			Map<Integer, Byte[]> chrMap = genotypes.get(words[0]);
			if (chrMap == null)
			{
				chrMap = new HashMap<Integer, Byte[]>();
				genotypes.put(words[0], chrMap);
			}
			
			List<Byte> list = new ArrayList<>(selected.values());
			Collections.reverse(list);
			chrMap.put(position, list.toArray(new Byte[0]));
		}

		reader.close();
	}

	/**
	 * read this person's genotypes at each position in the pileup file. These
	 * genotypes will be used as fall back genotypes in case no strong evidence
	 * for alternate genotypes is found. A persons genotype at each position is
	 * simply the base with max read count.
	 * 
	 * @param pileupFile
	 * @throws IOException
	 */
	/*
	 * private void readGenotypesOld(File pileupFile) throws IOException
	 * {
	 * genotypes = new HashMap<String, Map<Integer, Byte>>();
	 * 
	 * BufferedReader reader = new BufferedReader(new FileReader(pileupFile));
	 * 
	 * String line = null;
	 * String[] words = null;
	 * int[] counts = new int[4];
	 * while ((line = reader.readLine()) != null)
	 * {
	 * words = line.split("\t");
	 * int position = Integer.parseInt(words[1]);
	 * counts[0] = Integer.parseInt(words[4]);
	 * counts[1] = Integer.parseInt(words[5]);
	 * counts[2] = Integer.parseInt(words[6]);
	 * counts[3] = Integer.parseInt(words[7]);
	 * 
	 * int max = -1;
	 * int maxIndex = -1;
	 * for (int i = 0; i < counts.length; i++)
	 * {
	 * if (counts[i] > max)
	 * {
	 * max = counts[i];
	 * maxIndex = i;
	 * }
	 * }
	 * 
	 * // get the genotype
	 * byte genotype = 'N';
	 * if (maxIndex == 0)
	 * {
	 * genotype = 'A';
	 * }
	 * else if (maxIndex == 1)
	 * {
	 * genotype = 'C';
	 * }
	 * else if (maxIndex == 2)
	 * {
	 * genotype = 'G';
	 * }
	 * else if (maxIndex == 3)
	 * {
	 * genotype = 'T';
	 * }
	 * 
	 * // store the genotype in the map
	 * Map<Integer, Byte> chrMap = genotypes.get(words[0]);
	 * if (chrMap == null)
	 * {
	 * chrMap = new HashMap<Integer, Byte>();
	 * genotypes.put(words[0], chrMap);
	 * }
	 * 
	 * chrMap.put(position, genotype);
	 * }
	 * 
	 * reader.close();
	 * }
	 */

}
