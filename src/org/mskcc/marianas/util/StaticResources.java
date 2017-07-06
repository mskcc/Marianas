/**
 * 
 */
package org.mskcc.marianas.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
	private static IndexedFastaSequenceFile referenceFasta;
	private static FastaSequenceIndex referenceFastaIndex;
	private static Map<String, Map<Integer, Byte[]>> genotypes;
	private static String positionOfInterest;

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
			FastaSequenceIndex refFastaIndex, File pileupFile, String positionOfInterest) throws IOException
	{
		StaticResources.referenceFasta = referenceFasta;
		StaticResources.referenceFastaIndex = refFastaIndex;
		StaticResources.positionOfInterest = positionOfInterest;

		readGenotypes(pileupFile);
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

		BufferedReader reader = new BufferedReader(new FileReader(pileupFile));

		String line = null;
		String[] words = null;
		int[] counts = new int[4];
		List<Byte> selected = new ArrayList<Byte>();
		while ((line = reader.readLine()) != null)
		{
			words = line.split("\t");
			int position = Integer.parseInt(words[1]);
			int total = Integer.parseInt(words[3]);
			counts[0] = Integer.parseInt(words[4]);
			counts[1] = Integer.parseInt(words[5]);
			counts[2] = Integer.parseInt(words[6]);
			counts[3] = Integer.parseInt(words[7]);

			selected.clear();
			double proportionThreshold = 0.15;

			if ((counts[0] * 1.0) / total >= proportionThreshold)
			{
				selected.add((byte) 'A');
			}

			if ((counts[1] * 1.0) / total >= proportionThreshold)
			{
				selected.add((byte) 'C');
			}

			if ((counts[2] * 1.0) / total >= proportionThreshold)
			{
				selected.add((byte) 'G');
			}

			if ((counts[3] * 1.0) / total >= proportionThreshold)
			{
				selected.add((byte) 'T');
			}
			
			// store the genotype in the map
			Map<Integer, Byte[]> chrMap = genotypes.get(words[0]);
			if (chrMap == null)
			{
				chrMap = new HashMap<Integer, Byte[]>();
				genotypes.put(words[0], chrMap);
			}

			chrMap.put(position, selected.toArray(new Byte[0]));
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
