/**
 * 
 */
package org.mskcc.marianas.commands;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.mskcc.marianas.commands.util.SiteInfo;
import org.mskcc.marianas.commands.util.SpikeInInfo;

/**
 * @author Juber Patel
 *
 */
public class FindSpikedInAlleleFractions
{

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		File sitesFile = new File("005-086-different-genotypes.txt");
		File spikeInInfoFile = new File("spikein-info.txt");

		Map<String, SiteInfo> sites = loadSites(sitesFile);
		SpikeInInfo[] spikeInInfos = computeSpikeInInfo(sites, spikeInInfoFile);

		DecimalFormat df = new DecimalFormat("#.#######");

		// write the output
		BufferedWriter writer = new BufferedWriter(
				new FileWriter("005-086-spike-in-results.txt"));
		writer.write(
				"Name\tMethod\tSite\tGenotype1\tGenotype2\tAlleleOfInterest\tExpectedAF\tObservedAF\tAlleleCount\tTotalCount\n");

		for (SpikeInInfo spikeInInfo : spikeInInfos)
		{
			writer.write(spikeInInfo.getIdentifier() + "\t");
			writer.write(spikeInInfo.method + "\t");
			writer.write(spikeInInfo.getSite() + "\t");
			writer.write(spikeInInfo.getPersonOfInterestGenotype() + "\t");
			writer.write(spikeInInfo.getOtherPersonGenotype() + "\t");
			writer.write(spikeInInfo.alleleOfInterest + "\t");
			writer.write(df.format(spikeInInfo.expectedAF) + "\t");
			writer.write(df.format(spikeInInfo.observedAF) + "\t");
			writer.write(spikeInInfo.alleleOfInterestCount + "\t");
			writer.write(spikeInInfo.total + "\n");
		}

		writer.close();

	}

	private static SpikeInInfo[] computeSpikeInInfo(Map<String, SiteInfo> sites,
			File spikeInInfoFile) throws IOException
	{
		BufferedReader reader = new BufferedReader(
				new FileReader(spikeInInfoFile));

		List<SpikeInInfo> spikeInInfos = new ArrayList<SpikeInInfo>();
		String line;
		String[] words;
		int NR = 0;
		String[] persons = null;
		while ((line = reader.readLine()) != null)
		{
			words = line.split("\t");
			NR++;

			if (NR == 1)
			{
				persons = new String[2];
				System.arraycopy(words, 0, persons, 0, persons.length);
				continue;
			}

			// read the pileup file line by line and collect base counts for
			// each site
			// make a new SpikeInInfo for each site for each file
			BufferedReader pileupReader = new BufferedReader(
					new FileReader(words[3]));
			String pileupLine = null;
			String[] pileupWords = null;
			while ((pileupLine = pileupReader.readLine()) != null)
			{
				pileupWords = pileupLine.split("\t");
				SiteInfo site = sites
						.get(pileupWords[0] + "\t" + pileupWords[1]);
				if (site == null)
				{
					continue;
				}

				SpikeInInfo spikeInInfo = new SpikeInInfo(site, persons[0],
						persons[1], Integer.parseInt(words[0]),
						Integer.parseInt(words[1]), words[2]);

				spikeInInfo.addBaseCounts(pileupWords);
				spikeInInfos.add(spikeInInfo);
			}

			pileupReader.close();

		}

		reader.close();

		return spikeInInfos.toArray(new SpikeInInfo[0]);
	}

	private static Map<String, SiteInfo> loadSites(File sitesFile)
			throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(sitesFile));

		Map<String, SiteInfo> sites = new HashMap<String, SiteInfo>();
		String line;
		String[] words;
		int NR = 0;
		String[] persons = null;
		while ((line = reader.readLine()) != null)
		{
			words = line.split("\t");
			NR++;

			if (NR == 1)
			{
				persons = new String[words.length - 2];
				System.arraycopy(words, 2, persons, 0, persons.length);
				continue;
			}

			SiteInfo site = new SiteInfo(words[0], Integer.parseInt(words[1]));
			for (int i = 0; i < persons.length; i++)
			{
				site.addAllele1(persons[i], words[i + 2].charAt(0));
				site.addAllele2(persons[i], words[i + 2].charAt(1));
			}

			sites.put(words[0] + "\t" + words[1], site);
		}

		reader.close();

		return sites;

	}

}
