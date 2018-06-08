/**
 * 
 */
package org.mskcc.marianas.commands.util;

/**
 * @author Juber Patel
 * 
 *         A simple class that holds the spike-in info for a site.
 *
 */
public class SpikeInInfo
{
	public String personOfInterest;
	public String otherPerson;
	public int personOfInterestProportion;
	public int otherPersonProportion;
	public String method;
	public SiteInfo siteInfo;
	public int a;
	public int c;
	public int g;
	public int t;
	public int total;
	public char alleleOfInterest;
	public int alleleOfInterestCount;
	public double expectedAF;
	public double observedAF;

	public SpikeInInfo(SiteInfo siteInfo, String person1, String person2,
			int person1Proportion, int person2Proportion, String method)
	{
		this.siteInfo = siteInfo;

		// the person with lower proportion is the person of interest
		if (person1Proportion < person2Proportion)
		{
			personOfInterest = person1;
			otherPerson = person2;
			personOfInterestProportion = person1Proportion;
			otherPersonProportion = person2Proportion;
		}
		else
		{
			personOfInterest = person2;
			otherPerson = person1;
			personOfInterestProportion = person2Proportion;
			otherPersonProportion = person1Proportion;

		}

		this.method = method;

	}

	public void addBaseCounts(String[] pileupWords)
	{
		this.a = Integer.parseInt(pileupWords[4]);
		this.c = Integer.parseInt(pileupWords[5]);
		this.g = Integer.parseInt(pileupWords[6]);
		this.t = Integer.parseInt(pileupWords[7]);
		this.total = a + c + g + t;

		process();
	}

	/**
	 * 1. figure out allele of interest
	 * 2. calculate expected AF
	 * 3. calculate observed AF
	 */
	private void process()
	{
		setAlleleOfInterest();
		setExpectedAF();
		setObservedAF();
	}

	private void setObservedAF()
	{
		alleleOfInterestCount = -1;
		if (alleleOfInterest == 'A')
		{
			alleleOfInterestCount = a;
		}
		else if (alleleOfInterest == 'C')
		{
			alleleOfInterestCount = c;
		}
		else if (alleleOfInterest == 'G')
		{
			alleleOfInterestCount = g;
		}
		else if (alleleOfInterest == 'T')
		{
			alleleOfInterestCount = t;
		}

		observedAF = ((double) alleleOfInterestCount) / total;

	}

	private void setExpectedAF()
	{
		int dose1 = 1;
		int dose2 = 0;

		char allele1 = siteInfo.allele1.get(personOfInterest);
		char allele2 = siteInfo.allele2.get(personOfInterest);

		if (allele1 == alleleOfInterest && allele2 == alleleOfInterest)
		{
			dose1 = 2;
		}

		allele1 = siteInfo.allele1.get(otherPerson);
		allele2 = siteInfo.allele2.get(otherPerson);

		if (allele1 == alleleOfInterest)
		{
			dose2++;
		}

		if (allele2 == alleleOfInterest)
		{
			dose2++;
		}

		int alleleCount = (dose1 * personOfInterestProportion)
				+ (dose2 * otherPersonProportion);

		expectedAF = ((double) alleleCount)
				/ (2 * (personOfInterestProportion + otherPersonProportion));
	}

	private void setAlleleOfInterest()
	{
		// if this is not a spike in, just set the first allele from other
		// person as allele of interest
		if (personOfInterestProportion == 0)
		{
			alleleOfInterest = siteInfo.allele1.get(otherPerson);
			return;
		}

		// the allele from the person with lower proportion that is not present
		// in the other person is the allele of interest
		// assuming the site is biallelic
		char allele1I = siteInfo.allele1.get(personOfInterest);
		char allele2I = siteInfo.allele2.get(personOfInterest);
		char allele1O = siteInfo.allele1.get(otherPerson);
		char allele2O = siteInfo.allele2.get(otherPerson);

		// if the two alleles in the person of interest are same, you don't have
		// a choice
		if (allele1I == allele2I)
		{
			alleleOfInterest = allele1I;
			return;
		}

		if (allele1I != allele1O && allele1I != allele2O)
		{
			alleleOfInterest = allele1I;
			return;
		}

		if (allele2I != allele1O && allele2I != allele2O)
		{
			alleleOfInterest = allele2I;
			return;
		}
	}

	public String getIdentifier()
	{
		return personOfInterest + "-" + otherPerson + "-"
				+ personOfInterestProportion + "-" + otherPersonProportion;
	}

	public String getSite()
	{
		return siteInfo.contig + ":" + siteInfo.position;
	}

	public String getPersonOfInterestGenotype()
	{
		return new StringBuilder()
				.append(siteInfo.allele1.get(personOfInterest))
				.append(siteInfo.allele2.get(personOfInterest)).toString();
	}

	public String getOtherPersonGenotype()
	{
		return new StringBuilder().append(siteInfo.allele1.get(otherPerson))
				.append(siteInfo.allele2.get(otherPerson)).toString();
	}

}
