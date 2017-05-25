/**
 * 
 */
package org.mskcc.marianas.commands.util;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Juber Patel
 *
 *         A Simple class to hold information about a genomic site and genotypes
 *         at that site for different people
 */
public class SiteInfo
{
	String contig;
	int position;
	/**
	 * person -> allele mappings
	 */
	Map<String, Character> allele1;
	Map<String, Character> allele2;

	public SiteInfo(String contig, int position)
	{
		this.contig = contig;
		this.position = position;
		this.allele1 = new HashMap<String, Character>();
		this.allele2 = new HashMap<String, Character>();
	}

	public void addAllele1(String person, char allele)
	{
		allele1.put(person, allele);
	}

	public void addAllele2(String person, char allele)
	{
		allele2.put(person, allele);
	}

}
