/**
 * 
 */
package org.mskcc.marianas.util;

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

	public static IndexedFastaSequenceFile getReference()
	{
		return referenceFasta;
	}

	public StaticResources(IndexedFastaSequenceFile referenceFasta)
	{
		this.referenceFasta = referenceFasta;
	}

}
