/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.io.IOException;

/**
 * @author Juber Patel
 *
 */
public class NoiseModelBuilderTest
{

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		String pileupDirectory = "pileupFiles/polishing-test/";

		NoiseModelBuilder.main(new String[] { pileupDirectory });
	}

}
