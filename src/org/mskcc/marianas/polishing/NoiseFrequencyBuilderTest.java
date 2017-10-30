/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.io.IOException;

/**
 * @author Juber Patel
 *
 */
public class NoiseFrequencyBuilderTest
{

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		String pileupDirectory = "pileupFiles/post-collapsing/polishing-normals-fulcrum-1-1-mq20";

		NoiseFrequencyBuilder.main(new String[] { pileupDirectory });
	}

}
