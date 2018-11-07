/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.io.IOException;

/**
 * @author Juber Patel
 *
 */
public class PolisherTest
{

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		String mafFile = "/Users/patelj1/workspace/ACCESS/Lowery-Cholangiocarcinoma/genotypes.maf";

		String afFrequenciesFile = "/Users/patelj1/workspace/Marianas/polishing/duplex/frequencies-m1/af-frequencies.txt";
		String countFrequenciesFile = "/Users/patelj1/workspace/Marianas/polishing/duplex/frequencies-m1/count-frequencies.txt";

		Polisher.main(new String[] { mafFile, afFrequenciesFile,
				countFrequenciesFile });

	}

}
