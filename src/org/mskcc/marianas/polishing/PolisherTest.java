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

		String afFrequenciesFile = "/Users/patelj1/workspace/Marianas/polishing/novaseq-duplex/af-frequencies.txt";
		String countFrequenciesFile = "/Users/patelj1/workspace/Marianas/polishing/novaseq-duplex/count-frequencies.txt";

		String depthColumnName = "SD_t_depth_count_fragment";
		String altColumnName = "SD_t_alt_count_fragment";

		Polisher.main(new String[] { mafFile, depthColumnName, altColumnName,
				afFrequenciesFile, countFrequenciesFile });

	}

}
