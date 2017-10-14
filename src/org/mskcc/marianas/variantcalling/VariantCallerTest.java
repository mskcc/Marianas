/**
 * 
 */
package org.mskcc.marianas.variantcalling;

import java.io.IOException;

/**
 * @author Juber Patel
 *
 */
public class VariantCallerTest
{

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		String tumorPileupFile = "pileupFiles/PC55-PC41-1to1000-IGO-05500-CZ-8_bc216_Pool-05500-CZ-Tube1-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		String normalPileupFile = "pileupFiles/PC41-PC41-IGO-05500-CZ-1_bc209_Pool-05500-CZ-Tube1-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		String sampleName = "PC55-PC41-1-1000";
		String hotspotsFile = "hotspot-list-union-v1-v2.txt";
		String noiseFrequenciesFile = "noise-frequencies.txt";

		VariantCaller.main(new String[] { tumorPileupFile, normalPileupFile,
				sampleName, hotspotsFile, noiseFrequenciesFile });

	}

}
