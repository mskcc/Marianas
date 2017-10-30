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
		// String tumorPileupFile =
		// "pileupFiles/PC55-PC41-1to1000-IGO-05500-CZ-8_bc216_Pool-05500-CZ-Tube1-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		// String normalPileupFile =
		// "pileupFiles/PC41-PC41-IGO-05500-CZ-1_bc209_Pool-05500-CZ-Tube1-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";

		// String tumorPileupFile =
		// "pileupFiles/post-collapsing/lung/MSK-L-007-cf-IGO-05500-DY-22_bc213_5500-DY-4_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		// String normalPileupFile =
		// "pileupFiles/post-collapsing/lung/MSK-L-007-bc-IGO-05500-DY-5_bc217_5500-DY-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";

		String tumorPileupFile = "pileupFiles/post-collapsing/lung/"
				+ "MSK-L-017-cf-IGO-05500-DY-18_bc209_5500-DY-4_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		String normalPileupFile = "pileupFiles/post-collapsing/lung/"
				+ "MSK-L-017-bc-IGO-05500-DY-1_bc221_5500-DY-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";

		String sampleName = "MSK-L-017";
		String hotspotsFile = "hotspot-list-union-v1-v2.txt";
		String noiseFrequenciesFile = "noise-frequencies.txt";

		VariantCaller.main(new String[] { tumorPileupFile, normalPileupFile,
				sampleName, hotspotsFile, noiseFrequenciesFile });

	}

}
