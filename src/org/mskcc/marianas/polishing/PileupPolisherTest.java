/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.io.IOException;

/**
 * @author Juber Patel
 *
 */
public class PileupPolisherTest
{

	static final String[] pileups = {
			"ML-fgfr2-015-pl-T02_IGO_05500_FC_16_S16_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-015-pl-T07_IGO_05500_FC_17_S17_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-016-pl-T01_IGO_05500_FC_1_S1_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-016-pl-T11_IGO_05500_FC_2_S2_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-017-pl-T02_IGO_05500_FC_3_S3_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-017-pl-T03_IGO_05500_FC_4_S4_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-033-pl-T02_IGO_05500_FC_5_S5_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-033-pl-T07_IGO_05500_FC_6_S6_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-038-pl-T02_IGO_05500_FC_7_S7_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-038-pl-T04_IGO_05500_FC_8_S8_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-040-pl-T04_IGO_05500_FC_18_S18_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-040-pl-T07_IGO_05500_FC_9_S9_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-041-pl-T01_IGO_05500_FC_11_S11_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-041-pl-T07_IGO_05500_FC_10_S10_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-045-pl-T01_IGO_05500_FC_13_S13_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-045-pl-T06_IGO_05500_FC_12_S12_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-053-pl-T00_IGO_05500_FC_14_S14_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt",
			"ML-fgfr2-053-pl-T09_IGO_05500_FC_15_S15_001_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-pileup.txt" };

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		String pileupPath = "/Volumes/innovation/Innovation/ACCESS-Projects/5500-FC/"
				+ "5500-FC-LoweryM-Cholangiocarcinoma-FGFR2-spike-in-0.0.35-hotfix-3-g8fcc361/QC_Results/waltz_duplex_pool_a/";

		String afFrequenciesFile = "/Users/patelj1/workspace/Marianas/polishing/duplex/frequencies-m1/af-frequencies.txt";
		String countFrequenciesFile = "/Users/patelj1/workspace/Marianas/polishing/duplex/frequencies-m1/count-frequencies.txt";

		for (String f : pileups)
		{
			PileupPolisher.main(new String[] { pileupPath + f,
					afFrequenciesFile, countFrequenciesFile });
		}

	}

}
