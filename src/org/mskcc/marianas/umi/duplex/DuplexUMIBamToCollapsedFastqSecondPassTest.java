/**
 * 
 */
package org.mskcc.marianas.umi.duplex;

/**
 * @author Juber Patel
 *
 */
public class DuplexUMIBamToCollapsedFastqSecondPassTest
{

	/**
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args)
	{
		// String UMIProcessedBam =
		// "/Users/patelj1/workspace/Waltz/bamFiles/MSK-L-017-cf-IGO-05500-DY-18_bc209_5500-DY-4_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		// String UMIProcessedBam = "/Users/patelj1/workspace/Waltz/bamFiles/"
		// +
		// "MSK-L-051-cf-IGO-05500-DY-21_bc212_5500-DY-4_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		String UMIProcessedBam = "/Users/patelj1/workspace/Waltz/bamFiles/chr1.bam";

		///////////////////////////////////////

		String pileupFile = "/Users/patelj1/workspace/Marianas/pileupFiles/pre-collapsing/"
				+ "BL-tdm1-012-pl-T02-10ng-IGO-05500-ES-15_bc403_Pool-05500-ES-Tube3-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";

		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/bam-metrics/custom-panel-5500-CR-G-normal-bams/SK-PB-191-G-30-Loop-IGO-05500-CR-9_bc212_05500-CR_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/pileupFiles/PC41-PC41-IGO-05500-CZ-1_bc209_Pool-05500-CZ-Tube1-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/pileupFiles/SK-PB-191-G-30ng-1-5uM-IGO-05500-DL-12_bc212_Pool-05500-DL-Tube3-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/pileupFiles/MSK-L-035-bc-IGO-05500-DY-7_bc210_5500-DY-2_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/pileupFiles/pre-collapsing/"
		// +
		// "MSK-L-017-bc-IGO-05500-DY-1_bc221_5500-DY-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		// String pileupFile = "whatevs";

		// String outputFolder =
		// "/Users/patelj1/workspace/Marianas/collapsed-fastqs/test/";
		String minMappingQuality = "1";
		String minBaseQuality = "20";
		String UMIMismatches = "0";
		String wobble = "2";
		String minConsensusPercent = "90";
		String referenceFasta = "/Users/patelj1/resources/hg19-ncbi/Homo_sapiens_assembly19.fasta";
		String firstPassFile = "first-pass.mate-position-sorted.txt";

		try
		{
			DuplexUMIBamToCollapsedFastqSecondPass.main(new String[] {
					UMIProcessedBam, pileupFile, minMappingQuality,
					minBaseQuality, UMIMismatches, wobble, minConsensusPercent,
					referenceFasta, firstPassFile });

		}
		catch (Exception e)
		{
			System.err.println("Exception bubbled to top level:");
			e.printStackTrace();
		}

	}

}
