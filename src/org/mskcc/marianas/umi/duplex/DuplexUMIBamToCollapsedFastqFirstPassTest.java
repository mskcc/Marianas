/**
 * 
 */
package org.mskcc.marianas.umi.duplex;

/**
 * @author Juber Patel
 *
 */
public class DuplexUMIBamToCollapsedFastqFirstPassTest
{

	/**
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args)
	{

		// String UMIProcessedBam =
		// "/Users/patelj1/workspace/Waltz/bamFiles/hapmap-deletion-site.bam";

		// String UMIProcessedBam =
		// "/Users/patelj1/workspace/Waltz/bamFiles/hapmap-insertion-site.bam";

		String UMIProcessedBam = "/Users/patelj1/workspace/Waltz/bamFiles/chr1.bam";

		// String UMIProcessedBam = "/Users/patelj1/workspace/Waltz/bamFiles/"
		/// + "a.bam";

		// String UMIProcessedBam = "/Users/patelj1/workspace/Waltz/bamFiles/"
		// + "DS-fast-005-15-45468636.bam";

		// String UMIProcessedBam =
		// "/Volumes/innovation/Innovation/projects/Juber/HiSeq/5500-DY/"
		// + "run-5500-DY-4/FinalBams/"
		// +
		// "MSK-L-017-cf-IGO-05500-DY-18_bc209_5500-DY-4_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";

		// String UMIProcessedBam = "/Users/patelj1/workspace/Waltz/bamFiles/"
		// +
		// "MSK-L-017-cf-IGO-05500-DY-18_bc209_5500-DY-4_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";

		// String UMIProcessedBam = "/Users/patelj1/workspace/Waltz/bamFiles/"
		// +
		// "MSK-L-051-cf-IGO-05500-DY-21_bc212_5500-DY-4_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";

		//////////////////////////////////////////////////////////////////////////

		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/pileupFiles/pre-collapsing/"
		// +
		// "BL-tdm1-012-pl-T02-10ng-IGO-05500-ES-15_bc403_Pool-05500-ES-Tube3-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";

		String pileupFile = "/Users/patelj1/workspace/Marianas/pileupFiles/pre-collapsing/"
				+ "Hapmap-ctrl-1_S2_001_cl_aln_srt_MD_IR_FX_BR-pileup.txt";

		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/pileupFiles/pre-collapsing/"
		// +
		// "DS-fast086to005-1to100-IGO-05500-ET-29_bc479_Pool-05500-ET-Tube5-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";

		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/pileupFiles/pre-collapsing/"
		// +
		// "MSK-L-017-bc-IGO-05500-DY-1_bc221_5500-DY-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/pileupFiles/pre-collapsing/"
		// +
		// "MSK-L-051-bc-IGO-05500-DY-4_bc220_5500-DY-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";

		// String pileupFile = "whatevs";

		/////////////////////////////////////////////////////////////////////////////////

		// String outputFolder =
		// "/Users/patelj1/workspace/Marianas/collapsed-fastqs/Test/DL-12/";
		// String outputFolder =
		// "/Users/patelj1/workspace/Marianas/collapsed-fastqs/G-30/";
		String outputFolder = "/Users/patelj1/workspace/Marianas/collapsed-fastqs/test/";

		String minMappingQuality = "1";
		String minBaseQuality = "20";
		String UMIMismatches = "0";
		String wobble = "2";
		String minConsensusPercent = "90";

		String referenceFasta = "/Users/patelj1/resources/hg19-ncbi/Homo_sapiens_assembly19.fasta";

		try
		{
			DuplexUMIBamToCollapsedFastqFirstPass.main(new String[] {
					UMIProcessedBam, pileupFile, minMappingQuality,
					minBaseQuality, UMIMismatches, wobble, minConsensusPercent,
					referenceFasta, outputFolder });
		}
		catch (Exception e)
		{
			System.err.println("Exception bubbled to top level:");
			e.printStackTrace();
		}
	}

}
