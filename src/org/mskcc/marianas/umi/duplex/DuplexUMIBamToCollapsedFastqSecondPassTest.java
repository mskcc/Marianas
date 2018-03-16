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
		// "/Users/patelj1/workspace/Marianas/umi/umi-processed-bams/Test/FinalBams/Test-1.bam";
		// String UMIProcessedBam =
		// "/Users/patelj1/workspace/Marianas/umi/umi-processed-bams/run-5500-CR-G/FinalBams/SK-PB-191-G-30-Loop-IGO-05500-CR-9_bc212_05500-CR_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		// String
		// UMIProcessedBam="/Users/patelj1/workspace/Marianas/bamFiles/SK-PB-191-G-30ng-1-5uM-IGO-05500-DL-12_bc212_Pool-05500-DL-Tube3-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		// String
		// UMIProcessedBam="/Users/patelj1/workspace/Waltz/bamFiles/MSK-L-035-bc-IGO-05500-DY-7_bc210_5500-DY-2_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		// String UMIProcessedBam =
		// "/Users/patelj1/workspace/Waltz/bamFiles/MSK-L-017-cf-IGO-05500-DY-18_bc209_5500-DY-4_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		String UMIProcessedBam = "/Users/patelj1/workspace/Waltz/bamFiles/"
				+ "MSK-L-051-cf-IGO-05500-DY-21_bc212_5500-DY-4_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";

		// String outputFolder =
		// "/Users/patelj1/workspace/Marianas/umi/collapsed-fastqs/Test/Test-1/";
		// String outputFolder =
		// "/Users/patelj1/workspace/Marianas/umi/collapsed-fastqs/run-5500-CR-G/G-30/";
		// String outputFolder =
		// "/Users/patelj1/workspace/Marianas/umi/collapsed-fastqs/Test/DL-12/";
		String outputFolder = "/Users/patelj1/workspace/Marianas/collapsed-fastqs/test/";

		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/bam-metrics/custom-panel-5500-CR-G-normal-bams/SK-PB-191-G-30-Loop-IGO-05500-CR-9_bc212_05500-CR_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/pileupFiles/PC41-PC41-IGO-05500-CZ-1_bc209_Pool-05500-CZ-Tube1-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/pileupFiles/SK-PB-191-G-30ng-1-5uM-IGO-05500-DL-12_bc212_Pool-05500-DL-Tube3-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		// String pileupFile =
		// "/Users/patelj1/workspace/Marianas/pileupFiles/MSK-L-035-bc-IGO-05500-DY-7_bc210_5500-DY-2_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		//String pileupFile = "/Users/patelj1/workspace/Marianas/pileupFiles/pre-collapsing/"
		//		+ "MSK-L-017-bc-IGO-05500-DY-1_bc221_5500-DY-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		String pileupFile = "whatevs";
		
		
		String UMIMismatches = "1";
		String wobble = "2";
		String minConsensusPercent = "90";
		String referenceFasta = "/Users/patelj1/resources/hg19-ncbi/Homo_sapiens_assembly19.fasta";

		try
		{
			DuplexUMIBamToCollapsedFastqSecondPass.main(new String[] {
					UMIProcessedBam, pileupFile, UMIMismatches, wobble,
					minConsensusPercent, referenceFasta, outputFolder });

		}
		catch (Exception e)
		{
			System.err.println("Exception bubbled to top level:");
			e.printStackTrace();
		}

	}

}
