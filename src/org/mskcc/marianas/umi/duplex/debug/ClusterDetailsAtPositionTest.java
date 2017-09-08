package org.mskcc.marianas.umi.duplex.debug;

/**
 * @author Juber Patel
 *
 */
public class ClusterDetailsAtPositionTest
{
	/**
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args)
	{
		// String UMIProcessedBam =
		// "/Users/patelj1/workspace/Marianas/umi/umi-processed-bams/Test/FinalBams/Test-1.bam";
		String positionBam = "/Users/patelj1/workspace/Marianas/MSK-L-046-cf-IGO-05500-DY-25_bc221_5500-DY-5_L000_mrg_cl_aln_srt_MD_IR_FX_BR-7:55259515-55259515.bam";

		String pileupFile = "/Users/patelj1/workspace/Marianas/pileupFiles/MSK-L-046-bc-IGO-05500-DY-8_bc211_5500-DY-2_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";

		String position = "7:55259515-55259515";
		String UMIMismatches = "1";
		String wobble = "2";
		String referenceFasta = "/Users/patelj1/resources/hg19-ncbi/Homo_sapiens_assembly19.fasta";
		String outputFile = "MSK-L-046-cf-7:55259515-cluster-details.txt";

		try
		{
			ClusterDetailsAtPosition
					.main(new String[] { positionBam, pileupFile, position,
							UMIMismatches, wobble, referenceFasta, outputFile });
		}
		catch (Exception e)
		{
			System.err.println("Exception bubbled to top level:");
			e.printStackTrace();
		}
	}

}
