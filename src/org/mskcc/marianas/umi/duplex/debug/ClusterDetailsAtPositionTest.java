package org.mskcc.marianas.umi.duplex.debug;

import org.mskcc.marianas.umi.duplex.DuplexUMIBamToCollapsedFastqFirstPass;

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
		String positionBam = "/Users/patelj1/workspace/Marianas/bamFiles/PC55-10-42646082.bam";

		String pileupFile = "/Users/patelj1/workspace/Marianas/pileupFiles/PC55-PC55-IGO-05500-CZ-5_bc213_Pool-05500-CZ-Tube1-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";

		String position = "10:42646082-42646082";
		String UMIMismatches = "1";
		String wobble = "2";
		String referenceFasta = "/Users/patelj1/resources/hg19-ncbi/Homo_sapiens_assembly19.fasta";
		String outputFile = "PC55-10-42646082-cluster-details.txt";

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
