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
	
		//String UMIProcessedBam="/Users/patelj1/workspace/Waltz/bamFiles/PC41-PC41-IGO-05500-CZ-1_bc209_Pool-05500-CZ-Tube1-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		String UMIProcessedBam= "/Users/patelj1/workspace/Waltz/bamFiles/SK-PB-191-G-30ng-1-5uM-IGO-05500-DL-12_bc212_Pool-05500-DL-Tube3-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
				
				

		//String outputFolder = "/Users/patelj1/workspace/Marianas/collapsed-fastqs/Test/DL-12/";
		String outputFolder = "/Users/patelj1/workspace/Marianas/collapsed-fastqs/G-30/";

		
		//String pileupFile = "/Users/patelj1/workspace/Marianas/pileupFiles/PC41-PC41-IGO-05500-CZ-1_bc209_Pool-05500-CZ-Tube1-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		String pileupFile = "/Users/patelj1/workspace/Marianas/pileupFiles/SK-PB-191-G-30ng-1-5uM-IGO-05500-DL-12_bc212_Pool-05500-DL-Tube3-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR-pileup.txt";
		
		
		String UMIMismatches = "1";
		String wobble = "2";
		String referenceFasta = "/Users/patelj1/resources/hg19-ncbi/Homo_sapiens_assembly19.fasta";

		try
		{
			DuplexUMIBamToCollapsedFastqFirstPass.main(
					new String[] { UMIProcessedBam, pileupFile, UMIMismatches,
							wobble, referenceFasta, outputFolder });
		}
		catch (Exception e)
		{
			System.err.println("Exception bubbled to top level:");
			e.printStackTrace();
		}
	}

}
