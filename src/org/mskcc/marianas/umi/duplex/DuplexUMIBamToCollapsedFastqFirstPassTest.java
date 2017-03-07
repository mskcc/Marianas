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
		String UMIProcessedBam = "/Users/patelj1/workspace/Marianas/umi/umi-processed-bams/Test/FinalBams/Test-1.bam";
		String outputFolder = "/Users/patelj1/workspace/Marianas/umi/collapsed-fastqs/Test/Test-1/";

		String bedFile = "/Users/patelj1/workspace/Marianas/bedFiles/Sarath-10-genes.bed";
		String UMIMismatches = "1";
		String wobble = "2";
		String referenceFasta = "/Users/patelj1/resources/hg19-ncbi/Homo_sapiens_assembly19.fasta";

		try
		{
			DuplexUMIBamToCollapsedFastqFirstPass.main(
					new String[] { UMIProcessedBam, bedFile, UMIMismatches,
							wobble, referenceFasta, outputFolder });
		}
		catch (Exception e)
		{
			System.err.println("Exception bubbled to top level:");
			e.printStackTrace();
		}
	}

}
