/**
 * 
 */
package org.mskcc.marianas.umi.loeb.withcorrection;

/**
 * @author Juber Patel
 *
 */
public class ProcessDuplexUMIBamSecondPassTest
{

	/**
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception
	{
		// String bamFile =
		// "/Users/patelj1/workspace/Marianas/umi/run-5500-CD/FinalBams/100ng-test-1.bam";
		String bamFile = "/Users/patelj1/workspace/Marianas/umi/run-5500-CD/FinalBams/139417-100ng-IGO-05500-CD-4_bc207_5500-CD_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";

		// String bamFile = "bamFiles/t.bam";
		// String bamFile =
		// "bamFiles/BC09_bc09_Project-5500-AV_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";

		String bedFile = "/Users/patelj1/workspace/Marianas/bedFiles/Sarath-10-genes.bed";
		String UMIMismatches = "2";
		String wobble = "2";
		String referenceFasta = "/Users/patelj1/resources/hg19-ncbi/Homo_sapiens_assembly19.fasta";
		String outputFolder = "/Users/patelj1/workspace/Marianas/umi/umi-processed-fastqs/Project-CD/Sample_139417-100ng_IGO_05500_CD_4";

		ProcessDuplexUMIBamSecondPass.main(new String[] { bamFile, bedFile,
				UMIMismatches, wobble, referenceFasta, outputFolder });

	}

}
