/**
 * 
 */
package org.mskcc.marianas.metrics;

import java.io.IOException;

import org.mskcc.marianas.util.CustomCaptureException;

/**
 * @author Juber Patel
 *
 */
public class CountReadsTest
{

	/**
	 * @param args
	 * @throws IOException
	 * @throws CustomCaptureException
	 */
	public static void main(String[] args)
			throws IOException, CustomCaptureException
	{
		String geneList = "/Users/patelj1/resources/gene-list/juber-hg19-gene-list.bed";
		String coverageThreshold = "100";
		
		
		//String bam = "/Users/patelj1/workspace/PUMA/5500-AR/FinalBams/DS-puma-0027-PL-C3-IGO-05500-AR-33_bc67_5500-AR_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		//String bedFile = "bedFiles/ERBB2.bed";
		
		String bam = "bamFiles/MCC_P-0014336-T01_IGO_05500_DG_11_S79_L004.bam";
		String bedFile = "bedFiles/impact410-mcpyv-ebv-hpv.bed";

		CountReads.main(
				new String[] { bam, coverageThreshold, geneList, bedFile });
	}
}
