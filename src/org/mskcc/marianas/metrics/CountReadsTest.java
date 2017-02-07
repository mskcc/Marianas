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
		// String bam =
		// "/Users/patelj1/workspace/CustomCapture/MOMO_0089_AC6G07ANXX___P5500_Y___PC11_5_A___hg19___MD.sorted.bam";
		String bam = "/Users/patelj1/workspace/PUMA/5500-AR/FinalBams/DS-puma-0027-PL-C3-IGO-05500-AR-33_bc67_5500-AR_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		String coverageThreshold = "100";
		String geneList = "/Users/patelj1/resources/gene-list/juber-hg19-gene-list.bed";
		String bedFile = "bedFiles/ERBB2.bed";

		CountReads.main(
				new String[] { bam, coverageThreshold, geneList, bedFile });
	}
}
