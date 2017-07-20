/**
 * 
 */
package org.mskcc.marianas.metrics.old;

import java.io.IOException;

import org.mskcc.juber.util.CustomCaptureException;

/**
 * @author Juber Patel
 *
 */
public class ScanBamTest
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
		String bam = "/Users/patelj1/workspace/Shukla/FinalBams/ES-CTDNA-15-01-IGO-05500-AQ-4_bc42_5500-AQ_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		String coverageThreshold = "100";
		String geneList = "/Users/patelj1/resources/gene-list/genelist.with_aa.interval_list";

		ScanBam.main(new String[] { bam, coverageThreshold, geneList });
	}
}
