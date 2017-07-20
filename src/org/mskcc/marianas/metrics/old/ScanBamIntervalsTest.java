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
public class ScanBamIntervalsTest
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
		// "bamFiles/MOMO_0089_AC6G07ANXX___P5500_Y___PC11_5_A___hg19___MD.sorted.bam";
		// String bam =
		// "/Users/patelj1/workspace/Moynahan/FinalBams/1196-2-IGO-05500-AL-21_bc37_5500-AL_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		String bam = "/Users/patelj1/workspace/PUMA/5500-AR/FinalBams/DS-puma-0027-PL-C3-IGO-05500-AR-33_bc67_5500-AR_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";

		// String impactBed =
		// "/Users/patelj1/workspace/CustomCapture/bedFiles/oldBedFiles/impact.bed";
		// String customBed =
		// "/Users/patelj1/workspace/Marianas/bedFiles/BRAF.bed";
		String customBed = "bedFiles/ERBB2.bed";

		ScanBamIntervals.main(new String[] { bam, customBed });
	}
}
