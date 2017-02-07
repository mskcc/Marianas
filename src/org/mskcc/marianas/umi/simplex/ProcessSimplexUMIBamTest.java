/**
 * 
 */
package org.mskcc.marianas.umi.simplex;

import java.io.IOException;

import org.mskcc.marianas.util.CustomCaptureException;

/**
 * @author Juber Patel
 *
 */
public class ProcessSimplexUMIBamTest
{

	/**
	 * @param args
	 * @throws CustomCaptureException
	 * @throws IOException
	 */
	public static void main(String[] args)
			throws IOException, CustomCaptureException
	{
		String bamFile = "bamFiles/BC12_bc12_Project-5500-AV_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		// String bamFile = "bamFiles/t.bam";
		// String bamFile =
		// "bamFiles/BC09_bc09_Project-5500-AV_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		String bedFile = "bedFiles/AL-Switch.bed";
		String UMILength = "12";
		String UMIMismatches = "2";
		String wobble = "2";

		ProcessSimplexUMIBam.main(new String[] { bamFile, bedFile, UMILength,
				UMIMismatches, wobble });

	}

}
