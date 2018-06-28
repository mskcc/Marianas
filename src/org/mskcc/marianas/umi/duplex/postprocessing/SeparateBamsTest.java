/**
 * 
 */
package org.mskcc.marianas.umi.duplex.postprocessing;

import java.io.IOException;

/**
 * @author Juber Patel
 *
 */
public class SeparateBamsTest
{

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		String collapsedBam = "../Waltz/bamFiles/collapsed/"
				+ "MSK-L-017-cf-IGO-05500-DY-18_bc209_5500-DY-4_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";

		SeparateBamsSimplexOnly.main(new String[] { collapsedBam });
	}

}
