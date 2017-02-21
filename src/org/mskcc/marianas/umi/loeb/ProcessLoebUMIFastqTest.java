/**
 * 
 */
package org.mskcc.marianas.umi.loeb;

import java.io.IOException;

/**
 * @author Juber Patel
 *
 */
public class ProcessLoebUMIFastqTest
{

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		String R1fastq = "/Users/patelj1/workspace/Marianas/umi/original-fastqs/Project_05500_CD/Sample_139417-100ng_IGO_05500_CD_4/139417-100ng_IGO_05500_CD_4_S80_L005_R1_001.fastq.gz";
		String UMILength = "10";
		String constantRegion = "TGACT";
		String projectFolder = "/Users/patelj1/workspace/Marianas/umi/Project-CD-2";

		ProcessLoebUMIFastq.main(new String[] { R1fastq, UMILength,
				constantRegion, projectFolder });
	}

}
