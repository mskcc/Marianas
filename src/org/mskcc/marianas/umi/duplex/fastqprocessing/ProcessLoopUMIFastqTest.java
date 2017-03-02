/**
 * 
 */
package org.mskcc.marianas.umi.duplex.fastqprocessing;

import java.io.IOException;

/**
 * @author Juber Patel
 *
 */
public class ProcessLoopUMIFastqTest
{

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		String R1fastq = "/Users/patelj1/workspace/Marianas/umi/original-fastqs/Project_05500_CR-test/Sample_SK-PB-191-G-30-Loop_IGO_05500_CR_9/SK-PB-191-G-30-Loop_IGO_05500_CR_9_S75_L004_R1_001.fastq.gz";
		String UMILength = "3";
		String outputProjectFolder = "/Users/patelj1/workspace/Marianas/umi/umi-processed-fastqs/Project-CR";

		ProcessLoopUMIFastq
				.main(new String[] { R1fastq, UMILength, outputProjectFolder });
	}

}
