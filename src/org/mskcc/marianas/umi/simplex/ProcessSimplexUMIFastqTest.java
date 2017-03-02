/**
 * 
 */
package org.mskcc.marianas.umi.simplex;

import java.io.IOException;

import org.mskcc.marianas.umi.duplex.fastqprocessing.ProcessLoebUMIFastq;

/**
 * @author Juber Patel
 *
 */
public class ProcessSimplexUMIFastqTest
{

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		// String read1Fastq = "fastqFiles/BC10_S2_L001_R1_001.fastq.gz";
		String read1Fastq = "fastqFiles/BC12_S4_L001_R1_001.fastq.gz";
		//String UMILength = "12";
		String UMILength = "163";
		// String barcodeIndex = "ACAAGCTA";
		String barcodeIndex = "AGTACAAG";
		String projectFolder = "umi/Project_UMI";

		ProcessLoebUMIFastq.main(new String[] { read1Fastq, UMILength, barcodeIndex,
				projectFolder });
	}
}
