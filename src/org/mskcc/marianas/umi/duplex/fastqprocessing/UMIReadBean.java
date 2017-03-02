/**
 * 
 */
package org.mskcc.marianas.umi.duplex.fastqprocessing;

/**
 * @author Juber Patel
 *
 */
public class UMIReadBean
{
	String seq;
	String qual;
	String UMI;
	String UMIQual;

	public UMIReadBean(String seq, String qual, String UMI, String UMIQual)
	{
		this.seq = seq;
		this.qual = qual;
		this.UMI = UMI;
		this.UMIQual = UMIQual;

	}

}
