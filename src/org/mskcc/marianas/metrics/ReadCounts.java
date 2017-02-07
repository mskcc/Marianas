/**
 * 
 */
package org.mskcc.marianas.metrics;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

/**
 * @author Juber Patel
 * 
 *         Simple bean class to hold all the read counts for a bam file
 *
 */
public class ReadCounts
{
	protected final String bamFileName;
	protected final String sampleID;
	protected final String targetLabel;
	protected long totalReads;
	protected long unmappedReads;
	protected long totalMappedReads;
	protected long duplicateMappedReads;
	protected long uniqueMappedReads;
	protected long totalTargetReads;
	protected long uniqueTargetReads;

	public ReadCounts(String bamFileName, String sampleID, String targetLabel)
	{
		this.bamFileName = bamFileName;
		this.sampleID = sampleID;
		this.targetLabel = targetLabel;
	}

	public void write() throws IOException
	{
		// write bam-wide numbers
		BufferedWriter writer = new BufferedWriter(
				new FileWriter(bamFileName + ".read-counts"));
		writer.write(bamFileName + "\t");
		writer.write(totalReads + "\t");
		writer.write(unmappedReads + "\t");
		writer.write(totalMappedReads + "\t");
		writer.write(uniqueMappedReads + "\t");
		writer.write(((duplicateMappedReads * 1.0) / totalMappedReads) + "\t");

		double onTargetTotalFraction = (totalTargetReads * 1.0)
				/ totalMappedReads;
		double onTargetUniqueFraction = (uniqueTargetReads * 1.0)
				/ uniqueMappedReads;
		DecimalFormat decimalFormat = new DecimalFormat("#.##");

		// write target-specific numbers
		writer.write(totalTargetReads + "\t");
		writer.write(uniqueTargetReads + "\t");
		writer.write(decimalFormat.format(onTargetTotalFraction) + "\t");
		writer.write(decimalFormat.format(onTargetUniqueFraction) + "\n");
		
		writer.close();
	}

}
