/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.util.Iterator;

import org.apache.commons.math3.stat.Frequency;

/**
 * @author Juber Patel
 * 
 *         Noise model for a position for a nucleotide change
 *
 */
public class NoiseFrequencyCollector
{
	private double coverageAccumulator;
	private Frequency values;

	public NoiseFrequencyCollector()
	{
		values = new Frequency();
	}

	public void add(int count, int coverage)
	{
		double af = (count * 1.0) / coverage;
		values.addValue(af);
		coverageAccumulator += coverage;
	}

	public long getAverageCoverage()
	{
		return (long) (coverageAccumulator / values.getSumFreq());
	}

	public String getFrequencyString()
	{
		StringBuilder builder = new StringBuilder();

		Iterator iterator = values.valuesIterator();
		while (iterator.hasNext())
		{
			Double value = (Double) iterator.next();
			long freq = values.getCount(value);
			builder.append("\t").append(value);
			builder.append("\t").append(freq);
		}

		return builder.substring(1);
	}

}
