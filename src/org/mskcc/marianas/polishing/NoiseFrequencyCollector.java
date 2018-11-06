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
	private Frequency afFrequencies;
	private Frequency countFrequencies;

	public NoiseFrequencyCollector()
	{
		afFrequencies = new Frequency();
		countFrequencies = new Frequency();
	}

	public void add(int count, int coverage)
	{
		double af = (count * 1.0) / coverage;
		afFrequencies.addValue(af);
		countFrequencies.addValue(count);
		coverageAccumulator += coverage;
	}

	public long getAverageCoverage()
	{
		return (long) (coverageAccumulator / afFrequencies.getSumFreq());
	}

	public String getAFFrequencyString()
	{
		StringBuilder builder = new StringBuilder();

		Iterator iterator = afFrequencies.valuesIterator();
		while (iterator.hasNext())
		{
			Double value = (Double) iterator.next();
			long freq = afFrequencies.getCount(value);
			builder.append("\t").append(value);
			builder.append("\t").append(freq);
		}

		return builder.substring(1);
	}
	
	public String getCountFrequencyString()
	{
		StringBuilder builder = new StringBuilder();

		Iterator iterator = countFrequencies.valuesIterator();
		while (iterator.hasNext())
		{
			Long value = (Long) iterator.next();
			long freq = countFrequencies.getCount(value);
			builder.append("\t").append(value);
			builder.append("\t").append(freq);
		}

		return builder.substring(1);
	}

}
