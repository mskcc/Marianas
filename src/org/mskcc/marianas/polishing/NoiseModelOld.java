/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.stat.Frequency;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

/**
 * @author Juber Patel
 * 
 *         Noise model for a position for a nucleotide change
 *
 */
public class NoiseModelOld
{
	private static GaussianCurveFitter curveFitter = GaussianCurveFitter
			.create();

	private double coverageAccumulator;
	private DescriptiveStatistics values;

	public NoiseModelOld()
	{
		values = new DescriptiveStatistics();
	}

	public void add(int count, int coverage)
	{
		double af = (count * 1.0) / coverage;
		values.addValue(af);
		coverageAccumulator += coverage;
	}

	public long getAverageCoverage()
	{
		return (long) (coverageAccumulator / values.getN());
	}

	public String getModelString()
	{
		StringBuilder builder1 = new StringBuilder();

		double initialSigma = values.getStandardDeviation();
		double initialMean = values.getMean();
		double initialNormalizationFactor = 1
				/ (initialSigma * FastMath.sqrt(2 * Math.PI));

		// make a frequency table
		Frequency frequencyTable = new Frequency();
		double[] v = values.getValues();
		for (int i = 0; i < v.length; i++)
		{
			frequencyTable.addValue(v[i]);
		}

		// collect (x, y) points from the frequency table
		WeightedObservedPoints points = new WeightedObservedPoints();
		Iterator iterator = frequencyTable.valuesIterator();
		while (iterator.hasNext())
		{
			Double value = (Double) iterator.next();
			long freq = frequencyTable.getCount(value);
			points.add(value, freq);
			builder1.append("\t").append(value);
			builder1.append("\t").append(freq);
		}

		StringBuilder builder2 = new StringBuilder();

		if (initialSigma == 0)
		{
			builder2.append(initialNormalizationFactor).append("\t")
					.append(initialMean).append("\t").append(initialSigma);
		}
		else
		{
			// Fit the curve!!!
			curveFitter = curveFitter.withStartPoint(new double[] {
					initialNormalizationFactor, initialMean, initialSigma });
			double[] parameters = curveFitter.fit(points.toList());

			builder2.append(Double.toString(parameters[0]));
			for (int i = 1; i < parameters.length; i++)
			{
				builder2.append("\t").append(parameters[i]);
			}

		}

		return builder2.toString() + builder1.toString();
	}

}
