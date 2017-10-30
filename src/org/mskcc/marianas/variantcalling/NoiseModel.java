/**
 * 
 */
package org.mskcc.marianas.variantcalling;

import java.util.Iterator;

import org.apache.commons.math3.stat.Frequency;
import org.apache.commons.math3.stat.inference.TestUtils;

/**
 * @author Juber Patel
 *
 */
public class NoiseModel
{
	public final NoiseModelID id;
	public final Frequency frequencyTable;

	public NoiseModel(NoiseModelID id, Frequency frequencyTable)
	{
		this.id = id;
		this.frequencyTable = frequencyTable;
	}

	public long getNumObservations()
	{
		return frequencyTable.getSumFreq();
	}

	public double getPValueFor(double af)
	{
		// convert frequencies to double[]
		double[] values = new double[(int) frequencyTable.getSumFreq()];
		int index = 0;
		Iterator it = frequencyTable.valuesIterator();
		while (it.hasNext())
		{
			Double value = (Double) it.next();
			int freq = (int) frequencyTable.getCount(value);
			for (int i = 0; i < freq; i++)
			{
				values[index] = value;
				index++;
			}
		}

		// compute p-value of the t-test
		return TestUtils.tTest(af, values);
	}

}
