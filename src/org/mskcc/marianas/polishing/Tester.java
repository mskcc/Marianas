/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.util.Iterator;

import org.apache.commons.math3.stat.Frequency;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.hipparchus.stat.inference.MannWhitneyUTest;

/**
 * @author Juber Patel
 *
 */
public class Tester
{
	private static final MannWhitneyUTest MannWhitneyUTest = new MannWhitneyUTest();

	public static double MannWhitneyUTest(Frequency frequencyTable, double af)
	{
		double[] values = doubleArray(frequencyTable);

		// compute p-value of the Mann Whitney U test
		return MannWhitneyUTest.mannWhitneyUTest(values, new double[] { af });
	}

	public static double tTest(Frequency frequencyTable, double af)
	{
		double[] values = doubleArray(frequencyTable);

		// compute p-value of the t-test
		return TestUtils.tTest(af, values);

	}

	private static double[] doubleArray(Frequency frequencyTable)
	{
		// convert frequencies to double[]
		double[] values = new double[(int) frequencyTable.getSumFreq()];
		int index = 0;
		Iterator it = frequencyTable.valuesIterator();
		while (it.hasNext())
		{
			Object o = it.next();
			double value;
			int freq;

			if (o.getClass().equals(Double.class))
			{
				value = (Double) o;
				freq = (int) frequencyTable.getCount((Double) o);
			}
			else
			{
				value = new Double((Long) o);
				freq = (int) frequencyTable.getCount((Long) o);
			}

			for (int i = 0; i < freq; i++)
			{
				values[index] = value;
				index++;
			}
		}

		return values;
	}

}
