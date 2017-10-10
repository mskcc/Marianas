/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Juber Patel
 * 
 *         Noise model for a position for a nucleotide change
 *
 */
public class NoiseModel
{
	private List<Double> values;

	public NoiseModel()
	{
		values = new ArrayList<Double>();
	}

	public void add(double value)
	{
		values.add(value);

	}

}
