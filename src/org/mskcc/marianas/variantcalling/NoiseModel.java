/**
 * 
 */
package org.mskcc.marianas.variantcalling;

import org.apache.commons.math3.stat.Frequency;

/**
 * @author Juber Patel
 *
 */
public class NoiseModel
{
	private NoiseModelID id;
	private Frequency frequencyTable;

	public NoiseModel(NoiseModelID id, Frequency frequencyTable)
	{
		this.id = id;
		this.frequencyTable = frequencyTable;
	}

}
