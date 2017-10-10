/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.util.Map;
import java.util.TreeMap;

/**
 * @author Juber Patel
 * 
 *         collection of noise models at a genomic position, one for each
 *         possible substitution
 *
 */
public class PositionNoiseModels
{
	private Map<Substitution, NoiseModel> noiseModels;

	public PositionNoiseModels()
	{
		noiseModels = new TreeMap<Substitution, NoiseModel>();
	}

	public NoiseModel getModel(Substitution substitution)
	{
		NoiseModel model = noiseModels.get(substitution);

		if (model == null)
		{
			model = new NoiseModel();
			noiseModels.put(substitution, model);
		}

		return model;
	}

}
