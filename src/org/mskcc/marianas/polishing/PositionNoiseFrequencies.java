/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * @author Juber Patel
 * 
 *         collection of noise models at a genomic position, one for each
 *         possible substitution
 *
 */
public class PositionNoiseFrequencies
{
	private Map<Substitution, NoiseFrequencyCollector> noiseFrequencies;

	public PositionNoiseFrequencies()
	{
		noiseFrequencies = new LinkedHashMap<Substitution, NoiseFrequencyCollector>();
	}

	public NoiseFrequencyCollector getFrequencyCollectorFor(Substitution substitution)
	{
		NoiseFrequencyCollector model = noiseFrequencies.get(substitution);

		if (model == null)
		{
			model = new NoiseFrequencyCollector();
			noiseFrequencies.put(substitution, model);
		}

		return model;
	}

	public Map<Substitution, NoiseFrequencyCollector> getModelsMap()
	{
		return noiseFrequencies;
	}

}
