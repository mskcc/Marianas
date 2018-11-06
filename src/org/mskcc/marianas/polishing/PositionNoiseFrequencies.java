/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * @author Juber Patel
 * 
 *         collection of noise frequencies at a genomic position, one for each
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

	public NoiseFrequencyCollector getFrequencyCollectorFor(
			Substitution substitution)
	{
		NoiseFrequencyCollector freqCollector = noiseFrequencies
				.get(substitution);

		if (freqCollector == null)
		{
			freqCollector = new NoiseFrequencyCollector();
			noiseFrequencies.put(substitution, freqCollector);
		}

		return freqCollector;
	}

	public Map<Substitution, NoiseFrequencyCollector> getModelsMap()
	{
		return noiseFrequencies;
	}

}
