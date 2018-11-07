/**
 * 
 */
package org.mskcc.marianas.polishing;

/**
 * @author Juber Patel
 * 
 *         Position-substitution specific ID for a noise model
 *
 */
public class FreqID
{
	public final String chr;
	public final int position;
	public final char ref;
	public final char alt;

	public FreqID(String chr, int position, char ref, char alt)
	{
		this.chr = chr;
		this.position = position;
		this.ref = ref;
		this.alt = alt;
	}

	public String toString()
	{
		return chr + "\t" + position + "\t" + ref + "\t" + alt;
	}

	@Override
	public int hashCode()
	{
		final int prime = 31;
		int result = 1;
		result = prime * result + alt;
		result = prime * result + ((chr == null) ? 0 : chr.hashCode());
		result = prime * result + position;
		result = prime * result + ref;
		return result;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
		{
			return true;
		}

		if (obj == null)
		{
			return false;
		}

		if (getClass() != obj.getClass())
		{
			return false;
		}

		FreqID other = (FreqID) obj;
		if (alt != other.alt)
		{
			return false;
		}

		if (chr == null)
		{
			if (other.chr != null)
			{
				return false;
			}
		}
		else if (!chr.equals(other.chr))
		{
			return false;
		}

		if (position != other.position)
		{
			return false;
		}

		if (ref != other.ref)
		{
			return false;
		}

		return true;
	}

}
