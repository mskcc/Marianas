/**
 * 
 */
package org.mskcc.marianas.polishing;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Juber Patel
 * 
 *         base substitutions
 *
 */
public class Substitution implements Comparable<Substitution>
{
	private static final Map<Character, Map<Character, Substitution>> lookup = //
			new HashMap<Character, Map<Character, Substitution>>();

	private char ref;
	private char alt;

	Substitution(char ref, char alt)
	{
		this.ref = ref;
		this.alt = alt;

		addToLookup();
	}

	private void addToLookup()
	{
		Map<Character, Substitution> l2 = lookup.get(ref);
		if (l2 == null)
		{
			l2 = new HashMap<Character, Substitution>();
			lookup.put(ref, l2);
		}

		Substitution sub = l2.get(alt);
		if (sub == null)
		{
			sub = this;
			l2.put(alt, sub);
		}
	}

	public char getRef()
	{
		return ref;
	}

	public char getAlt()
	{
		return alt;
	}

	public static Substitution get(char ref, char alt)
	{
		Map<Character, Substitution> l2 = lookup.get(ref);
		Substitution sub = null;

		if (l2 != null)
		{
			sub = l2.get(alt);
		}

		if (sub == null || l2 == null)
		{
			sub = new Substitution(ref, alt);
		}

		return sub;

	}

	@Override
	public int compareTo(Substitution o)
	{
		if (ref != o.ref)
		{
			return ref - o.ref;
		}

		return alt - o.alt;
	}

	public String toString()
	{
		return ref + "to" + alt;
	}

}
