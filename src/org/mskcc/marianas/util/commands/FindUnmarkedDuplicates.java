/**
 * 
 */
package org.mskcc.marianas.util.commands;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.mskcc.juber.alignment.Fragment;
import org.mskcc.juber.util.JuberUtilException;
import org.mskcc.marianas.util.Constants;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;

/**
 * @author Juber Patel
 *
 */
public class FindUnmarkedDuplicates
{

	public static Map<String, Fragment> loadFragments(File bamFile,
			Interval interval) throws IOException, JuberUtilException
	{
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(bamFile);

		int start = interval.getStart();
		int end = interval.getEnd();

		// adding some padding to catch the reads which are on target or
		// near target
		SAMRecordIterator iterator = reader.query(interval.getContig(),
				start - 50, end + 50, false);

		// initiate fragment map
		Map<String, Fragment> fragments = new TreeMap<String, Fragment>();

		int records = 0;
		while (iterator.hasNext())
		{
			SAMRecord record = iterator.next();

			// read and its mate must be mapped
			if (record.getReadUnmappedFlag() || record.getMateUnmappedFlag())
			{
				continue;
			}

			// ignore the ones already marked duplicate
			if (record.getDuplicateReadFlag())
			{
				continue;
			}

			addToFragments(record, fragments);
			records++;
			if (records % 1000000 == 0)
			{
				System.out.println(records);
			}
		}

		System.out.println("non-duplicate records: " + records);
		System.out.println("non-duplicate fragments: " + fragments.size());

		int completeFragments = 0;
		for (String key : fragments.keySet())
		{
			if (fragments.get(key).hasCompleteInformation())
			{
				completeFragments++;
			}
		}

		System.out.println("non-duplicate fragments with both proper reads: "
				+ completeFragments);

		iterator.close();
		reader.close();
		iterator = null;
		reader = null;

		return fragments;
	}

	private static void addToFragments(SAMRecord record,
			Map<String, Fragment> fragments) throws JuberUtilException
	{
		String readName = record.getReadName();

		Fragment fragment = fragments.get(readName);

		// make new fragment
		if (fragment == null)
		{
			fragment = new Fragment(record);
			fragments.put(readName, fragment);
		}
		else
		{
			fragment.addRightRecord(record);
		}
	}

	private static void findUnmarkedDuplicates(Map<String, Fragment> fragments)
	{
		Map<Fragment, Fragment> map = new TreeMap<Fragment, Fragment>();
		int duplicates = 0;

		for (String key : fragments.keySet())
		{
			Fragment fragment = fragments.get(key);
			if (!fragment.hasCompleteInformation())
			{
				continue;
			}

			Fragment f = map.get(fragment);

			// map already has a fragment that equals fragment
			if (f != null)
			{
				// System.out.println(fragment);
				// System.out.println(f);
				duplicates++;
			}
			else
			{
				map.put(fragment, fragment);
			}
		}

		System.out.println("unmarked duplicate fragments: " + duplicates);
	}

	/**
	 * @param args
	 * @throws JuberUtilException
	 * @throws IOException
	 */
	public static void main(String[] args)
			throws IOException, JuberUtilException
	{

		File bamFile = new File(
				"/Users/patelj1/workspace/Moynahan/FinalBams/1196-2-IGO-05500-AL-21_bc37_5500-AL_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam");
		// ERBB2
		Interval interval = new Interval("17", 37855812, 37884297);

		Map<String, Fragment> fragments = loadFragments(bamFile, interval);

		findUnmarkedDuplicates(fragments);

	}

}
