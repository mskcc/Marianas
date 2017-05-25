package org.mskcc.marianas.commands;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.mskcc.marianas.util.Util;

/**
 * 
 * @author Juber Patel
 * 
 *         calculate edit distance between each pair of barcodes in the given
 *         list
 *
 */
public class FindSimilarUMIs
{
	public static void main(String[] args) throws IOException
	{
		String fileName = "umi/BC12_bc12_Project-5500-AV_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam.umi-summary";
		int mismatches = 2;
		int wobble = 2;

		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		List<String[]> lines = new ArrayList<String[]>();

		String line = null;
		while ((line = reader.readLine()) != null)
		{
			lines.add(line.split("\t"));
		}

		reader.close();

		for (int i = 0; i < lines.size(); i++)
		{
			String[] partsI = lines.get(i);
			int positionI = Integer.parseInt(partsI[2]);
			int chosenClusterCount = 0;
			String[] chosenJ = new String[] { "", "0", "0", "000000000000", "0" };
			for (int j = 0; j < lines.size(); j++)
			{
				if (i == j)
				{
					continue;
				}

				String[] partsJ = lines.get(j);
				int positionJ = Integer.parseInt(partsJ[2]);
				if (!partsI[1].equals(partsJ[1])
						|| Math.abs(positionI - positionJ) > wobble
						|| Util.distance(partsI[3], partsJ[3]) > mismatches)
				{
					continue;
				}

				int clusterCount = Integer.parseInt(partsJ[4]);
				if (clusterCount > chosenClusterCount)
				{
					chosenClusterCount = clusterCount;
					chosenJ = partsJ;
				}
			}

			// print
			StringBuilder builder = new StringBuilder(partsI[0]);
			for (int k = 1; k < partsI.length; k++)
			{
				builder.append("\t").append(partsI[k]);
			}

			builder.append("\t").append(chosenJ[2]);
			builder.append("\t").append(chosenJ[3]);
			builder.append("\t").append(chosenJ[4]);
			builder.append("\t")
					.append(Math.abs(positionI - Integer.parseInt(chosenJ[2])));
			builder.append("\t").append(Util.distance(partsI[3], chosenJ[3]));

			System.out.println(builder.toString());
		}
	}
}
