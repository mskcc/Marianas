package org.mskcc.marianas.util.commands;

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
public class CalculateEditDistance
{
	public static void main(String[] args) throws IOException
	{
		List<String> seqs = new ArrayList<String>();
		List<String> names = new ArrayList<String>();

		BufferedReader reader = new BufferedReader(
				new FileReader("barcodeKey96.txt"));
		System.out.println("Barcode1\tName1\tBarcode2\tName2\tDistance");

		String line = null;
		while ((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");

			// compare current seq with all the seqs that are already in the
			// list
			for (int i = 0; i < seqs.size(); i++)
			{
				int distance = Util.distance(seqs.get(i), parts[0]);
				if (distance < 3)
				{
					System.out.println(
							parts[0] + "\t" + parts[1] + "\t" + seqs.get(i)
									+ "\t" + names.get(i) + "\t" + distance);
				}
			}

			// add current seq to the list
			seqs.add(parts[0]);
			names.add(parts[1]);
		}

		reader.close();
	}
}
