package org.mskcc.marianas.umi.duplex.fastqprocessing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.mskcc.marianas.util.Util;

/**
 * 
 * @author Juber Patel
 * 
 *         post process the paired fastq files.
 * 
 *         0. Choose only those read pairs for which both reads have UMIs, for
 *         some definition of having.
 *         1. remove UMIs and constant region from the beginning of read1 and
 *         read2
 *         2. put the UMI combination in the read name for both read1 and read2
 *         3. create new sample folders and new fastq files
 *         4. copy sample sheets
 *
 */
public class ProcessLoopUMIFastq
{

	/**
	 * 
	 * @param args
	 *            read1 fastq, UMI length, constant region, project folder
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		File R1FastqFile = new File(args[0]);
		int UMILength = Integer.parseInt(args[1]);
		File outputProjectFolder = new File(args[2]);

		// no args[] beyond this point

		BufferedReader fastq1 = new BufferedReader(new InputStreamReader(
				new GZIPInputStream(new FileInputStream(R1FastqFile))));
		BufferedReader fastq2 = new BufferedReader(new InputStreamReader(
				new GZIPInputStream(new FileInputStream(new File(
						R1FastqFile.getParentFile(), R1FastqFile.getName()
								.replace("_R1_", "_R2_"))))));
		File sampleSheet = new File(R1FastqFile.getParentFile(),
				"SampleSheet.csv");

		int maxOccurrences = UMILength + 1;

		// make sample folder
		File sampleFolder = sampleFolder(outputProjectFolder, R1FastqFile);
		sampleFolder.mkdir();

		// copy samplesheet
		Files.copy(sampleSheet.toPath(),
				new File(sampleFolder, "SampleSheet.csv").toPath());

		// process fastq files
		String sampleName = sampleFolder.getName().substring(7);
		// make output fastq files
		File outFile1 = new File(sampleFolder, R1FastqFile.getName());
		File outFile2 = new File(sampleFolder,
				R1FastqFile.getName().replace("_R1_", "_R2_"));

		BufferedWriter out1 = new BufferedWriter(new OutputStreamWriter(
				new GZIPOutputStream(new FileOutputStream(outFile1))));
		BufferedWriter out2 = new BufferedWriter(new OutputStreamWriter(
				new GZIPOutputStream(new FileOutputStream(outFile2))));

		// iterate
		String name1 = null;
		String name2 = null;
		String seq1 = null;
		String seq2 = null;
		String garbage1 = null;
		String garbage2 = null;
		String qual1 = null;
		String qual2 = null;
		long[] polyUMICounts = new long[5];
		long totalReadPairs = 0;

		long[][] quals = new long[5][UMILength * 2 + 1];

		while ((name1 = fastq1.readLine()) != null)
		{
			seq1 = fastq1.readLine();
			garbage1 = fastq1.readLine();
			qual1 = fastq1.readLine();

			name2 = fastq2.readLine();
			seq2 = fastq2.readLine();
			garbage2 = fastq2.readLine();
			qual2 = fastq2.readLine();

			totalReadPairs++;

			LoopUMIProcessor readPair = new LoopUMIProcessor(seq1, qual1, seq2,
					qual2, UMILength);

			// invalid read pair, for whatever reason
			if (!(readPair.read1HasUMI() && readPair.read2HasUMI()))
			{
				continue;
			}

			String UMI = readPair.compositeUMI();
			String currentUMIQuals = readPair.compositeUMIQuals();

			if (Util.poly(UMI, 'A', maxOccurrences))
			{
				polyUMICounts[1]++;
				addQualities(quals[1], currentUMIQuals);
			}
			else if (Util.poly(UMI, 'C', maxOccurrences))
			{
				polyUMICounts[2]++;
				addQualities(quals[2], currentUMIQuals);
			}
			else if (Util.poly(UMI, 'G', maxOccurrences))
			{
				polyUMICounts[3]++;
				addQualities(quals[3], currentUMIQuals);
			}
			else if (Util.poly(UMI, 'T', maxOccurrences))
			{
				polyUMICounts[4]++;
				addQualities(quals[4], currentUMIQuals);
			}
			else
			{
				polyUMICounts[0]++;
				addQualities(quals[0], currentUMIQuals);
			}

			// write read1
			String[] parts = name1.split(" ");
			parts[0] = parts[0] + ":" + UMI;
			out1.write(parts[0] + " " + parts[1] + "\n");
			out1.write(readPair.seq1() + "\n");
			out1.write(garbage1 + "\n");
			out1.write(readPair.qual1() + "\n");

			// write read2
			parts = name2.split(" ");
			parts[0] = parts[0] + ":" + UMI;
			out2.write(parts[0] + " " + parts[1] + "\n");
			out2.write(readPair.seq2() + "\n");
			out2.write(garbage2 + "\n");
			out2.write(readPair.qual2() + "\n");
		}

		fastq1.close();
		fastq2.close();
		out1.close();
		out2.close();

		// write basic stats
		BufferedWriter writer = new BufferedWriter(
				new FileWriter(new File(sampleFolder, "info.txt")));

		long readPairsWithUMIs = polyUMICounts[0] + polyUMICounts[1]
				+ polyUMICounts[2] + polyUMICounts[3] + polyUMICounts[4];

		writer.write("Total read pairs: " + totalReadPairs + "\n");
		writer.write("Read pairs with UMIs: " + readPairsWithUMIs + "/"
				+ totalReadPairs + " ("
				+ (readPairsWithUMIs * 1.0 / totalReadPairs) + ")\n");
		writer.write("Poly-A: " + polyUMICounts[1] + "/" + readPairsWithUMIs
				+ " (" + (polyUMICounts[1] * 1.0 / readPairsWithUMIs) + ")\n");
		writer.write("Poly-G: " + polyUMICounts[2] + "/" + readPairsWithUMIs
				+ " (" + (polyUMICounts[2] * 1.0 / readPairsWithUMIs) + ")\n");
		writer.write("Poly-C: " + polyUMICounts[3] + "/" + readPairsWithUMIs
				+ " (" + (polyUMICounts[3] * 1.0 / readPairsWithUMIs) + ")\n");
		writer.write("Poly-T: " + polyUMICounts[4] + "/" + readPairsWithUMIs
				+ " (" + (polyUMICounts[4] * 1.0 / readPairsWithUMIs) + ")\n");
		writer.write("Good: " + polyUMICounts[0] + "/" + readPairsWithUMIs
				+ " (" + (polyUMICounts[0] * 1.0 / readPairsWithUMIs) + ")\n");

		writer.close();

		writer = new BufferedWriter(
				new FileWriter(new File(sampleFolder, "umi-frequencies.txt")));
		Map<String, Integer> map = LoopUMIProcessor.UMIFrequencies;
		for (String key : map.keySet())
		{
			int value = map.get(key);
			writer.write(key + "\t" + value + "\n");
		}

		writer.close();

		writer = new BufferedWriter(new FileWriter(
				new File(sampleFolder, "composite-umi-frequencies.txt")));
		map = LoopUMIProcessor.compositeUMIFrequencies;
		for (String key : map.keySet())
		{
			int value = map.get(key);
			writer.write(key + "\t" + value + "\n");
		}

		writer.close();

		// writeQualities(new File(sampleFolder, sampleName + ".qualities.txt"),
		// quals, UMICounts);
	}

	private static void addQualities(long[] quals, String UMIQuals)
	{
		for (int i = 0; i < UMIQuals.length(); i++)
		{
			quals[i] += (long) UMIQuals.charAt(i);
		}
	}

	private static void writeQualities(File qualitiesFile, long[][] quals,
			long[] UMICounts) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(
				new FileWriter(qualitiesFile));
		writer.write("Position\tAveragePhred\tGroup\n");

		// print position and averages
		for (int i = 0; i < quals[0].length; i++)
		{
			double sum = 0;
			double count = 0;
			for (int j = 0; j < quals.length; j++)
			{
				sum += quals[j][i];
				count += UMICounts[j];
			}

			double phred = (sum * 1.0) / count;

			writer.write((i + 1) + "\t");
			writer.write(phred + "\t");
			writer.write("All\n");
		}

		for (int i = 0; i < quals[0].length; i++)
		{
			for (int j = 0; j < quals.length; j++)
			{
				writer.write((i + 1) + "\t");
				writer.write((quals[j][i] * 1.0) / UMICounts[j] + "\t");
				writer.write(j + "\n");
			}
		}

		writer.close();
	}

	private static File sampleFolder(File projectFolder, File sampleFastq)
	{
		String name = sampleFastq.getParentFile().getName();

		return new File(projectFolder, name);
	}

	private static File sampleFolderOld(File projectFolder, File sampleFastq)
	{
		String[] parts = sampleFastq.getName().split("_");
		StringBuilder sampleFolderName = new StringBuilder("Sample");
		for (int i = 0; i < parts.length - 4; i++)
		{
			sampleFolderName.append("_" + parts[i]);
		}

		return new File(projectFolder, sampleFolderName.toString());
	}

}
