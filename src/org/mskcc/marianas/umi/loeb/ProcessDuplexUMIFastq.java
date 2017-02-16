package org.mskcc.marianas.umi.loeb;

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
public class ProcessDuplexUMIFastq
{

	/**
	 * 
	 * @param args
	 *            read1 fastq, UMI length, constant region, project folder
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		BufferedReader fastq1 = new BufferedReader(new InputStreamReader(
				new GZIPInputStream(new FileInputStream(args[0]))));
		BufferedReader fastq2 = new BufferedReader(
				new InputStreamReader(new GZIPInputStream(
						new FileInputStream(args[0].replace("_R1_", "_R2_")))));
		File sampleSheet = new File(new File(args[0]).getParentFile(),
				"SampleSheet.csv");

		int UMILength = Integer.parseInt(args[1]);
		String constantRegion = args[2];
		// int maxGs = (UMILength / 4) + 1;
		int maxGs = 80;

		// make sample folder
		File projectFolder = new File(args[3]);
		File sampleFolder = sampleFolder(projectFolder, new File(args[0]));
		sampleFolder.mkdir();
		
		// copy samplesheet
		Files.copy(sampleSheet.toPath(),
				new File(sampleFolder, "SampleSheet.csv").toPath());
		
		//process fastq files
		String sampleName = sampleFolder.getName().substring(7);
		// make output fastq files
		File outFile1 = new File(sampleFolder, new File(args[0]).getName());
		File outFile2 = new File(sampleFolder,
				new File(args[0]).getName().replace("_R1_", "_R2_"));

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
		long goodUMIs = 0;
		long polyGUMIs = 0;
		long[] goodQuals = new long[UMILength * 2 + 1];
		long[] polyGQuals = new long[UMILength * 2 + 1];

		while ((name1 = fastq1.readLine()) != null)
		{
			seq1 = fastq1.readLine();
			garbage1 = fastq1.readLine();
			qual1 = fastq1.readLine();

			name2 = fastq2.readLine();
			seq2 = fastq2.readLine();
			garbage2 = fastq2.readLine();
			qual2 = fastq2.readLine();

			DuplexUMIReadPair readPair = new DuplexUMIReadPair(seq1, qual1,
					seq2, qual2, UMILength, constantRegion);

			// invalid read pair, for whatever reason
			if (!(readPair.read1HasUMI() && readPair.read2HasUMI()))
			{
				continue;
			}

			String UMI = readPair.compositeUMI();
			String UMIQuals = readPair.compositeUMIQuals();

			if (Util.polyG(UMI, maxGs))
			{
				polyGUMIs++;
				addQualities(polyGQuals, UMIQuals);
			}
			else
			{
				goodUMIs++;
				addQualities(goodQuals, UMIQuals);
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

		writeQualities(sampleName, polyGQuals, polyGUMIs, goodQuals, goodUMIs);
	}

	private static void addQualities(long[] quals, String UMIQuals)
	{
		for (int i = 0; i < UMIQuals.length(); i++)
		{
			quals[i] += (long) UMIQuals.charAt(i);
		}
	}

	private static void writeQualities(String sampleName, long[] polyGQuals,
			long polyGUMIs, long[] goodQuals, long goodUMIs) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(
				new FileWriter(new File(sampleName + ".qualities.txt")));
		writer.write("Position\tAveragePhred\tGroup\n");

		// print position and averages
		for (int i = 0; i < polyGQuals.length; i++)
		{
			writer.write((i + 1) + "\t");
			double phred = ((goodQuals[i] + polyGQuals[i]) * 1.0)
					/ (goodUMIs + polyGUMIs);
			writer.write(phred + "\t");
			writer.write("All\n");
		}

		for (int i = 0; i < goodQuals.length; i++)
		{
			writer.write((i + 1) + "\t");
			double phred = (goodQuals[i] * 1.0) / goodUMIs;
			writer.write(phred + "\t");
			writer.write("Good\n");
		}

		for (int i = 0; i < polyGQuals.length; i++)
		{
			writer.write((i + 1) + "\t");
			double phred = (polyGQuals[i] * 1.0) / polyGUMIs;
			writer.write(phred + "\t");
			writer.write("PolyG\n");
		}

		writer.close();
	}

	private static File sampleFolder(File projectFolder, File sampleFastq)
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
