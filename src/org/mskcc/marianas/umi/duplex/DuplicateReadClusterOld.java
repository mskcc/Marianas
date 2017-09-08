/**
 * 
 */
package org.mskcc.marianas.umi.duplex;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.mskcc.juber.genotype.Genotype;
import org.mskcc.juber.genotype.GenotypeEventType;
import org.mskcc.juber.genotype.GenotypeID;
import org.mskcc.juber.util.Util;
import org.mskcc.marianas.util.StaticResources;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

/**
 * @author Juber Patel
 * 
 *         class that represents a pileup of reads in a very efficient way
 * 
 */
public class DuplicateReadClusterOld
{
	private static final IndexedFastaSequenceFile referenceFasta = StaticResources
			.getReference();
	private static final FastaSequenceIndex referenceFastaIndex = StaticResources
			.getReferenceIndex();

	private static final Map<String, Map<Integer, Byte[]>> genotypes = StaticResources
			.getGenotypes();
	private static final String pos = StaticResources.getPositionOfInterest();

	private static final int maxReadLength = 150;
	private static final int lastPileupIndex = maxReadLength - 1;

	private String contig;
	private int startPosition;
	private String UMI;

	// for debug module
	private String contigOfInterest;
	private int positionOfInterest;
	private List<StringBuilder> linesOfInterest;

	private int psReadCount;
	private int nsReadCount;
	// pileup for fragments coming from positive strand
	private PositionPileup[] psPositions;
	// pileup for fragments coming from negative strand
	private PositionPileup[] nsPositions;
	/**
	 * holds special genotypes: multi-base events and insertions
	 * multi-base substitution not handled yet.
	 */
	private Map<GenotypeID, Genotype> psSpecialGenotypes;
	private Map<GenotypeID, Genotype> nsSpecialGenotypes;

	private byte[] psConsensus;
	private byte[] nsConsensus;
	private byte[] consensus;
	private StringBuilder consensusSequenceBuilder;
	private byte[] referenceBases;

	private byte[] readBases;
	private byte[] baseQualities;
	private int readIndex;

	private Map<String, Integer> matePositions;

	public DuplicateReadClusterOld()
	{
		this.psPositions = new PositionPileup[maxReadLength];
		this.nsPositions = new PositionPileup[maxReadLength];

		// initialize position pileups
		for (int i = 0; i < psPositions.length; i++)
		{
			psPositions[i] = new PositionPileup();
			nsPositions[i] = new PositionPileup();
		}

		// initialize consensii
		this.psConsensus = new byte[maxReadLength];
		this.nsConsensus = new byte[maxReadLength];
		this.consensus = new byte[maxReadLength];
		this.consensusSequenceBuilder = new StringBuilder(maxReadLength + 50);

		psSpecialGenotypes = new HashMap<GenotypeID, Genotype>();
		nsSpecialGenotypes = new HashMap<GenotypeID, Genotype>();

		matePositions = new HashMap<String, Integer>();

		// if running in debug mode
		if (pos != null)
		{
			String[] words = pos.split(":");
			this.contigOfInterest = words[0];
			this.positionOfInterest = Integer.parseInt(words[1].split("-")[0]);
			this.linesOfInterest = new ArrayList<StringBuilder>();
		}
	}

	/**
	 * must call this method before using the region pileup
	 * 
	 * @param interval
	 */
	public void prepareFor(String contig, int startPosition, String UMI)
	{
		this.contig = contig;
		this.startPosition = startPosition;
		this.UMI = UMI;
		this.psReadCount = 0;
		this.nsReadCount = 0;

		int contigLength = referenceFasta.getSequenceDictionary()
				.getSequence(contig).getSequenceLength();
		int endPosition = startPosition + maxReadLength - 1;
		if (endPosition > contigLength)
		{
			endPosition = contigLength;
		}

		referenceBases = referenceFasta
				.getSubsequenceAt(contig, startPosition, endPosition)
				.getBases();

		// clean the pileup for reuse
		for (int i = 0; i < referenceBases.length; i++)
		{
			psPositions[i].reset(referenceBases[i]);
			nsPositions[i].reset(referenceBases[i]);
		}

		// in case referenceBases has shorter length than psPositions and
		// nsPositions.
		// This could happen at the end of the contig (see above adjustment)
		for (int i = referenceBases.length; i < psPositions.length; i++)
		{
			psPositions[i].reset((byte) 'N');
			nsPositions[i].reset((byte) 'N');
		}

		psSpecialGenotypes.clear();
		nsSpecialGenotypes.clear();

		matePositions.clear();
	}

	/**
	 * process the alignment record
	 * 
	 * 
	 * @param record
	 */
	public void add(SAMRecord record, boolean positiveStrand)
	{
		// example of CIGAR string: 146M3I4M5S
		// Assuming that CIGAR will only have these operators: M I D S H

		// TODO HANDLE NEGATIVE STRAND ??? (this is from original RegionPileup)

		PositionPileup[] positions = null;
		Map<GenotypeID, Genotype> specialGenotypes = null;

		recordMatePosition(record);

		if (positiveStrand)
		{
			psReadCount++;
			positions = psPositions;
			specialGenotypes = psSpecialGenotypes;
		}
		else
		{
			nsReadCount++;
			positions = nsPositions;
			specialGenotypes = nsSpecialGenotypes;
		}

		readBases = record.getReadBases();
		// tracks the current position in the read
		readIndex = 0;
		// points to the current position in the pileup
		int pileupIndex = record.getStart() - startPosition;
		int operatorLength = 0;
		Cigar cigar = record.getCigar();
		List<CigarElement> elements = cigar.getCigarElements();

		for (int i = 0; i < elements.size(); i++)
		{
			CigarElement e = elements.get(i);
			CigarOperator operator = e.getOperator();
			operatorLength = e.getLength();

			if (operator.equals(CigarOperator.MATCH_OR_MISMATCH))
			{
				// add the bases
				for (int j = 0; j < operatorLength; j++)
				{
					if (pileupIndex >= 0 && pileupIndex <= lastPileupIndex)
					{
						positions[pileupIndex].addBase(readBases[readIndex]);

						// if running in debug mode and at the right position
						if (contigOfInterest != null
								&& contigOfInterest.equals(record.getContig())
								&& positionOfInterest == startPosition
										+ pileupIndex)
						{
							StringBuilder line = new StringBuilder(
									record.getReadName());
							line.append("\t")
									.append((char) readBases[readIndex]);
							line.append("\t").append(contig).append(":")
									.append(startPosition).append(":")
									.append(UMI);
							linesOfInterest.add(line);
						}
					}

					// increment both pileup index and read index
					pileupIndex++;
					readIndex++;
				}
			}
			else if (operator.equals(CigarOperator.INSERTION))
			{
				// TODO replace this boundary check with proper tracking of
				// CIGAR and quitting when it goes out of the region
				if (pileupIndex >= 1 && pileupIndex <= lastPileupIndex + 1)
				{
					positions[pileupIndex].addInsertion();

					// add insertion to special genotypes map
					// make genotype id
					int precedingGenomicPosition = startPosition
							+ (pileupIndex - 1);
					byte[] ref = new byte[] { referenceBases[pileupIndex - 1] };
					byte[] alt = new byte[operatorLength + 1];
					alt[0] = ref[0];
					System.arraycopy(readBases, readIndex, alt, 1,
							operatorLength);
					// copy(readBases, readIndex - 1,
					// readIndex + operatorLength);
					GenotypeID genotypeID = new GenotypeID(
							GenotypeEventType.INSERTION, contig,
							precedingGenomicPosition, ref, alt);
					// add
					addSpecialGenotype(specialGenotypes, genotypeID);
				}

				// increment readIndex but not PileupIndex
				readIndex += operatorLength;
			}
			else if (operator.equals(CigarOperator.DELETION))
			{
				// add deletion to the special genotypes map iff it is a
				// multi-base deletion
				if (operatorLength > 1 && pileupIndex >= 1
						&& pileupIndex <= lastPileupIndex + 1)
				{
					// make genotype id
					int precedingGenomicPosition = startPosition
							+ (pileupIndex - 1);
					byte[] alt = new byte[] { referenceBases[pileupIndex - 1] };
					byte[] ref = referenceFasta
							.getSubsequenceAt(contig, precedingGenomicPosition,
									precedingGenomicPosition + operatorLength)
							.getBases();
					// byte[] ref = Arrays.copyOfRange(referenceBases,
					// pileupIndex - 1, pileupIndex + operatorLength);

					GenotypeID genotypeID = new GenotypeID(
							GenotypeEventType.DELETION, contig,
							precedingGenomicPosition, ref, alt);
					// TODO only adding insertions
					// TODO confirm that's the right move. Then delete this code
					// block.
					// addSpecialGenotype(specialGenotypes, genotypeID);
				}

				// add deletions to the pileup
				for (int j = 0; j < operatorLength; j++)
				{
					if (pileupIndex >= 0 && pileupIndex <= lastPileupIndex)
					{
						positions[pileupIndex].addDeletion();
					}

					// increment pileupIndex but don't increment readIndex
					pileupIndex++;
				}
			}
			else if (operator.equals(CigarOperator.SOFT_CLIP))
			{
				// soft clips do not count in alignment length, no change in
				// pileupIndex

				// but increment the read index since those bases are not
				// clipped from the read
				readIndex += e.getLength();
			}
			else if (operator.equals(CigarOperator.HARD_CLIP))
			{
				// hard clips do not count in alignment length, nor are the
				// clipped bases included in the read sequence
				// so no change in either pileupIndex or readIndex
			}
			else
			{
				// increment both pileupIndex and readIndex properly
				pileupIndex += operatorLength;
				readIndex += operatorLength;
			}
		}
	}

	private void recordMatePosition(SAMRecord record)
	{
		String matePosition = record.getMateReferenceIndex() + "\t"
				+ record.getMateReferenceName() + "\t"
				+ record.getMateAlignmentStart();

		Integer freq = matePositions.get(matePosition);
		if (freq == null)
		{
			matePositions.put(matePosition, 1);
		}
		else
		{
			matePositions.put(matePosition, freq + 1);
		}
	}

	/**
	 * 
	 * @param genotypeID
	 * 
	 *            increment the coverage count of the given genotype by 1
	 */
	private void addSpecialGenotype(Map<GenotypeID, Genotype> specialGenotypes,
			GenotypeID genotypeID)
	{
		Genotype genotype = specialGenotypes.get(genotypeID);

		if (genotype == null)
		{
			genotype = new Genotype("FromPileup", genotypeID);
			specialGenotypes.put(genotypeID, genotype);
		}

		genotype.totalSupportingCoverage++;
	}

	/**
	 * merge the other RegionPileup with this one
	 * 
	 * @param other
	 */
	public void merge(DuplicateReadClusterOld other)
	{
		psReadCount += other.psReadCount;
		nsReadCount += other.nsReadCount;

		// merge position pileups
		for (int i = 0; i < psPositions.length; i++)
		{
			psPositions[i].merge(other.psPositions[i]);
			nsPositions[i].merge(other.nsPositions[i]);
		}

		// merge special genotypes
		for (GenotypeID id : other.psSpecialGenotypes.keySet())
		{
			Genotype otherGenotype = other.psSpecialGenotypes.get(id);
			Genotype genotype = psSpecialGenotypes.get(id);
			if (genotype == null)
			{
				psSpecialGenotypes.put(id, otherGenotype);
			}
			else
			{
				genotype.totalSupportingCoverage += otherGenotype.totalSupportingCoverage;
			}
		}

		for (GenotypeID id : other.nsSpecialGenotypes.keySet())
		{
			Genotype otherGenotype = other.nsSpecialGenotypes.get(id);
			Genotype genotype = nsSpecialGenotypes.get(id);
			if (genotype == null)
			{
				nsSpecialGenotypes.put(id, otherGenotype);
			}
			else
			{
				genotype.totalSupportingCoverage += otherGenotype.totalSupportingCoverage;
			}
		}
	}

	/**
	 * build consensus sequence from the pileup. This involves building
	 * consensus sequences per strand first and then building the final cluster
	 * consensus sequence.
	 * 
	 * If running in debug mode, write to the writer and return null. Otherwise
	 * writer is null; return the consensus info
	 * 
	 * @param positiveStrand
	 *            does the consensus read map on the positive strand or negative
	 *            strand?
	 * @return
	 * @throws IOException
	 */
	public String consensusSequenceInfo(BufferedWriter writer,
			boolean positiveStrand) throws IOException
	{
		// TODO This method will undergo revision as we determine how exactly we
		// want to build the consensus sequence. There are many parameters to
		// play with.
		///////////////////

		// clear the consensii arrays and StringBuilder
		for (int i = 0; i < psConsensus.length; i++)
		{
			psConsensus[i] = -1;
		}

		System.arraycopy(psConsensus, 0, nsConsensus, 0, psConsensus.length);
		System.arraycopy(psConsensus, 0, consensus, 0, psConsensus.length);
		consensusSequenceBuilder.setLength(0);

		// compute per strand consensus sequence. Right now, just choosing the
		// base with highest count at each position.
		for (int i = 0; i < psPositions.length; i++)
		{
			psConsensus[i] = psPositions[i].getMaxCountBase();
			nsConsensus[i] = nsPositions[i].getMaxCountBase();
		}

		// build cluster consensus sequence.
		
		// record a base that is not person's genotype iff it is present in both
		// strands
		for (int i = 0; i < psPositions.length; i++)
		{
			if (psConsensus[i] == nsConsensus[i])
			{
				consensus[i] = psConsensus[i];
			}
			else
			{
				// find this person's genotype. If not available, use ref base
				Byte[] genotype = null;
				Map<Integer, Byte[]> chrMap = genotypes.get(contig);
				if (chrMap != null)
				{
					genotype = chrMap.get(startPosition + i);
				}

				if (genotype == null || genotype.length == 0)
				{
					genotype = new Byte[] { psPositions[i].getRefBase() };
				}

				// if either strand consensus is part of person's genotype,
				// accept it.

				// TODO Assuming it is impossible to have 2 different bases as
				// consensus where those bases are person's authentic genotype
				int j;
				for (j = 0; j < genotype.length; j++)
				{
					byte g = genotype[j];
					if (psConsensus[i] == g || nsConsensus[i] == g)
					{
						consensus[i] = g;
						break;
					}
				}

				// could not assign genotype so far, choose the first from
				// genotype array
				if (j == genotype.length)
				{
					consensus[i] = genotype[0];
				}

			}

			// what if one strand has alt and the other strand has N??
		}

		// build the consensus sequence
		for (int i = 0; i < consensus.length; i++)
		{
			byte b = consensus[i];
			if (b == -1)
			{
				b = 'N';
			}

			consensusSequenceBuilder.append((char) b);
		}

		// if running in debug mode
		if (contigOfInterest != null && contigOfInterest.equals(contig))
		{
			int baseOfInterestIndex = positionOfInterest - startPosition;

			byte b = psConsensus[baseOfInterestIndex];
			char psConsensusBase = (b == -1 ? 'N' : (char) b);

			b = nsConsensus[baseOfInterestIndex];
			char nsConsensusBase = (b == -1 ? 'N' : (char) b);

			char consensusBase = consensusSequenceBuilder
					.charAt(baseOfInterestIndex);

			// print the read lines
			for (StringBuilder line : linesOfInterest)
			{
				line.append("\t").append(psConsensusBase).append("\t")
						.append(psReadCount).append("\t")
						.append(nsConsensusBase).append("\t")
						.append(nsReadCount).append("\t").append(consensusBase)
						.append("\n");
				writer.write(line.toString());
			}

			return null;
		}

		// Include insertions
		for (GenotypeID genotypeID : psSpecialGenotypes.keySet())
		{
			Genotype negative = nsSpecialGenotypes.get(genotypeID);

			// must be present on both strands, must have at least 50% support
			// on each strand
			if (negative == null
					|| negative.totalSupportingCoverage < nsReadCount / 2
					|| psSpecialGenotypes.get(
							genotypeID).totalSupportingCoverage < psReadCount)
			{
				continue;
			}

			// insert the insertion
			int index = genotypeID.position - startPosition;
			char[] bases = new char[genotypeID.alt.length];
			for (int i = 0; i < bases.length; i++)
			{
				bases[i] = (char) genotypeID.alt[i];
			}

			consensusSequenceBuilder.insert(index, bases);
		}

		// trim trailing N's
		int trimmedLength = consensusSequenceBuilder.length();
		while (trimmedLength > 0)
		{
			if (consensusSequenceBuilder.charAt(trimmedLength - 1) != 'N')
			{
				break;
			}

			trimmedLength--;
		}
		
		// the consensus is all N's
		if(trimmedLength == 0)
		{
			return null;
		}

		consensusSequenceBuilder.setLength(trimmedLength);

		// remove deletion bases
		String sequence = consensusSequenceBuilder.toString().replace("D", "");

		// build consensus sequence qualities
		// right now, returning constant quality of 50
		// phred 50 + 33 = 83 = 'S'
		// this could be average quality of consensus bases or some other scheme
		// in future
		char[] qualities = new char[sequence.length()];
		for (int i = 0; i < qualities.length; i++)
		{
			qualities[i] = 'S';
		}

		// reverse the strand if necessary
		if (!positiveStrand)
		{
			sequence = Util.reverseComplement(sequence);
			Util.reverse(qualities);
		}

		StringBuilder info = new StringBuilder();

		info.append(contig).append("\t").append(startPosition).append("\t")
				.append(UMI).append("\t").append(psReadCount).append("\t")
				.append(nsReadCount).append("\t").append(getMatePosition())
				.append("\t").append(sequence).append("\t").append(qualities);

		return info.toString();
	}

	private String getMatePosition()
	{
		int max = 0;
		String position = null;
		for (String key : matePositions.keySet())
		{
			int freq = matePositions.get(key);
			if (freq > max)
			{
				max = freq;
				position = key;
			}
		}

		return position;
	}

	public String getContig()
	{
		return contig;
	}

	public int getStartPosition()
	{
		return startPosition;
	}

	public String getUMI()
	{
		return UMI;
	}

	public int getTotalCount()
	{
		return psReadCount + nsReadCount;
	}
}
