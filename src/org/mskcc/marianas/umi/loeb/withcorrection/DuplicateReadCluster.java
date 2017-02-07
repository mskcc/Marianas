/**
 * 
 */
package org.mskcc.marianas.umi.loeb.withcorrection;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.mskcc.juber.genotype.Genotype;
import org.mskcc.juber.genotype.GenotypeEventType;
import org.mskcc.juber.genotype.GenotypeID;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

/**
 * @author Juber Patel
 * 
 *         class that represents a pileup of reads in a very efficient way
 * 
 */
public class DuplicateReadCluster
{
	private static final IndexedFastaSequenceFile referenceFasta = ProcessDuplexUMIBamFirstPass
			.getReference();
	private static final int maxReadLength = 150;
	private static final int lastPileupIndex = maxReadLength - 1;

	private String contig;
	private int startPosition;
	private String UMI;

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
	private boolean duplicate = false;

	public DuplicateReadCluster()
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

		referenceBases = referenceFasta.getSubsequenceAt(contig, startPosition,
				startPosition + maxReadLength - 1).getBases();

		// clean the pileup for reuse
		for (int i = 0; i < psPositions.length; i++)
		{
			psPositions[i].reset(referenceBases[i]);
			nsPositions[i].reset(referenceBases[i]);
		}

		psSpecialGenotypes.clear();
		nsSpecialGenotypes.clear();
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
	public void merge(DuplicateReadCluster other)
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
	 * @return
	 */
	public String[] consensusSequence()
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
		// record a non-ref position iff it is present in both strands
		for (int i = 0; i < psPositions.length; i++)
		{
			if (psConsensus[i] == nsConsensus[i])
			{
				consensus[i] = psConsensus[i];
			}
			else
			{
				consensus[i] = psPositions[i].getRefBase();
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
		while (true)
		{
			if (consensusSequenceBuilder.charAt(trimmedLength - 1) != 'N')
			{
				break;
			}

			trimmedLength--;
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

		return new String[] { sequence, new String(qualities) };
	}

	/**
	 * 
	 * @return unique name for the collapsed sequence
	 */
	public String consensusSequenceName()
	{
		StringBuilder builder = new StringBuilder("@Marianas:");

		builder.append(contig).append(":").append(startPosition).append(":")
				.append(UMI).append(":").append(psReadCount).append(":")
				.append(nsReadCount);

		return builder.toString();
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
