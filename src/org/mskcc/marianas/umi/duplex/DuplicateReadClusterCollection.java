/**
 * 
 */
package org.mskcc.marianas.umi.duplex;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.pool2.ObjectPool;

import htsjdk.samtools.SAMRecord;

/**
 * @author Juber Patel
 * 
 *         UMI -> cluster map for all the UMIs with the same genomic start
 *         position ie 5'-most position on the positive strand
 *
 */
public class DuplicateReadClusterCollection
{
	private String contig;
	private int contigIndex;
	private int startPosition;
	private int minMappingQuality;
	private int minBaseQuality;
	private int minConsensusPercent;

	/**
	 * UMI -> cluster map
	 */
	private Map<String, DuplicateReadCluster> clusters;
	private DuplicateReadCluster[] processedClusters;
	private ObjectPool<DuplicateReadCluster> clusterPool;

	public DuplicateReadClusterCollection(
			ObjectPool<DuplicateReadCluster> clusterPool, int minMappingQuality,
			int minBaseQuality, int minConsensusPercent)
	{
		this.clusters = new HashMap<String, DuplicateReadCluster>();
		this.clusterPool = clusterPool;
		this.minMappingQuality = minMappingQuality;
		this.minBaseQuality = minBaseQuality;
		this.minConsensusPercent = minConsensusPercent;
	}

	public void prepareFor(String contig, int contigIndex, int startPosition)
			throws Exception
	{
		this.contig = contig;
		this.contigIndex = contigIndex;
		this.startPosition = startPosition;

		// return the DuplicateReadCluster objects to the pool
		// for (DuplicateReadCluster cluster : clusters.values())
		// {
		// clusterPool.returnObject(cluster);
		// }

		// now clear the map
		this.clusters.clear();
		this.processedClusters = null;
	}

	public void add(String UMI, SAMRecord record, boolean positiveStrand)
			throws Exception
	{
		DuplicateReadCluster cluster = clusters.get(UMI);

		// if (record.getReferenceName().equals("6")
		// && record.getAlignmentStart() == 117725360
		// && UMI.equals("CTC+GCA"))
		// {
		// int a = 5;
		// }

		if (cluster == null)
		{
			// TODO decide if you want to use Apache Pool !!!
			// cluster = clusterPool.borrowObject();
			cluster = new DuplicateReadCluster(minMappingQuality,
					minBaseQuality, minConsensusPercent);
			cluster.prepareFor(contig, startPosition, UMI);
			clusters.put(UMI, cluster);
		}

		try
		{
			cluster.add(record, positiveStrand);
		}
		catch (Exception e)
		{
			System.err.println("Problem processing record:");
			System.err.println(record.getSAMString());
			e.printStackTrace();
		}
	}

	public DuplicateReadCluster[] getProcessedClusters()
	{
		if (processedClusters == null)
		{
			processedClusters = clusters.values()
					.toArray(new DuplicateReadCluster[0]);
		}

		return processedClusters;

	}

	public String getContig()
	{
		return contig;
	}

	public int getContigIndex()
	{
		return contigIndex;
	}

	public int getStartPosition()
	{
		return startPosition;
	}

}
