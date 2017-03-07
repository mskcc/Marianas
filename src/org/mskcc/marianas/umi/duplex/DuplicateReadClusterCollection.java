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

	/**
	 * UMI -> cluster map
	 */
	private Map<String, DuplicateReadCluster> clusters;
	private DuplicateReadCluster[] processedClusters;
	private ObjectPool<DuplicateReadCluster> clusterPool;

	public DuplicateReadClusterCollection(
			ObjectPool<DuplicateReadCluster> clusterPool)
	{
		this.clusters = new HashMap<String, DuplicateReadCluster>();
		this.clusterPool = clusterPool;
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
		if (cluster == null)
		{
			// TODO decide if you want to use Apache Pool !!!
			// cluster = clusterPool.borrowObject();
			cluster = new DuplicateReadCluster();
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
