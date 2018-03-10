/**
 * 
 */
package org.mskcc.marianas.util;

import org.apache.commons.pool2.BasePooledObjectFactory;
import org.apache.commons.pool2.PooledObject;
import org.apache.commons.pool2.impl.DefaultPooledObject;
import org.mskcc.marianas.umi.duplex.DuplicateReadCluster;

/**
 * @author Juber Patel
 *
 */
public class DuplicateReadClusterFacory
		extends BasePooledObjectFactory<DuplicateReadCluster>
{
	private int minConsensusPercent;

	public DuplicateReadClusterFacory(int minConsensusPercent)
	{
		this.minConsensusPercent = minConsensusPercent;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.apache.commons.pool2.BasePooledObjectFactory#create()
	 */
	@Override
	public DuplicateReadCluster create() throws Exception
	{
		return new DuplicateReadCluster(minConsensusPercent);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.apache.commons.pool2.BasePooledObjectFactory#wrap(java.lang.Object)
	 */
	@Override
	public PooledObject<DuplicateReadCluster> wrap(DuplicateReadCluster cluster)
	{
		return new DefaultPooledObject<DuplicateReadCluster>(cluster);
	}
}
