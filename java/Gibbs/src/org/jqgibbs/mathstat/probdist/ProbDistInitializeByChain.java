package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.ChainLink;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Numeric;

public abstract class ProbDistInitializeByChain<T extends Numeric> extends ProbDist<T> {
	private Double2D data;
	
	protected Double2D getSamplerData() {
		return this.data;
	}
	
	protected abstract void initializeChainParms(ChainLink l);
	protected abstract void setUpFromChainParms();

	@Override
	protected void checkInitialized(Numeric... parms) {
		if(parms.length == 0) return;
		this.data = (Double2D)parms[1];
		this.initializeChainParms((ChainLink)parms[0]);
		this.setUpFromChainParms();
	}
	
}
