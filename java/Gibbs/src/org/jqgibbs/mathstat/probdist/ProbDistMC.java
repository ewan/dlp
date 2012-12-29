package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.ChainLink;
import org.jqgibbs.Flattenable;

public abstract class ProbDistMC<T extends Flattenable> extends ProbDist<T> {
	
	public abstract void setMCState(ChainLink l);
	
	public T variate(ChainLink l) {
		this.setMCState(l);
		return this.genVariate();
	}
	
	@Override
	protected double getLogDensity(T x) {
		// FIXME
		throw new UnsupportedOperationException("Too lazy, come back later");
	}
}
