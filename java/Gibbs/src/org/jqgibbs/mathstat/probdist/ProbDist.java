package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.mathstat.Numeric;

public abstract class ProbDist<T extends Numeric> {
	protected static final boolean CHECK_PARMS = true;
	protected boolean initialized;

	protected void checkInitialized() {
		if (!this.initialized) {
			throw new IllegalStateException(
					"use of uninitialized probability distribution");
		}
	}

	public double density(T pt) {
		this.checkInitialized();
		return this.getDensity(pt);
	}

	protected double getDensity(T pt) {
		return Math.exp(this.getLogDensity(pt));
	}
	
	public double logDensity(T pt) {
		this.checkInitialized();
		return this.getLogDensity(pt);
	}
	
	protected abstract double getLogDensity(T pt);

	public T variate() {
		this.checkInitialized();
		return this.genVariate();
	}
	
	protected abstract T genVariate();
}
