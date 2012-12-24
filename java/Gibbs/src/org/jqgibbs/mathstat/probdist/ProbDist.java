package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.mathstat.Numeric;

public abstract class ProbDist<T extends Numeric> {
	protected static final boolean CHECK_PARMS = true;
	protected boolean initialized;

	protected abstract double getDensity(T pt);
	
	protected double getLogDensity(T pt) {
		return Math.log(this.getDensity(pt));
	}

	protected abstract T genVariate();
	
	/*
	 * Initializes chain parms for ProbDistInitializeByChain,
	 * Fixed parms for other distributions
	 */
	protected abstract void checkInitialized(Numeric... parms);
	
	/*
	 * This is declared as final in order to guarantee that
	 * ProbDistParmException is never thrown on a zero-args call.
	 */
	public final double density(T pt, Numeric... parms)
			throws ProbDistParmException {
		this.checkInitialized(parms);
		return this.getDensity(pt);
	}
	
	public final double logDensity(T pt) {
		return this.getLogDensity(pt);
	}

	/*
	 * This is declared as final in order to guarantee that
	 * ProbDistParmException is never thrown on a zero-args call.
	 */
	public final T variate(Numeric... parms) {
		this.checkInitialized(parms);
		return this.genVariate();
	}

	public T variateFast(Numeric... parms) {
		return this.variate(parms);
	}

}
