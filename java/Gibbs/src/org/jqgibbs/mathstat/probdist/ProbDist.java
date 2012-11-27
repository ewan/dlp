package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.jqgibbs.mathstat.Numeric;

public abstract class ProbDist<T extends Numeric> {
	protected static final boolean CHECK_PARMS = true;
	protected boolean initialized;

	protected abstract double getDensity(T pt) throws ProbDistParmException;
	
	protected double getLogDensity(T pt) throws ProbDistParmException {
		return Math.log(this.getDensity(pt));
	}

	protected abstract T genVariate() throws ProbDistParmException;
	
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
	
	public final double logDensity(T pt) throws ProbDistParmException {
		return this.getLogDensity(pt);
	}

	/*
	 * This is declared as final in order to guarantee that
	 * ProbDistParmException is never thrown on a zero-args call.
	 */
	public final T variate(Numeric... parms) throws ProbDistParmException {
		this.checkInitialized(parms);
		return this.genVariate();
	}

	/*
	 * This is declared as final in order to guarantee that
	 * ProbDistParmException is never thrown on a zero-args call.
	 */
	@SuppressWarnings("unchecked")
	public final List<T> variatesIID(
			int n, Numeric... parms) throws ProbDistParmException {
		if (n < 1) {
			return new ArrayList<T>();
		}
		this.checkInitialized(parms);
		List<T> l = new LinkedList<T>();
		for (int i = 0; i < n; i++) {
			l.add(this.genVariate());
		}
		List<T> s = new ArrayList<T>();
		for(Numeric elem : l.get(0).rowVec()) {
			s.add((T)elem);
		}
		if (n > 1) {
			l.remove(0);
			s.addAll(l);
		}
		return s;
	}

	public T variateFast(Numeric... parms) throws ProbDistParmException {
		return this.variate(parms);
	}

}
