package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.mathstat.AbstractSequence;
import org.jqgibbs.mathstat.Numeric;

public class IIDDist<U extends AbstractSequence<U, T>, T extends Numeric<T>>
		extends ProbDist<U> {

	private ProbDist<T> dist;
	private int N;

	protected void setDist(ProbDist<T> dist) {
		this.dist = dist;
	}

	protected ProbDist<T> getDist() {
		return this.dist;
	}

	protected void setN(int n) {
		N = n;
	}

	protected int getN() {
		return N;
	}

	public IIDDist(Numeric<?>... parms) throws ProbDistParmException {
		// FIXME - this won't work because you'll still call super() and crash
		// consider making each PD a factory for a dist like this!
		throw new UnsupportedOperationException(
				"Cannot initialize this class using numeric parameters");
	}

	public IIDDist(ProbDist<T> dist, int N) throws ProbDistParmException {
		this.setDist(dist);
		this.setN(N);
	}

	@SuppressWarnings("unchecked")
	@Override
	protected U genVariate() throws ProbDistParmException {
		return (U) this.getDist().variatesIID(this.getN());
	}

	@Override
	protected double getDensity(U pt) throws ProbDistParmException {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}

	@Override
	protected void initializeParms(Numeric<?>... parms)
			throws ProbDistParmException {
		this.initialized = true;
	}

	@Override
	protected void initializeParmsInitial(Numeric<?>... parms) throws ProbDistParmException {
		this.initializeParms(parms);
	}

	@Override
	protected void installParmChecks() {
		return;
	}
}
