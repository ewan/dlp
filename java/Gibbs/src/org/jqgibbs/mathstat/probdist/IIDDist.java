package org.jqgibbs.mathstat.probdist;

//import org.jqgibbs.mathstat.AbstractSequence;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Numeric;

public class IIDDist<U extends Numeric, T extends Numeric>
		extends ProbDist<U> {

	private ProbDist<T> dist;
	private int N;
	
	
	// TODO: this seems wrong
	// does variate() ever get called on IIDDist with parameters?
	protected void checkInitialized(Numeric... parms) {
		if(parms.length == 0) return;
		this.setDist((ProbDist<T>)parms[0]);
		this.setN(((Integer0D)parms[1]).value());
	}

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

	public IIDDist(Numeric... parms) throws ProbDistParmException {
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
}
