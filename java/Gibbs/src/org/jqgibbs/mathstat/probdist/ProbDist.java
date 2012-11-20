package org.jqgibbs.mathstat.probdist;

import java.util.LinkedList;
import java.util.List;

import org.jqgibbs.mathstat.AbstractSequence;
import org.jqgibbs.mathstat.ListSequence;
import org.jqgibbs.mathstat.Numeric;

public abstract class ProbDist<T extends Numeric<T>> {
	protected boolean initialized;

	public ProbDist(Numeric<?>... parms) throws ProbDistParmException {
		this.installParmChecks();
		this.initializeParmsInitial(parms);
	}

	protected abstract double getDensity(T pt) throws ProbDistParmException;
	
	protected double getLogDensity(T pt) throws ProbDistParmException {
		return Math.log(this.getDensity(pt));
	}

	protected abstract T genVariate() throws ProbDistParmException;

	private void checkInitialized(Numeric<?>... parms)
			throws ProbDistParmException {
		// FIXME - Should assert initializedInitial somewhere - should probably even test...
		if (parms.length == 0) {
			if (!this.initialized) {
				throw new IllegalStateException(
						"Tried to perform illegal operation with "
								+ "uninitialized ProbDist");
			}
		} else {
			this.initializeParms(parms);
		}
	}

	/*
	 * This is declared as final in order to guarantee that
	 * ProbDistParmException is never thrown on a zero-args call.
	 */
	public final double density(T pt, Numeric<?>... parms)
			throws ProbDistParmException {
		this.checkInitialized(parms);
		return this.getDensity(pt);
	}
	
	public final double logDensity(T pt, Numeric<?>... parms) throws ProbDistParmException {
		this.checkInitialized(parms);
		return this.getLogDensity(pt);
	}

	/*
	 * This is declared as final in order to guarantee that
	 * ProbDistParmException is never thrown on a zero-args call.
	 */
	public final T variate(Numeric<?>... parms) throws ProbDistParmException {
		this.checkInitialized(parms);
		return this.genVariate();
	}

	/*
	 * This is declared as final in order to guarantee that
	 * ProbDistParmException is never thrown on a zero-args call.
	 */
	public final AbstractSequence<? extends AbstractSequence<?, T>, T> variatesIID(
			int n, Numeric<?>... parms) throws ProbDistParmException {
		if (n < 1) {
			return new ListSequence();
		}
		this.checkInitialized(parms);
		List<T> l = new LinkedList<T>();
		for (int i = 0; i < n; i++) {
			l.add(this.genVariate());
		}
		AbstractSequence<? extends AbstractSequence<?, T>, T> s;
		s = l.get(0).sequence();
		if (n > 1) {
			l.remove(0);
			s.addAll(l);
		}
		return s;
	}

	protected abstract void installParmChecks();

	protected abstract void initializeParms(Numeric<?>... parms)
			throws ProbDistParmException;

	protected abstract void initializeParmsInitial(Numeric<?>... parms)
			throws ProbDistParmException;

	protected void checkParm(String parmName,
			Class<? extends Numeric<?>> parmType, Numeric<?> parm,
			ProbDistParmCheck... cs) throws ProbDistParmException {
		if (parmType == null && cs == null) {
			return; // OK
		} else if (parm == null) {
			throw new NullPointerException(
					"Parameter object cannot be null (" + parmName + ")");
		} else if (!parmType.isInstance(parm)) {
			throw new IllegalArgumentException(
					"Parameter object of wrong class: expected "
							+ parmType.toString());
		}
		if (cs != null) {
			for (ProbDistParmCheck c : cs) {
				if (!c.test(parm)) {
					throw new ProbDistParmException("Invalid parameter "
							+ parmName + ": " + c.message());
				}
			}
		}
	}

	public T variateFast(Numeric<?>... parms) throws ProbDistParmException {
		return this.variate(parms);
	}

}
