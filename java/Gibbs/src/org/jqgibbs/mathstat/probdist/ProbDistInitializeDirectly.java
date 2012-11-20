package org.jqgibbs.mathstat.probdist;

import java.util.List;

import org.jqgibbs.mathstat.Numeric;

public abstract class ProbDistInitializeDirectly<T extends Numeric<T>> extends ProbDist<T> {

	protected Numeric<?>[] parms;
	
	public ProbDistInitializeDirectly(Numeric<?>... parms)
			throws ProbDistParmException {
		super(parms);
	}

	protected abstract void installParmChecks();
	protected abstract void setUpFromParms() throws ProbDistParmException;
	
	@Override
	/*
	 * Classes must be defined such that parameters unspecified in parameter list are
	 * at the end of getParmNames(), getParmClasses() [as null], and getParmCheck()
	 * [as null].
	 */
	protected void initializeParms(Numeric<?>... parms)
			throws ProbDistParmException {
		// Call all parameter checks and assignments
		this.parms = new Numeric<?>[this.getParmCheck().size()];
		if (parms.length > this.getParmCheck().size()) {
			throw new IllegalArgumentException("Too many parameters");
		} else if (parms.length == 0 && this.getParmCheck().size() > 0) {
			return; // FIXME?? 
		}
		String n;
		Class<? extends Numeric<?>> c;
		Numeric<?> p;
		ProbDistParmCheck[] pcs;
		for (int i=0; i<parms.length; i++) {
			n = this.getParmNames().get(i);
			c = this.getParmClasses().get(i);
			p = parms[i];
			pcs = this.getParmCheck().get(i);
			this.checkParm(n, c, p, pcs);
			if (c != null) {
				this.parms[i] = c.cast(parms[i]);
			} else {
				this.parms[i] = parms[i];
			}
		}
		// Do remaining work
		this.setUpFromParms();
		// Set flag
		this.initialized = true;
	}

	@Override
	protected void initializeParmsInitial(Numeric<?>... parms) throws ProbDistParmException {
		this.initializeParms(parms);
	}
	
	protected abstract List<String> getParmNames();

	protected abstract List<Class<? extends Numeric<?>>> getParmClasses();

	protected abstract List<ProbDistParmCheck[]> getParmCheck();

	@Override
	public T variateFast(Numeric<?>... parms) throws ProbDistParmException {
		this.checkInitializedFast(parms);
		return this.genVariate();
	}

	private void checkInitializedFast(Numeric<?>... parms)
	throws ProbDistParmException {
		if (parms.length == 0) {
			if (!this.initialized) {
				throw new IllegalStateException(
						"Tried to perform illegal operation with "
								+ "uninitialized ProbDist");
			}
		} else {
			this.initializeParmsFast(parms);
		}
	}
	
	protected void initializeParmsFast(Numeric<?>... parms)
	throws ProbDistParmException {
		// Call all parameter checks and assignments
		this.parms = new Numeric<?>[this.getParmCheck().size()];
		String n;
		Class<? extends Numeric<?>> c;
		Numeric<?> p;
		ProbDistParmCheck[] pcs;
		for (int i=0; i<parms.length; i++) {
			p = parms[i];
			this.parms[i] = parms[i];
		}
		// Do remaining work
		this.setUpFromParms();
		// Set flag
		this.initialized = true;
	}	
}
