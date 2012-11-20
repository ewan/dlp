package org.jqgibbs.mathstat.probdist;

import java.util.List;

import org.jqgibbs.ChainLink;
import org.jqgibbs.SamplerData;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.RandomVar;

public abstract class ProbDistInitializeByChain<T extends Numeric<T>> extends
		ProbDist<T> {

	private Double2D data;
	private ChainLink l;

	protected Numeric<?>[] fixedParms;
	protected List<String> fixedParmNames;
	protected List<Class<? extends Numeric<?>>> fixedParmClasses;
	protected List<ProbDistParmCheck[]> fixedParmCheck;
	protected boolean initializedFixed; // FIXME - see below

	protected Numeric<?>[] chainParms;
	protected List<String> chainParmNames;
	protected List<Class<? extends Numeric<?>>> chainParmClasses;
	protected List<ProbDistParmCheck[]> chainParmCheck;
	protected boolean initializedChain; // Should probably fold these two
										// together somehow better than through
										// initializeParms

	public ProbDistInitializeByChain(Numeric<?>... fixed)
			throws ProbDistParmException {
		super(fixed);
	}

	protected Double2D getSamplerData() {
		return this.data;
	}

	private void setSamplerData(Double2D data) {
		this.data = data;
	}
	
	protected ChainLink getChainLink() {
		return this.l;
	}
	
	private void setChainLink(ChainLink l) {
		this.l = l;
	}

	/*
	 * Must be exactly the right length.
	 */
	protected void initializeParmsInitial(Numeric<?>... fixed)
			throws ProbDistParmException {
		this.fixedParms = new Numeric<?>[this.fixedParmCheck.size()];
		if (!(fixed.length == this.fixedParmCheck.size())) {
			throw new IllegalArgumentException(
					"Wrong number of fixed parameters");
		}
		String s;
		Class<? extends Numeric<?>> c;
		Numeric<?> p;
		ProbDistParmCheck[] pcs;
		for (int i = 0; i < this.fixedParmCheck.size(); i++) {
			s = this.fixedParmNames.get(i);
			c = this.fixedParmClasses.get(i);
			p = fixed[i];
			pcs = this.fixedParmCheck.get(i);
			this.checkParm(s, c, p, pcs);
			this.fixedParms[i] = c.cast(fixed[i]);
		}
		this.initializedFixed = true;
	}

	protected abstract void installParmChecks();

	protected abstract void setUpFromChainParms();

	/*
	 * Must be exactly the right length.
	 */
	private void initializeChainParms(ChainLink l, Double2D d)
			throws ProbDistParmException {
		this.chainParms = new Numeric<?>[this.chainParmCheck.size()];
		for (int i = 0; i < this.chainParmCheck.size(); i++) {
			if (!(l.contains(this.chainParmNames.get(i)))) {
				throw new IllegalArgumentException(
						"Tried to initialize ProbDist from chain missing variable "
								+ this.chainParmNames.get(i));
			}
			RandomVar<?> rv = (RandomVar<?>) l.get(this.chainParmNames.get(i));
			this.checkParm(this.chainParmNames.get(i), this.chainParmClasses
					.get(i), rv.getNumericValue(), this.chainParmCheck.get(i));
			this.chainParms[i] = rv.getNumericValue();
		}
		this.setChainLink(l);
		this.setSamplerData(d);
		this.setUpFromChainParms();
		this.initializedChain = true;
	}

	@Override
	protected void initializeParms(Numeric<?>... parms)
			throws ProbDistParmException {
		assert this.initializedFixed;
		if (!(parms.length == 2 && parms[0] instanceof ChainLink && parms[1] instanceof Double2D)) {
			throw new IllegalArgumentException(
					"ProbDist initialized by chain must be initialized "
							+ "using a ChainLink and a Double2D object");
		}
		ChainLink l = (ChainLink) parms[0];
		Double2D d = (Double2D) parms[1];
		this.initializeChainParms(l, d);
		assert this.initializedChain;
		this.initialized = true;
	}

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
		assert this.initializedFixed;
		ChainLink l = (ChainLink) parms[0];
		Double2D d = (Double2D) parms[1];
		this.initializeChainParmsFast(l, d);
		assert this.initializedChain;
		this.initialized = true;
	}

	private void initializeChainParmsFast(ChainLink l, Double2D d)
	throws ProbDistParmException {
		this.chainParms = new Numeric<?>[this.chainParmCheck.size()];
		for (int i = 0; i < this.chainParmCheck.size(); i++) {
			RandomVar<?> rv = (RandomVar<?>) l.get(this.chainParmNames.get(i));
			this.chainParms[i] = rv.getNumericValue();
		}
		this.setChainLink(l);
		this.setSamplerData(d);
		this.setUpFromChainParms();
		this.initializedChain = true;
	}	
}
