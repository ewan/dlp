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
		if(this.data != null) return;
		this.data = (Double2D)parms[1];
		this.initializeChainParms((ChainLink)parms[0]);
		this.setUpFromChainParms();
	}
	
}

//package org.jqgibbs.mathstat.probdist;
//
//import java.util.List;
//
//import org.jqgibbs.ChainLink;
//import org.jqgibbs.mathstat.Double2D;
//import org.jqgibbs.mathstat.Numeric;
//import org.jqgibbs.mathstat.RandomVar;
//
//public abstract class ProbDistInitializeByChain<T extends Numeric> extends
//		ProbDist<T> {
//
//	private Double2D data;
//	private ChainLink l;
//
//	protected Numeric[] chainParms;
//	protected List<String> chainParmNames;
//	protected List<Class<? extends Numeric>> chainParmClasses;
//	protected List<ProbDistParmCheck[]> chainParmCheck;
//	protected boolean initializedChain; // Should probably fold these two
//										// together somehow better than through
//										// initializeParms
//
//	protected Double2D getSamplerData() {
//		return this.data;
//	}
//
//	private void setSamplerData(Double2D data) {
//		this.data = data;
//	}
//	
//	protected ChainLink getChainLink() {
//		return this.l;
//	}
//	
//	private void setChainLink(ChainLink l) {
//		this.l = l;
//	}
//
//	protected abstract void installParmChecks();
//
//	protected abstract void setUpFromChainParms();
//
//	/*
//	 * Must be exactly the right length.
//	 */
//	private void initializeChainParms(ChainLink l, Double2D d)
//			throws ProbDistParmException {
//		this.chainParms = new Numeric[this.chainParmCheck.size()];
//		for (int i = 0; i < this.chainParmCheck.size(); i++) {
//			if (!(l.contains(this.chainParmNames.get(i)))) {
//				throw new IllegalArgumentException(
//						"Tried to initialize ProbDist from chain missing variable "
//								+ this.chainParmNames.get(i));
//			}
//			RandomVar<?> rv = (RandomVar<?>) l.get(this.chainParmNames.get(i));
//			//this.checkParm(this.chainParmNames.get(i), this.chainParmClasses
//			//		.get(i), (Numeric) rv.getNumericValue(), this.chainParmCheck.get(i));
//			this.chainParms[i] = (Numeric) rv.getNumericValue();
//		}
//		this.setChainLink(l);
//		this.setSamplerData(d);
//		this.setUpFromChainParms();
//		this.initializedChain = true;
//	}
//
//	@Override
//	public T variateFast(Numeric... parms) throws ProbDistParmException {
//		this.checkInitializedFast(parms);		
//		return this.genVariate();
//	}
//	
//	private void checkInitializedFast(Numeric... parms)
//	throws ProbDistParmException {
//		if (parms.length == 0) {
//			if (!this.initialized) {
//				throw new IllegalStateException(
//						"Tried to perform illegal operation with "
//								+ "uninitialized ProbDist");
//			}
//		} else {
//			this.initializeParmsFast(parms);
//		}
//	}
//
//	protected void initializeParmsFast(Numeric... parms)
//	throws ProbDistParmException {
//		ChainLink l = (ChainLink) parms[0];
//		Double2D d = (Double2D) parms[1];
//		this.initializeChainParmsFast(l, d);
//		assert this.initializedChain;
//		this.initialized = true;
//	}
//
//	private void initializeChainParmsFast(ChainLink l, Double2D d)
//	throws ProbDistParmException {
//		List<ProbDistParmCheck[]> x = this.chainParmCheck;
//		this.chainParms = new Numeric[this.chainParmCheck.size()];
//		for (int i = 0; i < this.chainParmCheck.size(); i++) {
//			RandomVar<?> rv = (RandomVar<?>) l.get(this.chainParmNames.get(i));
//			this.chainParms[i] = (Numeric) rv.getNumericValue();
//		}
//		this.setChainLink(l);
//		this.setSamplerData(d);
//		this.setUpFromChainParms();
//		this.initializedChain = true;
//	}
//	
//	protected void initializeParms(Numeric... parms)
//			throws ProbDistParmException {
//		if (!(parms.length == 2 && parms[0] instanceof ChainLink && parms[1] instanceof Double2D)) {
//			throw new IllegalArgumentException(
//					"ProbDist initialized by chain must be initialized "
//							+ "using a ChainLink and a Double2D object");
//		}
//		ChainLink l = (ChainLink) parms[0];
//		Double2D d = (Double2D) parms[1];
//		this.initializeChainParms(l, d);
//		assert this.initializedChain;
//		this.initialized = true;
//	}
//}
