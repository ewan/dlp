package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.random.Empirical;
import cern.jet.random.EmpiricalWalker;

public class CategoricalDist extends ProbDist<Integer0D> {

	private EmpiricalWalker empiricalGen;
	protected Integer0D k;
	protected Double1D p;

	public CategoricalDist(Double1D p, boolean checkParms) {
		this.setParms(p, checkParms);
	}

	public CategoricalDist(Double1D p) {
		this(p, CHECK_PARMS);
	}

	@Override
	protected void checkInitialized(Numeric... parms) {
		if (parms.length < 1) {
			if (!this.initialized) {
				throw new IllegalStateException(
						"use of uninitialized probability distribution");
			}
		} else {
			this.setParms((Double1D) parms[0]);
		}
	}

	public void setParms(Double1D p, boolean checkParms) {
		Integer0D k = new Integer0D(p.size());
		this.k = k;
		this.p = p;
		this.setUpFromParms(checkParms);
		this.initialized = true;
	}

	public void setParms(Double1D p) {
		this.setParms(p, CHECK_PARMS);
	}

	private void checkParms() {
		for (Double0D d : this.p) {
			if (d.value() < 0) {
				throw new IllegalArgumentException("Negative value in p");
			}
		}
	}

	protected void setUpFromParms(boolean checkParms) {
		if (checkParms) {
			this.checkParms();
		}
		// *** FIXME ***
		// *** DEBUGGING ONLY *** currently the posterior computation for
		// P is broken for numerical reasons
		Double1D p = this.p.plus(Double.MIN_VALUE);
		// *** FIXME***
		if (this.empiricalGen == null) {
			this.empiricalGen = new EmpiricalWalker(p.value(),
					Empirical.NO_INTERPOLATION,
					RandomEngineSelector.getEngine());
		} else {
			this.empiricalGen.setState(p.value(), Empirical.NO_INTERPOLATION);
			this.empiricalGen.setState2(p.value());
		}
	}

	/**
	 * @ewan sez fix this
	 */
	@Override
	protected Integer0D genVariate() {
		int[] sample = new int[20];
		for (int i = 0; i < 20; i++) {
			sample[i] = this.empiricalGen.nextInt();
		}
		return new Integer0D(this.empiricalGen.nextInt());
	}

	@Override
	protected double getDensity(Integer0D pt) {
		return this.p.get(pt.value()).value();
	}

}