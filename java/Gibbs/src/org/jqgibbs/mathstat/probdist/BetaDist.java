package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.random.Beta;

public class BetaDist extends ProbDist<Double0D> {
	private Double0D shape1;
	private Double0D shape2;

	private Beta betaGen;

	protected void checkInitialized(Numeric... parms) {
		if (parms.length < 2) {
			if (!this.initialized) {
				throw new IllegalStateException(
						"use of uninitialized probability distribution");
			}
		} else {
			this.setParms((Double0D) parms[0], (Double0D) parms[1]);
		}
	}

	public BetaDist(Double0D shape1, Double0D shape2, boolean checkParms) {
		this.setParms(shape1, shape2, checkParms);
	}

	public BetaDist(Double0D shape1, Double0D shape2) {
		this(shape1, shape2, CHECK_PARMS);
	}

	public void setParms(Double0D shape1, Double0D shape2, boolean checkParms) {
		this.shape1 = shape1;
		this.shape2 = shape2;
		this.setUpFromParms(checkParms);
		this.initialized = true;
	}

	public void setParms(Double0D shape1, Double0D shape2) {
		this.setParms(shape1, shape2, CHECK_PARMS);
	}

	private void checkParms() {
		if (this.shape1.value() <= 0) {
			throw new IllegalArgumentException("shape1 must be positive");
		}
		if (this.shape2.value() <= 0) {
			throw new IllegalArgumentException("shape2 must be positive");
		}
	}

	private void setUpFromParms(boolean checkParms) {
		if (checkParms) {
			this.checkParms();
		}
		if (this.betaGen == null) {
			this.betaGen = new Beta(this.shape1.value(), this.shape2.value(),
					RandomEngineSelector.getEngine());
		} else {
			this.betaGen.setState(this.shape1.value(), this.shape2.value());
		}
	}

	@Override
	protected Double0D genVariate() {
		return new Double0D(this.betaGen.nextDouble());
	}

	@Override
	protected double getDensity(Double0D pt) {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}
}
