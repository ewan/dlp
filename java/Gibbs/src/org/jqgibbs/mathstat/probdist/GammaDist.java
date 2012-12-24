package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.random.Gamma;

public class GammaDist extends ProbDist<Double0D> {
	private Double0D shape;
	private Double0D rate;

	private Gamma gammaGen;

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

	public GammaDist(Double0D shape, Double0D rate, boolean checkParms) {
		this.setParms(shape, rate, checkParms);
	}

	public GammaDist(Double0D shape, Double0D rate) {
		this(shape, rate, CHECK_PARMS);
	}

	public void setParms(Double0D shape, Double0D rate, boolean checkParms) {
		this.shape = shape;
		this.rate = rate;
		this.setUpFromParms(checkParms);
		this.initialized = true;
	}

	public void setParms(Double0D shape, Double0D rate) {
		this.setParms(shape, rate, CHECK_PARMS);
	}

	private void checkParms() {
		if (this.shape.value() <= 0) {
			throw new IllegalArgumentException("shape must be positive");
		}
		if (this.rate.value() <= 0) {
			throw new IllegalArgumentException("rate must be positive");
		}
	}

	private void setUpFromParms(boolean checkParms) {
		if (checkParms) {
			this.checkParms();
		}
		if (this.gammaGen == null) {
			this.gammaGen = new Gamma(this.shape.value(), this.rate.value(),
					RandomEngineSelector.getEngine());
		} else {
			this.gammaGen.setState(this.shape.value(), this.rate.value());
		}
	}

	@Override
	protected double getDensity(Double0D pt) {
		return this.gammaGen.pdf(pt.value());
	}

	public double getLogDensity(Double0D pt) {
		return Math.log(this.getDensity(pt));
	}

	@Override
	protected Double0D genVariate() {
		return new Double0D(this.gammaGen.nextDouble());
	}
}
