package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;

import cern.jet.random.Gamma;

public class GammaDist extends ProbDist<Double0D> {
	private Double0D shape;
	private Double0D rate;
	
	private Gamma gammaGen;

	public GammaDist(Double0D shape, Double0D rate, boolean checkParms) throws ProbDistParmException {
		this.shape = shape;
		this.rate = rate;
		if(checkParms) {
			checkParms();
		}
		setUpFromParms();
	}
	
	public GammaDist(Double0D shape, Double0D rate) throws ProbDistParmException {
		this(shape, rate, CHECK_PARMS);
	}
	
	private void checkParms() throws ProbDistParmException {
		if(shape.value() <= 0) {
			throw new ProbDistParmException("shape must be positive");
		}
		if(rate.value() <= 0) {
			throw new ProbDistParmException("rate must be positive");
		}
	}
	
	protected Gamma getGammaGen() {
		return this.gammaGen;
	}
	
	private Double0D getShape() {
		return shape;
	}
	
	private Double0D getRate() {
		return rate;
	}	

	private void setUpFromParms() {
		if (this.getGammaGen() == null) {
			this.gammaGen = new Gamma(this.getShape().value(), this.getRate().value(), RandomEngineSelector.getEngine());
		} else {
			this.getGammaGen().setState(this.getShape().value(), this.getRate().value());
		}
	}

	@Override
	protected double getDensity(Double0D pt) {
		return this.getGammaGen().pdf(pt.value());
	}	
	
	public double getLogDensity(Double0D pt) {
		return Math.log(this.getDensity(pt));
	}
	
	@Override
	protected Double0D genVariate() throws ProbDistParmException {
		assert this.initialized;
		return new Double0D(this.getGammaGen().nextDouble());
	}
}
