package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;

import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;

public class TruncatedNormalDist extends ProbDist<Double0D> {

	private final double ALPHA = 2;
	
	private GammaDist gammaDist = new GammaDist(new Double0D(1), new Double0D(this.ALPHA));
	private RandomEngine uniformGen = RandomEngineSelector.getEngine();
	private Normal normalGen = new Normal(0, 1, RandomEngineSelector.getEngine());
	
	private Double0D mu;
	private Double0D sg;
	private Double0D min;
	
	public TruncatedNormalDist() {
		super();
	}
	public TruncatedNormalDist(Double0D mu, Double0D sg, Double0D min, boolean checkParms) {
		this.setParms(mu, sg, min, checkParms);
	}
	
	public TruncatedNormalDist(Double0D Mu, Double0D Sg, Double0D Min) {
		this(Mu, Sg, Min, CHECK_PARMS);
	}
	
	public void setParms(Double0D mu, Double0D sg, Double0D min, boolean checkParms) {
		this.mu = mu;
		this.sg = sg;
		this.min = min;
		this.setUpFromParms(checkParms);
		this.initialized = true;
	}
	
	public void setParms(Double0D mu, Double0D sg, Double0D min) {
		this.setParms(mu, sg, min, CHECK_PARMS);
	}
	
	private void checkParms() {
		if(this.sg.value() < 0) {
			throw new IllegalArgumentException("Expected non-negative value for Sg");
		}
	}

	protected void setUpFromParms(boolean checkParms) {
		if (checkParms) {
			this.checkParms();
		}
	}

	@Override
	protected Double0D genVariate() {
		double min = (this.min.value() - this.mu.value())/this.sg.value();
		double z;
		if (min < 0) {
			do {
				z = this.normalGen.nextDouble();
			} while (z <= min);
		} else {
			double u, rho;
			do {
				z = this.gammaDist.variate().plus(min).value();
				if (min < this.ALPHA) {
					rho = Math.exp(-Math.pow(this.ALPHA-z, 2)/2);
				} else {
					rho = Math.exp((Math.pow(min-this.ALPHA, 2)-Math.pow(this.ALPHA-z, 2))/2);
				}
				u = this.uniformGen.nextDouble();
			} while (u > rho);
		}
		return new Double0D(z*this.sg.value()+this.mu.value());
	}

	@Override
	protected double getLogDensity(Double0D pt) {
		// FIXME
		throw new UnsupportedOperationException();
	}	

}
