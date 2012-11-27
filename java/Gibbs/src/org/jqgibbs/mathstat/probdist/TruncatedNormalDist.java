package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.random.Normal;
import cern.jet.random.Uniform;

public class TruncatedNormalDist extends ProbDist<Double0D> {

	private final double ALPHA = 2;
	
	private GammaDist gammaDist;
	private Uniform uniformGen;
	private Normal normalGen;
	
	private Double0D Mu;
	private Double0D Sg;
	private Double0D Min;
	
	protected void checkInitialized(Numeric... parms) {
		if(parms.length == 0) return;
		this.Mu = (Double0D)parms[0];
		this.Sg = (Double0D)parms[1];
		this.Min = (Double0D)parms[2];
		try {
			setUpFromParms();
		} catch (ProbDistParmException e) {
			e.printStackTrace();
		}
	}

	public TruncatedNormalDist(Double0D Mu, Double0D Sg, Double0D Min, boolean checkParms)
			throws ProbDistParmException {
		this.Mu = Mu;
		this.Sg = Sg;
		this.Min = Min;
		if(checkParms) {
			this.checkParms();
		}
		setUpFromParms();
	}
	
	public TruncatedNormalDist(Double0D Mu, Double0D Sg, Double0D Min)
	throws ProbDistParmException {
		this(Mu, Sg, Min, CHECK_PARMS);
	}
	
	
	private void checkParms() throws ProbDistParmException {
		if(this.Sg.value() < 0) {
			throw new ProbDistParmException("Expected non-negative value for Sg");
		}
	}

	private Double0D getMu() {
		return this.Mu;
	}	
	
	private Double0D getSg() {
		return this.Sg;
	}
	
	private Double0D getMin() {
		return this.Min;
	}
	
	protected void setUpFromParms() throws ProbDistParmException {
		if (this.gammaDist == null) {
			this.gammaDist = new GammaDist(new Double0D(1), new Double0D(this.ALPHA));
		}
		if (this.uniformGen == null) {
			this.uniformGen = new Uniform(0, 1, (int) System.currentTimeMillis());
		}
		if (this.normalGen == null) {
			this.normalGen = new Normal(0, 1, RandomEngineSelector.getEngine());
		}
	}

	@Override
	protected Double0D genVariate() throws ProbDistParmException {
		double min = (this.getMin().value() - this.getMu().value())/this.getSg().value();
		double z;
		if (min < 0) {
			do {
				z = this.normalGen.nextDouble();
			} while (z <= min);
		} else {
			double u, rho;
			do {
				z = this.gammaDist.variateFast().plus(min).value();
				if (min < this.ALPHA) {
					rho = Math.exp(-Math.pow(this.ALPHA-z, 2)/2);
				} else {
					rho = Math.exp((Math.pow(min-this.ALPHA, 2)-Math.pow(this.ALPHA-z, 2))/2);
				}
				u = this.uniformGen.nextDouble();
			} while (u > rho);
		}
		return new Double0D(z*this.getSg().value()+this.getMu().value());
	}

	@Override
	protected double getDensity(Double0D pt) {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}	

}
