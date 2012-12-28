package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;

import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;

/**
 * 
 * Note: creating a NormalDist object advances the pseudo-random number sequence
 * of RandomEngine.getSequence() by one
 * 
 * @author ewan
 * @author josh
 */
public class NormalDist extends ProbDist<Double0D> {

	public static final double log2Pi = Math.log(2*Math.PI);
	
	private RandomEngine randomEngine;
	private Normal normalGen;

	private Double0D mu;
	private Double0D sg;
	
	private double logNormConst;

	public NormalDist() {
		super();
	}

	public NormalDist(Double0D mu, Double0D sg, boolean checkParms) {
		this.setParms(mu, sg, checkParms);
	}

	public NormalDist(Double0D mu, Double0D sg) {
		this(mu, sg, CHECK_PARMS);
	}

	public void setParms(Double0D mu, Double0D sg, boolean checkParms) {
		this.mu = mu;
		this.sg = sg;
		this.setUpFromParms(checkParms);
		this.initialized = true;
	}

	public void setParms(Double0D mu, Double0D sg) {
		this.setParms(mu, sg, CHECK_PARMS);
	}

	private void checkParms() {
		if (this.sg.value() < 0) {
			throw new IllegalArgumentException(
					"Expected non-negative value for sg");
		}
	}

	private void setUpFromParms(boolean checkParms) {
		if (checkParms) {
			this.checkParms();
		}
		if (this.normalGen == null) {
			this.randomEngine = RandomEngineSelector.getEngine();
			this.normalGen = new Normal(this.mu.value(), this.sg.value(),
					this.randomEngine);
		} else {
			this.normalGen.setState(this.mu.value(), this.sg.value());
		}
		this.logNormConst = 0.5*NormalDist.log2Pi + Math.log(this.sg.value());
	}

	@Override
	protected Double0D genVariate() {
		return new Double0D(this.normalGen.nextDouble());
	}

	@Override
	protected double getLogDensity(Double0D x) {
		return -0.5*x.minus(this.mu).divide(this.sg).pow(2).value() - this.logNormConst;
	}
} 