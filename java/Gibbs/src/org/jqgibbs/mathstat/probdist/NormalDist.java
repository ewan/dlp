package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;

import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;

public class NormalDist extends ProbDist<Double0D> {	

	private RandomEngine randomEngine;
	private Normal normalGen;
	
	private Double0D Mu;
	private Double0D Sg;
	
	public void initializeParms(Double0D Mu, Double0D Sg) throws ProbDistParmException {
		this.Mu = Mu;
		this.Sg = Sg;
		setUpFromParms();
	}

	public NormalDist(Double0D Mu, Double0D Sg, boolean checkParms) throws ProbDistParmException {
		this.Mu = Mu;
		this.Sg = Sg;
		if(checkParms) {
			checkParms();
		}
		setUpFromParms();
	}
	
	public NormalDist(Double0D Mu, Double0D Sg) throws ProbDistParmException {
		this(Mu, Sg, CHECK_PARMS);
	}

	private Double0D getMu() {
		return Mu;
	}	
	
	private Double0D getSg() {
		return Sg;
	}

	private Normal getNormalGen() {
		return normalGen;
	}
	
	private RandomEngine getRandomEngine() {
		return this.randomEngine;
	}
	
	private void checkParms() throws ProbDistParmException {
		if(Sg.value() < 0) {
			throw new ProbDistParmException("Expected non-negative value for Sg");
		}
	}
	
	private void setUpFromParms() throws ProbDistParmException {
		if (this.getNormalGen() == null) {
			this.randomEngine = RandomEngineSelector.getEngine();
			this.normalGen = new Normal(this.getMu().value(), this.getSg().value(),
					this.getRandomEngine());
		} else {
			assert this.getRandomEngine() != null;
			this.getNormalGen().setState(this.getMu().value(), this.getSg().value());
		}
	}

	@Override
	protected Double0D genVariate() throws ProbDistParmException {
		return new Double0D(this.getNormalGen().nextDouble());
	}

	@Override
	protected double getDensity(Double0D pt) {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}
}
