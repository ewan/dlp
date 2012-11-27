package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Integer0D;

import cern.jet.random.ChiSquare;

public class ChiDist extends ProbDist<Double0D> {

	private Integer0D DOF;

	private ChiSquare chisqGen;

	public ChiDist(Integer0D DOF, boolean checkParms) throws ProbDistParmException {
		this.DOF = DOF;
		if(checkParms) {
			checkParms();
		}
		setUpFromParms();
	}
	
	public ChiDist(Integer0D DOF) throws ProbDistParmException {
		this(DOF, CHECK_PARMS);
	}
	
	private void checkParms() throws ProbDistParmException {
		if(DOF.value() <= 0) {
			throw new ProbDistParmException("DOF must be positive");
		}
	}

	private Integer0D getDof() {
		return DOF;
	}

	private ChiSquare getChisqGen() {
		return chisqGen;
	}

	private void setUpFromParms() {
		if (this.getChisqGen() == null) {
			this.chisqGen = new ChiSquare(this.getDof().value(), RandomEngineSelector.getEngine());
		} else {
			this.getChisqGen().setState(this.getDof().value());
		}
	}

	@Override
	protected Double0D genVariate() throws ProbDistParmException {
		assert this.initialized;
		return new Double0D(Math.sqrt(this.getChisqGen().nextDouble()));
	}

	@Override
	protected double getDensity(Double0D pt) {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}
}
