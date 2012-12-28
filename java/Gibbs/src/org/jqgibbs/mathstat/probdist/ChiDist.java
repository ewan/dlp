package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Integer0D;

import cern.jet.random.ChiSquare;
import cern.jet.stat.Gamma;

public class ChiDist extends ProbDist<Double0D> {

	private Integer0D dof;
	private ChiSquare chisqGen;

	private double logNormConst;

	public ChiDist() {
		super();
	}

	public ChiDist(Integer0D dof, boolean checkParms) {
		this.setParms(dof, checkParms);
	}

	public ChiDist(Integer0D dof) {
		this(dof, CHECK_PARMS);
	}

	public void setParms(Integer0D dof, boolean checkParms) {
		this.dof = dof;
		this.setUpFromParms(checkParms);
		this.initialized = true;
	}

	public void setParms(Integer0D dof) {
		this.setParms(dof, CHECK_PARMS);
	}

	private void checkParms() {
		if (this.dof.value() <= 0) {
			throw new IllegalArgumentException("DOF must be positive");
		}
	}

	private void setUpFromParms(boolean checkParms) {
		if (checkParms) {
			this.checkParms();
		}
		if (this.chisqGen == null) {
			this.chisqGen = new ChiSquare(this.dof.value(),
					RandomEngineSelector.getEngine());
		} else {
			this.chisqGen.setState(this.dof.value());
		}
		this.logNormConst = Gamma.logGamma(this.dof.value() / 2)
				- (1 - this.dof.value() / 2) * Math.log(2);
	}

	@Override
	protected Double0D genVariate() {
		return new Double0D(Math.sqrt(this.chisqGen.nextDouble()));
	}

	@Override
	protected double getLogDensity(Double0D pt) {
		double x = pt.value();
		return (this.dof.value() - 1) * Math.log(x) - Math.pow(x, 2) / 2
				- this.logNormConst;
	}
}
