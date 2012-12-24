package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.random.ChiSquare;

public class ChiDist extends ProbDist<Double0D> {

	protected void checkInitialized(Numeric... parms) {
		if (parms.length < 1) {
			if (!this.initialized) {
				throw new IllegalStateException(
						"use of uninitialized probability distribution");
			}
		} else {
			this.setParms((Integer0D) parms[0]);
		}
	}

	private Integer0D dof;

	private ChiSquare chisqGen;

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
	}

	@Override
	protected Double0D genVariate() {
		return new Double0D(Math.sqrt(this.chisqGen.nextDouble()));
	}

	@Override
	protected double getDensity(Double0D pt) {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}
}
