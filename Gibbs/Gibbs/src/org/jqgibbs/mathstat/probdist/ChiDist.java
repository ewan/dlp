package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.List;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.random.ChiSquare;
import cern.jet.random.engine.RandomEngine;

public class ChiDist extends ProbDistInitializeDirectly<Double0D> {

	private List<ProbDistParmCheck[]> parmCheck;
	private List<String> parmNames;
	private List<Class<? extends Numeric<?>>> parmClasses;

	private ChiSquare chisqGen;

	public ChiDist(Numeric<?>... parms) throws ProbDistParmException {
		super(parms);
	}

	@Override
	protected void installParmChecks() {
		// Names
		this.parmNames = new ArrayList<String>(1);
		this.parmNames.add("DOF");
		// Checks
		this.parmCheck = new ArrayList<ProbDistParmCheck[]>(1);
		this.parmCheck.add(new ProbDistParmCheck[] { new ProbDistParmCheck() {
			public boolean test(Numeric<?> o) {
				Integer0D i = (Integer0D) o;
				return (i.value() > 0);
			}

			public String message() {
				return "Must be positive";
			}
		} });
		// Classes
		this.parmClasses = new ArrayList<Class<? extends Numeric<?>>>(1);
		this.parmClasses.add(Integer0D.class);
	}

	private Integer0D getDof() {
		return (Integer0D) this.parms[0];
	}

	private ChiSquare getChisqGen() {
		return chisqGen;
	}

	@Override
	protected void setUpFromParms() {
		if (this.getChisqGen() == null) {
			this.chisqGen = new ChiSquare(this.getDof().value(), RandomEngine
					.makeDefault());
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
	protected List<ProbDistParmCheck[]> getParmCheck() {
		return this.parmCheck;
	}

	@Override
	protected List<Class<? extends Numeric<?>>> getParmClasses() {
		return this.parmClasses;
	}

	@Override
	protected List<String> getParmNames() {
		return this.parmNames;
	}

	@Override
	protected double getDensity(Double0D pt) {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}
}
