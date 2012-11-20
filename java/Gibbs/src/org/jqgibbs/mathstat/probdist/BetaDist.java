package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.List;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.random.Beta;
import cern.jet.random.engine.RandomEngine;

public class BetaDist extends ProbDistInitializeDirectly<Double0D> {
	private List<ProbDistParmCheck[]> parmCheck;
	private List<String> parmNames;
	private List<Class<? extends Numeric<?>>> parmClasses;	
	
	private Beta betaGen;

	public BetaDist(Numeric<?>... parms) throws ProbDistParmException {
		super(parms);
	}
	
	@Override
	protected void installParmChecks() {
		// Names
		this.parmNames = new ArrayList<String>(2);
		this.parmNames.add("shape1");
		this.parmNames.add("shape2");
		// Checks
		this.parmCheck = new ArrayList<ProbDistParmCheck[]>(2);
		ProbDistParmCheck posDouble = new ProbDistParmCheck() {
			public boolean test(Numeric<?> o) {
				Double0D i = (Double0D) o;
				return (i.value() > 0);
			}

			public String message() {
				return "Must be positive";
			}
		};
		this.parmCheck.add(new ProbDistParmCheck[] { posDouble });
		this.parmCheck.add(new ProbDistParmCheck[] { posDouble });
		// Classes
		this.parmClasses = new ArrayList<Class<? extends Numeric<?>>>(2);
		this.parmClasses.add(Double0D.class);
		this.parmClasses.add(Double0D.class);
	}
	
	protected Beta getBetaGen() {
		return this.betaGen;
	}
	
	private Double0D getShape1() {
		return (Double0D) this.parms[0];
	}
	
	private Double0D getShape2() {
		return (Double0D) this.parms[1];
	}		

	@Override
	protected void setUpFromParms() {
		if (this.getBetaGen() == null) {
			this.betaGen = new Beta(this.getShape1().value(), this.getShape2().value(), RandomEngine
					.makeDefault());
		} else {
			this.getBetaGen().setState(this.getShape1().value(), this.getShape2().value());
		}
	}	

	@Override
	protected Double0D genVariate() throws ProbDistParmException {
		assert this.initialized;
		return new Double0D(this.getBetaGen().nextDouble());
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
