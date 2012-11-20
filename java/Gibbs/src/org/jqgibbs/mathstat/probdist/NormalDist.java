package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.List;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;

public class NormalDist extends ProbDistInitializeDirectly<Double0D> {	

	private RandomEngine randomEngine;
	private Normal normalGen;
	
	private List<ProbDistParmCheck[]> parmCheck;
	private List<String> parmNames;
	private List<Class<? extends Numeric<?>>> parmClasses;

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

	public NormalDist(Numeric<?>... parms) throws ProbDistParmException {
		super(parms);
	}	

	private Double0D getMu() {
		return (Double0D) this.parms[0];
	}	
	
	private Double0D getSg() {
		return (Double0D) this.parms[1];
	}

	private Normal getNormalGen() {
		return normalGen;
	}
	
	private RandomEngine getRandomEngine() {
		return this.randomEngine;
	}
	
	@Override
	protected void installParmChecks() {
		// Names
		this.parmNames = new ArrayList<String>(2);
		this.parmNames.add("Mu");
		this.parmNames.add("Sg");
		// Checks
		this.parmCheck = new ArrayList<ProbDistParmCheck[]>(2);
		this.parmCheck.add(null);
		this.parmCheck.add(new ProbDistParmCheck[] { new ProbDistParmCheck() {
			public boolean test(Numeric<?> o) {
				Double0D sg = (Double0D) o;
				return (sg.value() >= 0);
			}
			public String message() {
				return "Expected positive value";
			}
		} });
		// Classes
		this.parmClasses = new ArrayList<Class<? extends Numeric<?>>>(2);
		this.parmClasses.add(Double0D.class);
		this.parmClasses.add(Double0D.class);
	}
	
	@Override
	protected void setUpFromParms() throws ProbDistParmException {
		if (this.getNormalGen() == null) {
			this.randomEngine = RandomEngine.makeDefault();
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
