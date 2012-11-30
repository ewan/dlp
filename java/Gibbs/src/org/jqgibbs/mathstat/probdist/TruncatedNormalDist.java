package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.List;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.random.Normal;
import cern.jet.random.Uniform;
import cern.jet.random.engine.RandomEngine;

public class TruncatedNormalDist extends ProbDistInitializeDirectly<Double0D> {

	private final double ALPHA = 2;
	
	private GammaDist gammaDist;
	private Uniform uniformGen;
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

	public TruncatedNormalDist(Numeric<?>... parms) throws ProbDistParmException {
		super(parms);
	}	

	private Double0D getMu() {
		return (Double0D) this.parms[0];
	}	
	
	private Double0D getSg() {
		return (Double0D) this.parms[1];
	}
	
	private Double0D getMin() {
		return (Double0D) this.parms[2];
	}
	
	@Override
	protected void installParmChecks() {
		// Names
		this.parmNames = new ArrayList<String>(3);
		this.parmNames.add("Mu");
		this.parmNames.add("Sg");
		this.parmNames.add("Min");
		// Checks
		this.parmCheck = new ArrayList<ProbDistParmCheck[]>(3);
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
		this.parmCheck.add(null);
		// Classes
		this.parmClasses = new ArrayList<Class<? extends Numeric<?>>>(3);
		this.parmClasses.add(Double0D.class);
		this.parmClasses.add(Double0D.class);
		this.parmClasses.add(Double0D.class);
	}
	
	@Override
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
