package org.jqgibbs.models;

import java.util.ArrayList;
import java.util.List;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistInitializeDirectly;
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;

import cern.jet.random.Gamma;
import cern.jet.random.engine.RandomEngine;

public class Deterministic<T extends Numeric<T>> extends ProbDistInitializeDirectly<T> {
	private List<ProbDistParmCheck[]> parmCheck;
	private List<String> parmNames;
	private List<Class<? extends Numeric<?>>> parmClasses;	
	
	public Deterministic(Numeric<?>... parms) throws ProbDistParmException {
		super(parms);
	}
	
	@Override
	protected void installParmChecks() {
		// Names
		this.parmNames = new ArrayList<String>(1);
		this.parmNames.add("answer");
		// Checks
		this.parmCheck = new ArrayList<ProbDistParmCheck[]>(1);
		this.parmCheck.add(null);
		// Classes
		this.parmClasses = new ArrayList<Class<? extends Numeric<?>>>(1);
		this.parmClasses.add(null);
	}
	@SuppressWarnings("unchecked")
	private T getAnswer() {
		return (T) this.parms[0];
	}
	
	@Override
	protected void setUpFromParms() {
		return;
	}	

	@Override
	protected T genVariate() throws ProbDistParmException {
		return this.getAnswer();
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
	protected double getDensity(T pt) {
		if (pt.equals(this.getAnswer())) {
			return 1;
		} else {
			return 0;
		}
	}	
}
