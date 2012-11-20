package org.jqgibbs.models;

import java.util.ArrayList;
import java.util.Map;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Integer1D;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.RandomVar;
import org.jqgibbs.mathstat.probdist.CategoricalDistInitializeByP;
import org.jqgibbs.mathstat.probdist.MVNormalDist;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistInitializeByChain;
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;

public class MpgGModel extends MpgFModel {

	public MpgGModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
			Map<String, Numeric<? extends Numeric<?>>> init, int dims)
			throws ProbDistParmException {
		super(hypers, init, dims);
		
		
		ProbDist<Double0D> priorNuB = new Deterministic<Double0D>(init.get("nub"));
		ProbDistInitializeByChain<Double0D> postNuB = new ProbDistInitializeByChain<Double0D>() {
			private Double0D getNuB() {
				return (Double0D) this.chainParms[0];
			}
			
			@Override
			protected void installParmChecks() {
				// Names
				this.chainParmNames = new ArrayList<String>(1);
				this.chainParmNames.add("nub");
				this.fixedParmNames = new ArrayList<String>(0);
				// Checks
				this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
				this.chainParmCheck.add(null); // FIXME
				this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
				// Classes
				this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(1);
				this.chainParmClasses.add(Double0D.class);
				this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(0);
			}

			@Override
			protected void setUpFromChainParms() {
				return;
			}

			@Override
			protected Double0D genVariate() throws ProbDistParmException {
				return this.getNuB();
			}

			@Override
			protected double getDensity(Double0D pt) {
				throw new UnsupportedOperationException(
						"Too lazy, come back later");
			}
		};
		Double0D defaultNuB = (Double0D) init.get("nub"); // FIXME
		RandomVar<Double0D> rvNuB = new RandomVar<Double0D>("nub", priorNuB, postNuB,
				defaultNuB);
		this.params.put("nub", rvNuB);		
	}

}
