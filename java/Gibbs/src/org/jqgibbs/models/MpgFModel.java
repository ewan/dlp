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
import org.jqgibbs.mathstat.probdist.CategoricalDist;
import org.jqgibbs.mathstat.probdist.MVNormalDist;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistInitializeByChain;
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;

public class MpgFModel extends MpgModel {

	public MpgFModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
			Map<String, Numeric<? extends Numeric<?>>> init, int dims)
			throws ProbDistParmException {
		super(hypers, init, dims);
		
		ProbDist<Integer1D> priorB = new Deterministic<Integer1D>(init.get("B"));
		ProbDistInitializeByChain<Integer1D> postB = new ProbDistInitializeByChain<Integer1D>() {
			private Integer1D getB() {
				return (Integer1D) this.chainParms[0];
			}
			
			@Override
			protected void installParmChecks() {
				// Names
				this.chainParmNames = new ArrayList<String>(1);
				this.chainParmNames.add("B");
				this.fixedParmNames = new ArrayList<String>(0);
				// Checks
				this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
				this.chainParmCheck.add(null); // FIXME
				this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
				// Classes
				this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(1);
				this.chainParmClasses.add(Integer1D.class);
				this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(0);
			}

			@Override
			protected void setUpFromChainParms() {
				return;
			}

			@Override
			protected Integer1D genVariate() throws ProbDistParmException {
				return this.getB();
			}

			@Override
			protected double getDensity(Integer1D pt) {
				throw new UnsupportedOperationException(
						"Too lazy, come back later");
			}
		};
		Integer1D defaultB = (Integer1D) init.get("B"); // FIXME
		RandomVar<Integer1D> rvB = new RandomVar<Integer1D>("B", priorB, postB,
				defaultB);
		this.params.put("B", rvB);
		
	}

}
