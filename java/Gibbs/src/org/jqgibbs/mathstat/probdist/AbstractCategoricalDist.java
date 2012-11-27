package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.random.Empirical;
import cern.jet.random.EmpiricalWalker;

public abstract class AbstractCategoricalDist extends
		ProbDist<Integer0D> {

	private EmpiricalWalker empiricalGen;
	protected Integer0D K;
	protected Double1D P;

	protected void setEmpiricalGen(EmpiricalWalker empiricalGen) {
		this.empiricalGen = empiricalGen;
	}

	protected EmpiricalWalker getEmpiricalGen() {
		return empiricalGen;
	}

	protected Integer0D getK() {
		return this.K;
	}

	protected Double1D getP() {
		return this.P;
	}

	protected void setUpEmpiricalGen() {
		// *** FIXME ***
		// *** DEBUGGING ONLY *** currently the posterior computation for
		// P is broken for numerical reasons
		Double1D p = this.getP();
		p = p.plus(Double.MIN_VALUE);
		// *** FIXME***
		if (this.getEmpiricalGen() == null) {
			this.setEmpiricalGen(new EmpiricalWalker(p.value(),
					Empirical.NO_INTERPOLATION, RandomEngineSelector.getEngine()));
		} else {
			assert this.initialized == true;
			this.getEmpiricalGen().setState(p.value(),
					Empirical.NO_INTERPOLATION);
			this.getEmpiricalGen().setState2(p.value());
		}
	}

	@Override
	protected Integer0D genVariate() {
		assert this.initialized == true;
		int[] sample = new int[20];
		for (int i=0; i<20; i++) {
			sample[i] = this.getEmpiricalGen().nextInt();
		}
		return new Integer0D(this.getEmpiricalGen().nextInt());
	}

	@Override
	protected double getDensity(Integer0D pt) {
		return this.getP().get(pt.value()).value();
	}
}
