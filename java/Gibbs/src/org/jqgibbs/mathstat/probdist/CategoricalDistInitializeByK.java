package org.jqgibbs.mathstat.probdist;

import java.util.Arrays;

import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Numeric;

public class CategoricalDistInitializeByK extends AbstractCategoricalDist {

	public CategoricalDistInitializeByK(Integer0D K, boolean checkParms)
			throws ProbDistParmException {
		double p[] = new double[K.value()];
		Arrays.fill(p, 1 / K.value());
		Double1D P = new Double1D(p);
		this.P = P;
		this.K = K;
		if(checkParms) {
			checkParms();
		}
	}
	
	@Override
	protected void checkInitialized(Numeric... parms) {
		if(parms.length > 0) {
			this.K = (Integer0D)parms[0];
			double p[] = new double[this.K.value()];
			Arrays.fill(p, 1 / this.K.value());
			this.P = new Double1D(p);
			this.setUpEmpiricalGen();
		}
	}
	
	public CategoricalDistInitializeByK(Integer0D K) 
			throws ProbDistParmException {
		this(K, CHECK_PARMS);
	}
	
	private void checkParms() throws ProbDistParmException{
		if(K.value() <= 0) {
			throw new ProbDistParmException("Expected positive value for K");
		}
	}
}
