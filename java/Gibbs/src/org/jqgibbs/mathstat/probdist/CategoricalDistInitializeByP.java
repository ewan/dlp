package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Numeric;

public class CategoricalDistInitializeByP extends
		AbstractCategoricalDist {
	
	public CategoricalDistInitializeByP(Double1D P, boolean checkParms)
			throws ProbDistParmException {
		Integer0D K = new Integer0D(P.size());
		this.K = K;
		this.P = P;
		if(checkParms) {
			checkParms();
		}
		this.setUpEmpiricalGen();
	}
	
	@Override
	protected void checkInitialized(Numeric... parms) {
		if(parms.length > 0) {
			this.P = (Double1D)parms[0];
			this.K = new Integer0D(this.P.size());
			this.setUpEmpiricalGen();
		}
	}
	
	private void checkParms() throws ProbDistParmException {
		for (Double0D d : this.getP()) {
			if(d.value() < 0) {
				throw new ProbDistParmException("Negative value in P");
			}
		}
	}
	
	public CategoricalDistInitializeByP(Double1D P)
			throws ProbDistParmException {
		this(P, CHECK_PARMS);
	}
}