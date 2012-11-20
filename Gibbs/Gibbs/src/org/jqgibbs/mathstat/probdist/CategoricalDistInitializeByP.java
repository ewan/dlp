package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.List;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Numeric;

public class CategoricalDistInitializeByP extends
		AbstractCategoricalDist {

	private List<ProbDistParmCheck[]> parmCheck;
	private List<Class<? extends Numeric<?>>> parmClasses;
	private List<String> parmNames;
	
	public CategoricalDistInitializeByP(Numeric<?>... parms)
			throws ProbDistParmException {
		super(parms);
	}
	
	@Override
	protected void installParmChecks() {
		// Names
		this.parmNames = new ArrayList<String>(2);
		this.parmNames.add("P");
		this.parmNames.add("K");
		// Checks
		this.parmCheck = new ArrayList<ProbDistParmCheck[]>(2);
		this.parmCheck.add(
				new ProbDistParmCheck[] {
					new ProbDistParmCheck() {
						public boolean test(Numeric<?> o) {
							Double1D ds = (Double1D) o;
							for (Double0D d : ds) {
								if (d.value() < 0) {
									return false;
								}
							}
							return true;
						}
						public String message() {
							return "Must be nonnegative";
						}
					}
				}
		);
		this.parmCheck.add(null);
		// Classes
		this.parmClasses = new ArrayList<Class<? extends Numeric<?>>>(1);
		this.parmClasses.add(Double1D.class);
		this.parmClasses.add(null);
	}

	@Override
	protected Double1D getP() {
		return (Double1D) this.parms[0];
	}
	
	@Override
	protected Integer0D getK() {
		return (Integer0D) this.parms[1];
	}
	
	protected void setK(Integer0D K) {
		this.parms[1] = K;
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
	protected void setUpFromParms() {
		this.setK(new Integer0D(this.getP().size()));
		this.setUpEmpiricalGen();
	}

}
