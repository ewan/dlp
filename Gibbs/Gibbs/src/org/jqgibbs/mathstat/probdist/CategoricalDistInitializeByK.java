package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Numeric;

public class CategoricalDistInitializeByK extends AbstractCategoricalDist {

	private List<ProbDistParmCheck[]> parmCheck;
	private List<Class<? extends Numeric<?>>> parmClasses;
	private List<String> parmNames;

	public CategoricalDistInitializeByK(Numeric<?>... parms)
			throws ProbDistParmException {
		super(parms);
	}

	@Override
	protected void installParmChecks() {
		// Names
		this.parmNames = new ArrayList<String>(2);
		this.parmNames.add("K");
		this.parmNames.add("P");
		// Checks
		this.parmCheck = new ArrayList<ProbDistParmCheck[]>(2);
		this.parmCheck.add(new ProbDistParmCheck[] { new ProbDistParmCheck() {
			public boolean test(Numeric<?> o) {
				Integer0D d = (Integer0D) o;
				return (d.value() > 0);
			}

			public String message() {
				return "Must be positive";
			}
		} });
		this.parmCheck.add(null);
		// Classes
		this.parmClasses = new ArrayList<Class<? extends Numeric<?>>>(1);
		this.parmClasses.add(Integer0D.class);
		this.parmClasses.add(null);
	}

	@Override
	protected Integer0D getK() {
		return (Integer0D) this.parms[0];
	}

	@Override
	protected Double1D getP() {
		return (Double1D) this.parms[1];
	}

	private void setP(Double1D p) {
		this.parms[1] = p;
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
		double p[] = new double[this.getK().value()];
		Arrays.fill(p, 1 / this.getK().value());
		this.setP(new Double1D(p));
		this.setUpEmpiricalGen();
	}
}
