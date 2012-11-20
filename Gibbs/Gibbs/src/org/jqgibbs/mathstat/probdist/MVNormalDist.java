package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.List;

import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Numeric;

import cern.colt.matrix.DoubleMatrix2D;

import umontreal.iro.lecuyer.probdistmulti.MultiNormalDist;
import umontreal.iro.lecuyer.randvar.NormalGen;
import umontreal.iro.lecuyer.randvarmulti.MultinormalCholeskyGen;
import umontreal.iro.lecuyer.rng.BasicRandomStreamFactory;
import umontreal.iro.lecuyer.rng.MRG32k3a;
import umontreal.iro.lecuyer.rng.RandomStream;
import umontreal.iro.lecuyer.rng.RandomStreamFactory;

public class MVNormalDist extends ProbDistInitializeDirectly<Double1D> {

	private NormalGen normalGen;
	private MultinormalCholeskyGen mvnGen;
	private MultiNormalDist mvnDist;
	
	private List<ProbDistParmCheck[]> parmCheck;
	private List<String> parmNames;
	private List<Class<? extends Numeric<?>>> parmClasses;
	
	private double log2Pi = Math.log(2*Math.PI);
	private double logDetSg;
	private Double2D sgInv;
	
	protected double getLog2Pi() {
		return this.log2Pi;
	}
	
	protected void setLogDetSg(Double2D sg) {
		this.logDetSg = Math.log(sg.det());
	}
	
	protected double getLogDetSg() {
		return this.logDetSg;
	}
	
	protected void setSgInv(Double2D sg) {
		this.sgInv = sg.inverse();
	}
	
	protected Double2D getSgInv() {
		return this.sgInv;
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

	public MVNormalDist(Numeric<?>... parms) throws ProbDistParmException {
		super(parms);
	}

	private Double1D getMu() {
		return (Double1D) this.parms[0];
	}

	private Double2D getSg() {
		return (Double2D) this.parms[1];
	}	
	
	private int getDims() {
		return this.getMu().size();
	}
	
	private NormalGen getNormalGen() {
		return this.normalGen;
	}

	private MultinormalCholeskyGen getMVNormalGen() {
		return this.mvnGen;
	}	

	private MultiNormalDist getMVNormalDist() {
		return this.mvnDist;
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
				Double2D sg = (Double2D) o;
				return (sg.square() && sg.numCols() == MVNormalDist.this
						.getDims());
			}

			public String message() {
				return "Expected square matrix";
			}
		} });
		// Classes
		this.parmClasses = new ArrayList<Class<? extends Numeric<?>>>(2);
		this.parmClasses.add(Double1D.class);
		this.parmClasses.add(Double2D.class);
	}
	
	@Override
	protected void setUpFromParms() throws ProbDistParmException {
		if (this.getMVNormalGen() == null) {
			Class<MRG32k3a> c = MRG32k3a.class;
			RandomStreamFactory rsf = new BasicRandomStreamFactory(c);
			RandomStream rs = rsf.newInstance();
			this.normalGen = new NormalGen(rs);
			this.mvnGen = new MultinormalCholeskyGen(this.getNormalGen(), this
					.getMu().value(), this.getSg().toColt());
		} else {
			assert this.getNormalGen() != null;
			this.getMVNormalGen().setMu(this.getMu().value());
			this.getMVNormalGen().setSigma(this.getSg().toColt());
		}
		if (this.getMVNormalDist() == null) {
			this.mvnDist = new MultiNormalDist(this.getMu().value(), this.getSg().value());
		} else {
			assert this.getMVNormalDist() != null;
			this.getMVNormalDist().setParams(this.getMu().value(), this.getSg().value());
		}
		this.setLogDetSg(this.getSg());
		this.setSgInv(this.getSg());
	}

	@Override
	protected Double1D genVariate() throws ProbDistParmException {
		double[] p = new double[this.getDims()];
		this.getMVNormalGen().nextPoint(p);
		return new Double1D(p);
	}

	@Override
	protected double getDensity(Double1D pt) throws ProbDistParmException {
		return this.getMVNormalDist().density(pt.value());
	}
	
	@Override
	protected double getLogDensity(Double1D pt) throws ProbDistParmException {
		Double1D dev = pt.minus(this.getMu());
		double mahal = dev.mult(this.getSgInv()).mult(dev);
		if (Double.isInfinite(mahal)) {
			System.err.println("Warning: infinite Mahalanobis distance");
			return -Double.MAX_VALUE;
		}
		if (Double.isNaN(mahal)) {
			System.err.println("Warning: NaN Mahalanobis distance");
			System.err.println("pt: " + pt);
			System.err.println("dev: " + dev);
			System.err.println("mu: " + this.getMu());
			System.err.println("sg: " + this.getSg());
			System.err.println("sgInv: " + this.getSgInv());
		}
		return -0.5*(this.getDims()*this.getLog2Pi() + this.getLogDetSg() + mahal);
	}
}
