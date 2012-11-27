package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;

import umontreal.iro.lecuyer.probdistmulti.MultiNormalDist;
import umontreal.iro.lecuyer.randvar.NormalGen;
import umontreal.iro.lecuyer.randvarmulti.MultinormalCholeskyGen;
import umontreal.iro.lecuyer.rng.BasicRandomStreamFactory;
import umontreal.iro.lecuyer.rng.MRG32k3a;
import umontreal.iro.lecuyer.rng.RandomStream;
import umontreal.iro.lecuyer.rng.RandomStreamFactory;

public class MVNormalDist extends ProbDist<Double1D> {

	private NormalGen normalGen;
	private MultinormalCholeskyGen mvnGen;
	private MultiNormalDist mvnDist;
	
	private Double1D Mu;
	private Double2D Sg;
	
	private double log2Pi = Math.log(2*Math.PI);
	private double logDetSg;
	private Double2D sgInv;
	
	public MVNormalDist() {
		// Empty (fix?)
	}
	
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

	public MVNormalDist(Double1D Mu, Double2D Sg) throws ProbDistParmException {
		this.Mu = Mu;
		this.Sg = Sg;
		setUpFromParms();
	}

	private Double1D getMu() {
		return Mu;
	}

	private Double2D getSg() {
		return Sg;
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
	
	public void checkParms() throws ProbDistParmException {
		if(! (Sg.square() && Sg.numCols() == getDims())) {
			throw new ProbDistParmException("Expected square matrix for Sg");
		}
	}
	
	private void setUpFromParms() throws ProbDistParmException {
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
