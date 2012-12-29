package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;

import umontreal.iro.lecuyer.randvar.NormalGen;
import umontreal.iro.lecuyer.randvarmulti.MultinormalCholeskyGen;
import umontreal.iro.lecuyer.rng.RandomStream;

public class MVNormalDist extends ProbDist<Double1D> {

	private NormalGen normalGen;
	private MultinormalCholeskyGen mvnGen;

	private Double1D mu;
	private Double2D Sg;

	private Double2D SgInv;

	private int p;
	private double logNormConst;
	
	public MVNormalDist() {
		super();
	}

	public MVNormalDist(Double1D mu, Double2D Sg, boolean checkParms) {
		this.setParms(mu, Sg, CHECK_PARMS);
	}

	public MVNormalDist(Double1D mu, Double2D Sg) {
		this(mu, Sg, CHECK_PARMS);
	}

	public void setParms(Double1D mu, Double2D Sg, boolean checkParms) {
		this.mu = mu;
		this.Sg = Sg;
		this.setUpFromParms(checkParms);
		this.initialized = true;
	}

	public void setParms(Double1D mu, Double2D Sg) {
		this.setParms(mu, Sg, CHECK_PARMS);
	}

	public void checkParms() {
		if (!this.Sg.square()) {
			throw new IllegalArgumentException("Expected square matrix for Sg");
		}
		if (!(this.Sg.numRows() == this.mu.size())) {
			throw new IllegalArgumentException("mu is of length "
					+ this.mu.size() + ", but Sg is " + this.Sg.numRows()
					+ " by " + this.Sg.numRows());
		}
	}

	private void setUpFromParms(boolean checkParms) {
		if (checkParms) {
			this.checkParms();
		}
		this.p = this.mu.size();
		if (!this.initialized) {
			RandomStream rs = RandomEngineSelector.getStream();
			this.normalGen = new NormalGen(rs);
			this.mvnGen = new MultinormalCholeskyGen(this.normalGen,
					this.mu.value(), this.Sg.toColt());
		} else {
			this.mvnGen.setMu(this.mu.value());
			this.mvnGen.setSigma(this.Sg.toColt());
		}
		this.SgInv = this.Sg.inverse();
		
		this.logNormConst = 0.5*(Math.log(this.Sg.det()) + this.p*NormalDist.log2Pi);
	}

	@Override
	protected Double1D genVariate() {
		double[] pt = new double[this.p];
		this.mvnGen.nextPoint(pt);
		return new Double1D(pt);
	}
	
	@Override
	protected double getLogDensity(Double1D x) {
		Double1D dev = x.minus(this.mu);
		double scaledDist2 = dev.mult(this.SgInv).mult(dev);
		// FIXME
		if (Double.isInfinite(scaledDist2)) {
			System.err.println("Warning: infinite Mahalanobis distance");
			return -Double.MAX_VALUE;
		}
		if (Double.isNaN(scaledDist2)) {
			System.err.println("Warning: NaN Mahalanobis distance");
			System.err.println("x: " + x);
			System.err.println("dev: " + dev);
			System.err.println("mu: " + this.mu);
			System.err.println("sg: " + this.Sg);
			System.err.println("sgInv: " + this.SgInv);
		}
		return -0.5*scaledDist2 - this.logNormConst;
	}
}
