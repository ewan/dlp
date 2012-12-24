package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Numeric;

import umontreal.iro.lecuyer.probdistmulti.MultiNormalDist;
import umontreal.iro.lecuyer.randvar.NormalGen;
import umontreal.iro.lecuyer.randvarmulti.MultinormalCholeskyGen;
import umontreal.iro.lecuyer.rng.BasicRandomStreamFactory;
import umontreal.iro.lecuyer.rng.MRG32k3a;
import umontreal.iro.lecuyer.rng.RandomStream;
import umontreal.iro.lecuyer.rng.RandomStreamFactory;

public class MVNormalDist extends ProbDist<Double1D> {

	private static double log2Pi = Math.log(2 * Math.PI);

	private NormalGen normalGen;
	private MultinormalCholeskyGen mvnGen;
	private MultiNormalDist mvnDist;

	private Double1D mu;
	private Double2D Sg;

	private double logDetSg;
	private Double2D SgInv;

	private int p;

	protected void checkInitialized(Numeric... parms) {
		if (parms.length < 2) {
			if (!this.initialized) {
				throw new IllegalStateException(
						"use of uninitialized probability distribution");
			}
		} else {
			this.setParms((Double1D) parms[0], (Double2D) parms[1]);
		}
	}

	public MVNormalDist() {
		// Empty (fix?)
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
			this.mvnDist = new MultiNormalDist(this.mu.value(), this.Sg.value());
		} else {
			this.mvnGen.setMu(this.mu.value());
			this.mvnGen.setSigma(this.Sg.toColt());
			this.mvnDist.setParams(this.mu.value(), this.Sg.value());
		}
		this.logDetSg = Math.log(this.Sg.det());
		this.SgInv = this.Sg.inverse();
	}

	@Override
	protected Double1D genVariate() {
		double[] pt = new double[this.p];
		this.mvnGen.nextPoint(pt);
		return new Double1D(pt);
	}

	@Override
	protected double getDensity(Double1D pt) {
		return this.mvnDist.density(pt.value());
	}

	@Override
	protected double getLogDensity(Double1D pt) {
		Double1D dev = pt.minus(this.mu);
		double mahal = dev.mult(this.SgInv).mult(dev);
		if (Double.isInfinite(mahal)) {
			System.err.println("Warning: infinite Mahalanobis distance");
			return -Double.MAX_VALUE;
		}
		if (Double.isNaN(mahal)) {
			System.err.println("Warning: NaN Mahalanobis distance");
			System.err.println("pt: " + pt);
			System.err.println("dev: " + dev);
			System.err.println("mu: " + this.mu);
			System.err.println("sg: " + this.Sg);
			System.err.println("sgInv: " + this.SgInv);
		}
		return -0.5 * (this.p * MVNormalDist.log2Pi + this.logDetSg + mahal);
	}
}
