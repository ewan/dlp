package org.jqgibbs.mathstat.probdist;

import java.util.LinkedList;
import java.util.List;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;

import cern.jet.stat.Gamma;

public class WishartDist extends ProbDist<Double2D> {

	public static final double logPi = Math.log(Math.PI);

	private static final Double0D D0 = new Double0D(0);
	private static final Double0D D1 = new Double0D(1);

	private Double2D Psi;
	private Double0D k;
	private int p;
	
	private double logNormConst;

	private GammaDist[] gammaDists;
	private NormalDist normalDist;

	public static double logMvGamma(int p, double a) {
		double result = 0.25 * p * (p - 1) * WishartDist.logPi;
		result += Gamma.logGamma(a);
		for (int i = 1; i < p; i++) {
			result += Gamma.logGamma(a - 0.5 * i);
		}
		return result;
	}

	public WishartDist() {
		super();
	}

	public WishartDist(Double2D Psi, Double0D k, boolean checkParms) {
		this.setParms(Psi, k, checkParms);
	}

	public WishartDist(Double2D Psi, Double0D k) {
		this(Psi, k, CHECK_PARMS);
	}

	public void setParms(Double2D Psi, Double0D k, boolean checkParms) {
		this.Psi = Psi;
		this.k = k;
		this.setUpFromParms(checkParms);
		this.initialized = true;
	}

	public void setParms(Double2D Psi, Double0D k) {
		this.setParms(Psi, k, CHECK_PARMS);
	}

	public void checkParms() {
		if (!this.Psi.square()) {
			throw new IllegalArgumentException("Expected square matrix for Psi");
		}
		// Positive definiteness will be tested when we try to sample
		if (this.k.value() <= ((double) this.p - 1.0)) {
			String error = "K: Must be >= number of dimensions - 1 ("
					+ String.valueOf(WishartDist.this.p - 1) + ")";
			throw new IllegalArgumentException(error);
		}
	}

	private void setUpFromParms(boolean checkParms) {
		if (checkParms) {
			this.checkParms();
		}
		
		this.p = this.Psi.numRows();
		this.logNormConst = this.k.value() / 2
				* (p * Math.log(2) + Math.log(this.Psi.det()));
		this.logNormConst += WishartDist.logMvGamma(p, this.k.value() / 2);
		
		if (!this.initialized) {
			this.gammaDists = new GammaDist[this.p];
			this.normalDist = new NormalDist(WishartDist.D0, WishartDist.D1);
			for (int i = 0; i < this.p; i++) {
				this.gammaDists[i] = new GammaDist(new Double0D(
						(this.k.value() - i) / 2), new Double0D(0.5));
			}
		} else {
			for (int i = 0; i < this.p; i++) {
				this.gammaDists[i].setParms(new Double0D(
						(this.k.value() - i) / 2), new Double0D(0.5));
			}
		}
	}

	/*
	 * Modified Bartlett decomposition (Ku and Bloomfield). Originally derived
	 * from MCMCpack 1.0-6 for R. Copyright 2009 Andrew D. Martin
	 * <admartin@wustl.edu>, Kevin M. Quinn <kevin_quinn@harvard.edu>, Jong Hee
	 * Park <jhp@uchicago.edu>.
	 */
	protected Double2D genVariate() {
		// Bartlett decomposition: step 1: Cholesky decomposition of scale
		// matrix
		Double2D L;
		try {
			L = this.Psi.cholesky();
		} catch (UnsupportedOperationException e) {
			throw new IllegalArgumentException("Invalid psi parameter", e);
		}
		Double2D A = null;
		boolean aOk = false;
		while (!aOk) {
			// Bartlett decomposition: step 2: construct Z, lower-triangular
			// with the diagonal sqrt of independent gamma with
			// dof v...(v-d+1), rate 0.5, the below-diagonals independent N(0,1)

			// Indices of a d x d matrix that will be sqrt gamma (diagonal),
			// then indices that will be normal (below-diagonals)
			List<int[]> indices = new LinkedList<int[]>();
			for (int i = 0; i < this.p; i++) {
				indices.add(new int[] { i, i });
			} // Sqrt gamma
			for (int i = 0; i < (this.p - 1); i++) {
				for (int j = (i + 1); j < this.p; j++) {
					indices.add(new int[] { i, j });
				} // Normal
			}
			// Sqrt gamma variates, then normal variates
			double[] v = new double[this.p + (this.p - 1) * this.p / 2];
			for (int i = 0; i < this.p; i++) {
				v[i] = this.gammaDists[i].variate().sqrt().value();
			}
			for (int i = 0; i < (this.p - 1) * this.p / 2; i++) {
				v[p + i] = this.normalDist.variate().value();
			}
			Double1D values = new Double1D(v);
			Double2D Z = new Double2D(this.p, this.p, indices, values);
			// Bartlett decomposition: step 3: return LZZ'L'
			Double2D LZ = L.mult(Z);
			A = LZ.mult(LZ.transpose());
			aOk = A.isWellConditioned();
			if (!aOk) {
				System.err
						.println("Warning: threw out ill-conditioned Wishart sample");
			}
		}
		return A;
	}

	@Override
	protected double getLogDensity(Double2D X) {
		return 0.5
				* (-this.Psi.inverse().mult(X).trace().value() + (this.k
						.value() - this.p - 1)
						* Math.log(X.det())) - this.logNormConst;
	}
}
