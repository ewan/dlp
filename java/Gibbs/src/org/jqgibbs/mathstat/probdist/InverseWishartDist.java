package org.jqgibbs.mathstat.probdist;

import java.util.List;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.stat.Gamma;

public class InverseWishartDist extends ProbDist<Double2D> {

	private static InverseWishartDist staticIWDist;

	private WishartDist wishartDist;
	private Double2D Psi;
	private Double0D K;

	protected void checkInitialized(Numeric... parms) {
		if (parms.length < 2) {
			if (!this.initialized) {
				throw new IllegalStateException(
						"use of uninitialized probability distribution");
			}
		} else {
			this.setParms((Double2D) parms[0], (Double0D) parms[1]);
		}
	}

	public static Double3D variates(List<Double2D> postPsi,
			List<Double0D> postKappa) {
		Double3D sequence = new Double3D();
		int i;
		if (InverseWishartDist.staticIWDist == null) {
			InverseWishartDist.staticIWDist = new InverseWishartDist(
					postPsi.get(0), postKappa.get(0));
			sequence.add(InverseWishartDist.staticIWDist.variate());
			i = 1;
		} else {
			i = 0;
		}
		for (; i < postPsi.size(); i++) {
			sequence.add(InverseWishartDist.staticIWDist.variate(
					postPsi.get(i), postKappa.get(i)));
		}
		return sequence;
	}

	public InverseWishartDist(Double2D Psi, Double0D K, boolean checkParms) {
		this.setParms(Psi, K);
	}

	public InverseWishartDist(Double2D Psi, Double0D K) {
		this(Psi, K, CHECK_PARMS);
	}

	public void setParms(Double2D Psi, Double0D K, boolean checkParms) {
		this.Psi = Psi;
		this.K = K;
		this.setUpFromParms(checkParms);
		this.initialized = true;
	}

	public void setParms(Double2D Psi, Double0D K) {
		this.setParms(Psi, K, CHECK_PARMS);
	}

	private void setUpFromParms(boolean checkParms) {
		if (checkParms) {
			this.checkParms();
		}
		if (this.wishartDist == null) {
			this.wishartDist = new WishartDist(this.Psi.inverse(), this.K,
					checkParms);
		} else {
			this.wishartDist.setParms(this.Psi.inverse(), this.K, checkParms);
		}
	}

	private void checkParms() {
		// In principle we ought to check for singularity here
		// Remaining checks will be done inside WishartDist
	}

	@Override
	protected Double2D genVariate() {
		assert this.initialized;
		Double2D W = this.wishartDist.variate();
		return W.inverse();
	}

	@Override
	protected double getDensity(Double2D pt) {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}

	public double getLogDensity(Double2D pt) {
		int d = this.Psi.numCols();
		double k = this.K.value();
		double ld = -k * d * Math.log(2) / 2;
		double detPt = Math.pow(pt.det(), (k + d + 1) / 2);
		double detPsi = Math.pow(this.Psi.det(), k / 2);
		double gamK2 = 0;
		double a = k / 2;
		for (int i = d; i > 1; i--) {
			gamK2 += Math.log(Math.pow(Math.PI, ((double) i - 1) / 2));
			gamK2 += Math.log(Gamma.gamma(a));
			a -= 0.5;
		}
		gamK2 += Math.log(Gamma.gamma(a));
		ld -= Math.log(detPt);
		ld += Math.log(detPsi);
		ld -= gamK2;
		Double2D PsiPtInv = this.Psi.mult(pt.inverse());
		ld -= 0.5 * PsiPtInv.trace().value();
		return ld;
	}
}