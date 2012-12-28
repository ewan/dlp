package org.jqgibbs.mathstat.probdist;

import java.util.List;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;

public class InverseWishartDist extends ProbDist<Double2D> {

	private static InverseWishartDist staticIWDist = new InverseWishartDist();

	private WishartDist wishartDist = new WishartDist();
	private Double2D Psi;
	private Double0D k;
	private int p;

	private double logNormConst;

	public InverseWishartDist() {
		super();
	}

	public InverseWishartDist(Double2D Psi, Double0D k, boolean checkParms) {
		this.setParms(Psi, k);
	}

	public InverseWishartDist(Double2D Psi, Double0D k) {
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

	private void checkParms() {
		// In principle we ought to check for singularity here
		// Remaining checks will be done inside WishartDist
	}

	private void setUpFromParms(boolean checkParms) {
		if (checkParms) {
			this.checkParms();
		}
		this.p = this.Psi.numRows();

		this.wishartDist.setParms(this.Psi.inverse(), this.k, checkParms);

		this.logNormConst = (this.k.value() / 2)
				* (p * Math.log(2) - Math.log(this.Psi.det()));
		this.logNormConst += WishartDist.logMvGamma(p, this.k.value() / 2);
	}

	@Override
	protected Double2D genVariate() {
		Double2D W = this.wishartDist.variate();
		return W.inverse();
	}

	public static Double3D variates(List<Double2D> postPsi,
			List<Double0D> postKappa) {
		Double3D sequence = new Double3D();
		for (int i = 0; i < postPsi.size(); i++) {
			InverseWishartDist.staticIWDist.setParms(postPsi.get(i),
					postKappa.get(i));
			sequence.add(InverseWishartDist.staticIWDist.variate());
		}
		return sequence;
	}

	@Override
	protected double getLogDensity(Double2D X) {
		return -0.5
				* (this.Psi.mult(X.inverse()).trace().value() + (this.k.value()
						+ this.p + 1)
						* Math.log(X.det())) - this.logNormConst;
	}
}