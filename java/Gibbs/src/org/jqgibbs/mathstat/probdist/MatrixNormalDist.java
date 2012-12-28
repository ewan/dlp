package org.jqgibbs.mathstat.probdist;

import java.util.List;

import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;

public class MatrixNormalDist extends ProbDist<Double2D> {

	private static MatrixNormalDist staticMNDist = new MatrixNormalDist();

	private MVNormalDist mvnDist = new MVNormalDist();
	private Double2D M;
	private Double2D Sg;
	private Double2D Omega;
	private Double2D SgInv;
	private Double2D OmegaInv;
	private int pOmega;
	private int pSg;
	private Double1D vecM;
	private Double2D kronSg;
	
	private double logNormConst;

	public static Double3D variates(List<Double2D> RepM, Double3D ActiveSg,
			List<Double2D> Omega) {
		boolean singleOmega = Omega.size() < RepM.size(); // using constant
															// Omega
		Double3D sequence = new Double3D();

		for (int i = 0; i < RepM.size(); i++) {
			MatrixNormalDist.staticMNDist.setParms(RepM.get(i),
					ActiveSg.get(i), Omega.get(singleOmega ? 0 : i));
			sequence.add(MatrixNormalDist.staticMNDist.variate());
		}
		return sequence;
	}

	public MatrixNormalDist(Double2D M, Double2D Sg, Double2D Omega,
			boolean checkParms) {
		this.setParms(M, Sg, Omega, checkParms);
	}

	public MatrixNormalDist(Double2D M, Double2D Sg, Double2D Omega) {
		this(M, Sg, Omega, CHECK_PARMS);
	}

	public MatrixNormalDist() {
		super();
	}

	public void setParms(Double2D M, Double2D Sg, Double2D Omega,
			boolean checkParms) {
		this.M = M;
		this.Sg = Sg;
		this.Omega = Omega;
		this.setUpFromParms(checkParms);
		this.initialized = true;
	}

	public void setParms(Double2D M, Double2D Sg, Double2D Omega) {
		this.setParms(M, Sg, Omega, CHECK_PARMS);
	}

	private void checkParms() {
		if (!this.Sg.square()) {
			throw new IllegalArgumentException("Expected square matrix for Sg");
		}
		if (!this.Omega.square()) {
			throw new IllegalArgumentException(
					"Expected square matrix for Omega");
		}
		if (!(this.Sg.numCols() == this.M.numCols())) {
			throw new IllegalArgumentException("M has " + this.M.numCols()
					+ " columns but Sg has " + this.Sg.numCols()
					+ " columns (need to match)");
		}
		if (!(this.Omega.numCols() == this.M.numRows())) {
			throw new IllegalArgumentException("M has " + this.M.numRows()
					+ " rows but Omega has " + this.Omega.numCols()
					+ " columns (need to match)");
		}
	}

	private void setUpFromParms(boolean checkParms) {
		if (checkParms) {
			this.checkParms();
		}
		this.pOmega = this.M.numRows();
		this.pSg = this.M.numCols();
		this.vecM = this.M.colVec();
		this.kronSg = this.Sg.kron(this.Omega);
		this.mvnDist.setParms(this.vecM, this.kronSg);
		this.SgInv = this.Sg.inverse();
		this.OmegaInv = this.Omega.inverse();
		
		this.logNormConst = 0.5*(this.pOmega*this.pSg*NormalDist.log2Pi + this.pOmega*Math.log(this.Sg.det()) + this.pSg*Math.log(this.Omega.det()));
	}

	@Override
	protected Double2D genVariate() {
		return this.mvnDist.variate().toDouble2D(this.pOmega);
	}
	
	@Override
	protected double getLogDensity(Double2D X) {
		Double2D Dev = X.minus(this.M);
		double scaledDist2 = this.SgInv.mult(Dev).mult(this.OmegaInv).mult(Dev).trace().value();
		// FIXME
		if (Double.isInfinite(scaledDist2)) {
			System.err.println("Warning: infinite scaled distance");
			return -Double.MAX_VALUE;
		}
		if (Double.isNaN(scaledDist2)) {
			System.err.println("Warning: NaN scaled distance");
			System.err.println("X: " + X);
			System.err.println("Dev: " + Dev);
			System.err.println("M: " + this.M);
			System.err.println("Sg: " + this.Sg);
			System.err.println("SgInv: " + this.SgInv);
			System.err.println("Omega: " + this.Omega);
			System.err.println("OmegaInv: " + this.OmegaInv);
		}		
		return -0.5*scaledDist2 - this.logNormConst;
	}
}
