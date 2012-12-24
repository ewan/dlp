package org.jqgibbs.mathstat.probdist;

import java.util.List;

import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;
import org.jqgibbs.mathstat.Numeric;

public class MatrixNormalDist extends ProbDist<Double2D> {

	private static MatrixNormalDist staticMNDist;

	private MVNormalDist mvnDist;
	private Double2D M;
	private Double2D Sg;
	private Double2D Omega;
	private int h;
	private Double1D vecM;
	private Double2D kronSg;

	protected void checkInitialized(Numeric... parms) {
		if (parms.length < 3) {
			if (!this.initialized) {
				throw new IllegalStateException(
						"use of uninitialized probability distribution");
			}
		} else {
			this.setParms((Double2D) parms[0], (Double2D) parms[1],
					(Double2D) parms[2]);
		}
	}

	public static Double3D variates(List<Double2D> RepM, Double3D ActiveSg,
			List<Double2D> Omega) {
		boolean singleOmega = Omega.size() < RepM.size(); // using constant
															// Omega across
															// entries
		Double3D sequence = new Double3D();

		int i;
		if (MatrixNormalDist.staticMNDist == null) {
			MatrixNormalDist.staticMNDist = new MatrixNormalDist(RepM.get(0),
					ActiveSg.get(0), Omega.get(0));
			sequence.add(MatrixNormalDist.staticMNDist.variate());
			i = 1;
		} else {
			i = 0;
		}
		for (; i < RepM.size(); i++) {
			sequence.add(MatrixNormalDist.staticMNDist.variate(RepM.get(i),
					ActiveSg.get(i), Omega.get(singleOmega ? 0 : i)));
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

	public MatrixNormalDist() throws ProbDistParmException {
		// this(new Double2D(), new Double2D(), new Double2D());
		// Empty (fix?)
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
		this.h = this.M.numRows();
		this.vecM = this.M.colVec();
		this.kronSg = this.Sg.kron(this.Omega);
	}

	@Override
	protected Double2D genVariate() {
		Double1D vecVariate;
		if (this.mvnDist == null) {
			this.mvnDist = new MVNormalDist(this.vecM, this.kronSg);
			vecVariate = this.mvnDist.variate();
		} else {
			vecVariate = this.mvnDist.variate(this.vecM, this.kronSg);
		}
		return vecVariate.toDouble2D(this.h);
	}

	@Override
	protected double getDensity(Double2D pt) {
		if (this.mvnDist == null) {
			this.mvnDist = new MVNormalDist(this.vecM, this.kronSg);
		}
		return this.mvnDist.getDensity(pt.colVec());
	}

	@Override
	protected double getLogDensity(Double2D pt) {
		if (this.mvnDist == null) {
			this.mvnDist = new MVNormalDist(this.vecM, this.kronSg);
		}
		return this.mvnDist.getLogDensity(pt.colVec());
	}
}
