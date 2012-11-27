package org.jqgibbs.mathstat.probdist;

import java.util.List;

import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;

public class MatrixNormalDist extends ProbDist<Double2D> {
	
	private MVNormalDist mvnDist;
	private Double2D M;
	private Double2D Sg;
	private Double2D Omega;
	
	public static Double3D variates(List<Double2D> RepM, Double3D ActiveSg, List<Double2D> Omega) throws ProbDistParmException {
		boolean singleOmega = Omega.size() < RepM.size(); // using constant Omega across entries
		Double3D sequence = new Double3D();
		sequence.add(new MatrixNormalDist(RepM.get(0), ActiveSg.get(0), Omega.get(0)).variate());
		for(int i = 1; i < RepM.size(); i++) {
			sequence.add(new MatrixNormalDist(RepM.get(i), ActiveSg.get(i), Omega.get(singleOmega ? 0 : i), false).variate());
		}
		return sequence;
	}
	
	public MatrixNormalDist(Double2D M, Double2D Sg, Double2D Omega, boolean checkParms) throws ProbDistParmException {
		this.M = M;
		if(this.getM() == null) {
			System.out.println("Why would you do that to me?");
		}
		this.Sg = Sg;
		this.Omega = Omega;
		if(checkParms) {
			checkParms();
		}
	}
	
	public MatrixNormalDist(Double2D M, Double2D Sg, Double2D Omega) throws ProbDistParmException {
		this(M, Sg, Omega, CHECK_PARMS);
	}
	
	public MatrixNormalDist() throws ProbDistParmException {
		//this(new Double2D(), new Double2D(), new Double2D());
		// Empty (fix?)
	}
	
	private MVNormalDist getMVNormalDist() {
		return this.mvnDist;
	}
	
	private Double2D getM() {
		if (this.M == null) {
			System.out.println("how did that happen?");
		}
		return this.M;
	}

	private Double1D getVecM() {
		return this.getM().colVec();
	}
	
	private int getDims() {
		return this.getM().numCols();
	}
		
	private int getH() {
		return this.getM().numRows();
	}
	
	private Double2D getSg() {
		return this.Sg;
	}
	
	private Double2D getOmega() {
		return this.Omega;
	}
	
	private Double2D getKronSg() {
		return this.getSg().kron(this.getOmega());
	}
	
	public void checkParms() throws ProbDistParmException {
		if(!(getSg().square() && getSg().numCols() == MatrixNormalDist.this.getDims())) {
			throw new ProbDistParmException("Expected square matrix for Sg");
		}
		if(!(getOmega().square() && getOmega().numCols() == MatrixNormalDist.this.getH())) {
			throw new ProbDistParmException("Expected square matrix for Omega");
		}
	}
	
	@Override
	protected Double2D genVariate() throws ProbDistParmException {
		Double1D vecVariate;
		if (this.getMVNormalDist() == null) {
			this.mvnDist = new MVNormalDist(this.getVecM(), this.getKronSg());
			vecVariate = this.mvnDist.variate();
		} else {
			vecVariate = this.mvnDist.variate(this.getVecM(), this.getKronSg());
		}
		return vecVariate.toDouble2D(this.getH());
	}
	
	@Override
	protected double getDensity(Double2D pt) throws ProbDistParmException {
		if (this.getMVNormalDist() == null) {
			this.mvnDist = new MVNormalDist(this.getVecM(), this.getKronSg());
		} 
		return this.getMVNormalDist().getDensity(pt.colVec());
	}
	
	@Override
	protected double getLogDensity(Double2D pt) throws ProbDistParmException {
		if (this.getMVNormalDist() == null) {
			this.mvnDist = new MVNormalDist(this.getVecM(), this.getKronSg());
		} 
		return this.getMVNormalDist().getLogDensity(pt.colVec());
	}
}
