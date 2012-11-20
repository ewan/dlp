package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.List;

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

public class MatrixNormalDist extends ProbDistInitializeDirectly<Double2D> {

	private List<ProbDistParmCheck[]> parmCheck;
	private List<String> parmNames;
	private List<Class<? extends Numeric<?>>> parmClasses;	
	
	private MVNormalDist mvnDist;
	
	public MatrixNormalDist(Numeric<?>... parms) throws ProbDistParmException {
		super(parms);
	}
	
	private MVNormalDist getMVNormalDist() {
		return this.mvnDist;
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
	
	private Double2D getM() {
		return (Double2D) this.parms[0];
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
		return (Double2D) this.parms[1];
	}
	
	private Double2D getOmega() {
		return (Double2D) this.parms[2];
	}
	
	private Double2D getKronSg() {
		return this.getSg().kron(this.getOmega());
	}
	
	@Override
	protected void installParmChecks() {
		// Names
		this.parmNames = new ArrayList<String>(3);
		this.parmNames.add("M");
		this.parmNames.add("Sg");
		this.parmNames.add("Omega");
		// Checks
		this.parmCheck = new ArrayList<ProbDistParmCheck[]>(3);
		this.parmCheck.add(null);
		this.parmCheck.add(new ProbDistParmCheck[] {
			new ProbDistParmCheck() {
				public boolean test(Numeric<?> o) {
					Double2D s = (Double2D) o;
					return (s.square() && s.numCols() == MatrixNormalDist.this.getDims());
				}

				public String message() {
					return "Expected square matrix";
				}
			}
		});
		this.parmCheck.add(new ProbDistParmCheck[] {
			new ProbDistParmCheck() {
				public boolean test(Numeric<?> o) {
					Double2D s = (Double2D) o;
					return (s.square() && s.numCols() == MatrixNormalDist.this.getH());
				}

				public String message() {
					return "Expected square matrix";
				}
			}
		});		
		// Classes
		this.parmClasses = new ArrayList<Class<? extends Numeric<?>>>(3);
		this.parmClasses.add(Double2D.class);
		this.parmClasses.add(Double2D.class);
		this.parmClasses.add(Double2D.class);
	}	

	@Override
	protected void setUpFromParms() throws ProbDistParmException {
		return;
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
