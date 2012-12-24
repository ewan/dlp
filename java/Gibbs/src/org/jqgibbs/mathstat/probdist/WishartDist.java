package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.AbstractSequence;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Numeric;

import cern.colt.matrix.linalg.LUDecomposition;
import cern.jet.random.engine.RandomEngine;

public class WishartDist extends ProbDistInitializeDirectly<Double2D> {

	private static Double0D D0 = new Double0D(0);
	private static Double0D D1 = new Double0D(1);

	private List<String> parmNames;
	private List<ProbDistParmCheck[]> parmCheck;
	private List<Class<? extends Numeric<?>>> parmClasses;

	private GammaDist[] gammaDists;
	private NormalDist normalDist;

	public WishartDist(Numeric<?>... parms) throws ProbDistParmException {
		super(parms);
	}

	@Override
	protected void installParmChecks() {
		// Names
		this.parmNames = new ArrayList<String>(2);
		this.parmNames.add("Psi");
		this.parmNames.add("K");
		// Checks
		this.parmCheck = new ArrayList<ProbDistParmCheck[]>(2);
		this.parmCheck.add(new ProbDistParmCheck[] { new ProbDistParmCheck() {
			public boolean test(Numeric<?> o) {
				Double2D psi = (Double2D) o;
				return (psi.square());
			}

			public String message() {
				return "Expected square matrix";
			}
		} });
		this.parmCheck.add(new ProbDistParmCheck[] { new ProbDistParmCheck() {
			public boolean test(Numeric<?> o) {
				Double0D K = (Double0D) o;
				return (K.value() > ((double) WishartDist.this.getDims() - 1.0));
			}

			public String message() {
				return "Must be at > number of dimensions - 1 ("
						+ String.valueOf(WishartDist.this.getDims() - 1) + ")";
			}
		} });
		// Classes
		this.parmClasses = new ArrayList<Class<? extends Numeric<?>>>(2);
		this.parmClasses.add(Double2D.class);
		this.parmClasses.add(Double0D.class);
	}

	@Override
	protected void setUpFromParms() {
		try {
			if (this.getGammaDists() == null) {
				this.gammaDists = new GammaDist[this.getDims()];
				for (int i=0; i<this.getDims(); i++) {
					this.gammaDists[i] = new GammaDist(new Double0D((this.getK().value()-i)/2), new Double0D(0.5));
				}
			} else {
				for (int i=0; i<this.getDims(); i++) {
					this.gammaDists[i] = new GammaDist(new Double0D((this.getK().value()-i)/2), new Double0D(0.5));
				}				
			}
			
			if (this.getNormalDist() == null) {
				this.normalDist = new NormalDist(WishartDist.D0, WishartDist.D1);
			} else {
				this.getNormalDist().initializeParms(WishartDist.D0, WishartDist.D1);
			}
		} catch (ProbDistParmException e) {
			throw new RuntimeException("Error initializing Wishart generator", e);
		}
	}

	private Double2D getPsi() {
		return (Double2D) this.parms[0];
	}

	private Double0D getK() {
		return (Double0D) this.parms[1];
	}

	private int getDims() {
		return this.getPsi().numCols();
	}

	private GammaDist[] getGammaDists() {
		return gammaDists;
	}

	private NormalDist getNormalDist() {
		return normalDist;
	}

	/*
	 * Modified Bartlett decomposition (Ku and Bloomfield).
	 * Originally derived from MCMCpack 1.0-6 for R. Copyright 2009
	 * Andrew D. Martin <admartin@wustl.edu>, Kevin M. Quinn
	 * <kevin_quinn@harvard.edu>, Jong Hee Park <jhp@uchicago.edu>.
	 */
	protected Double2D genVariate() throws ProbDistParmException {
		// Bartlett decomposition: step 1: Cholesky decomposition of scale
		// matrix
		Double2D L;
		try {
			L = this.getPsi().cholesky();
		} catch (UnsupportedOperationException e) {
			throw new ProbDistParmException("Invalid psi parameter", e);
		}
		Double2D A = null;
		boolean aOk = false;
		while (!aOk) {
			// Bartlett decomposition: step 2: construct Z, lower-triangular with
			// the diagonal sqrt of independent gamma
			// with dof v...(v-d+1), rate 0.5, the below-diagonals independent N(0,1)

			// Indices of a d x d matrix that will be sqrt gamma (diagonal), then indices
			// that will be normal
			// (below-diagonals)
			int p = this.getDims();
			List<int[]> indices = new LinkedList<int[]>();
			for (int i = 0; i < p; i++) {
				indices.add(new int[] { i, i });
			} // Sqrt gamma
			for (int i = 0; i < (p-1); i++) {
				for (int j = (i + 1); j < p; j++) {
					indices.add(new int[] { i, j });
				} // Normal
			}
			// Sqrt gamma variates, then normal variates
			double[] v = new double[p + (p-1)*p/2];
			for (int i = 0; i < p; i++) {
				v[i] = this.getGammaDists()[i].variate().sqrt().value();
			}
			for (int i = 0; i < (p-1)*p/2; i++) {
				v[p+i] = this.getNormalDist().variate().value();
			}
			Double1D values = new Double1D(v);
			Double2D Z = new Double2D(p, p, indices, values);
			// Bartlett decomposition: step 3: return LZZ'L'
			Double2D LZ = L.mult(Z);
			A = LZ.mult(LZ.transpose());
			aOk = A.isWellConditioned();
			if (!aOk) {
				System.err.println("Warning: threw out ill-conditioned Wishart sample");
			}
		}
		return A;
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

	@Override
	protected double getDensity(Double2D pt) throws ProbDistParmException {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}
}
